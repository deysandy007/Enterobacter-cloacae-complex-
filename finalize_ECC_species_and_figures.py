
from pathlib import Path
from datetime import datetime
import re
import shutil
import subprocess
import sys
import tempfile
import hashlib

import pandas as pd

# ============================================================
# ECC FINAL SPECIES NORMALIZATION AND FIGURE REBUILD
#
# 1. Validate final 2000-2026 master
# 2. Preserve original species annotation
# 3. Normalize to species level
# 4. Add taxonomy and exclusion audit worksheets
# 5. Replace current final master with normalized version
# 6. Rebuild Figures 1A and 1B
# 7. Update supplementary folder with master only
#
# No genome is removed during taxonomic normalization.
# ============================================================

ROOT = Path(r"D:\Enterobacter")
PROJECT = ROOT / "ECC_NPJ_AMR_FIGURES"

MASTER_DIR = (
    PROJECT / "00_FINAL_COUNTRY_REPORTED_MASTER"
)

MASTER = (
    MASTER_DIR / "ECC_FINAL_2000_2026_MASTER.xlsx"
)

BACKUP = (
    MASTER_DIR
    / "ECC_FINAL_2000_2026_BEFORE_TAXONOMY.xlsx"
)

REBUILD = (
    MASTER_DIR / "rebuild_FINAL_FIGURES_1A_1B.py"
)

SUPPLEMENT = (
    PROJECT / "ECC_SUPPLEMENTARY_AND_SOURCE_DATA"
)

SOURCE_COUNTRY_MASTER = (
    MASTER_DIR / "ECC_FINAL_COUNTRY_REPORTED_MASTER.xlsx"
)

FIG_A = (
    PROJECT / "Fig_01A_Global_Geographical_Distribution"
)

FIG_B = (
    PROJECT / "Fig_01B_Temporal_Distribution"
)

EXPECTED_FIGURES = [
    FIG_A / "Fig_01A_Global_Geographical_Distribution.jpg",
    FIG_A / "Fig_01A_Country_Income_Database.xlsx",
    FIG_B / "Fig_01B_Temporal_Distribution.jpg",
    FIG_B / "Fig_01B_Temporal_Data.xlsx",
]

VALID_SOURCES = {
    "Human",
    "Animal",
    "Environment",
    "Food",
    "Plant-associated"
}

# ============================================================
# NORMALIZATION FUNCTION
# ============================================================

def normalize_species(value):

    if pd.isna(value):
        raw = ""
    else:
        raw = " ".join(str(value).strip().split())

    if not raw or raw.casefold() in {
        "nan", "none", "unknown",
        "unclassified", "not reported"
    }:
        return (
            "Enterobacter cloacae complex - unresolved",
            "Enterobacter",
            "",
            "Unresolved",
            "Missing or unresolved taxonomy"
        )

    s = raw.casefold()

    # --------------------------------------------------------
    # Unresolved species-complex labels
    # --------------------------------------------------------

    if (
        "enterobacter cloacae complex" in s
        or "enterobacter cloacae species complex" in s
    ):
        return (
            "Enterobacter cloacae complex - unresolved",
            "Enterobacter",
            "",
            "Unresolved",
            "Complex-level annotation; no species assigned"
        )

    # --------------------------------------------------------
    # Named Enterobacter species
    #
    # Accept taxonomic qualifiers such as subsp. and
    # strain identifiers, but preserve original strings.
    # --------------------------------------------------------

    match = re.match(
        r"^enterobacter\s+([a-z][a-z0-9_-]*)\b",
        s
    )

    if match:

        epithet = match.group(1)

        # Avoid treating genus-level informal labels
        # as established species.
        if epithet in {
            "sp", "spp", "species",
            "cloacae-complex",
            "unclassified"
        }:
            return (
                "Enterobacter cloacae complex - unresolved",
                "Enterobacter",
                "",
                "Unresolved",
                "Genus-level or informal species annotation"
            )

        parent = "Enterobacter " + epithet

        subspecies = ""

        submatch = re.search(
            r"\bsubsp\.?\s+([a-z][a-z0-9_-]*)",
            s
        )

        if submatch:
            subspecies = (
                "Enterobacter "
                + epithet
                + " subsp. "
                + submatch.group(1)
            )

            resolution = "Subspecies-derived"
            rule = "Collapse named subspecies to parent species"

        elif s == parent.casefold():
            resolution = "Species"
            rule = "Original species annotation retained"

        else:
            resolution = "Species-derived"
            rule = (
                "Remove strain/qualifier from species annotation"
            )

        return (
            parent,
            "Enterobacter",
            subspecies,
            resolution,
            rule
        )

    # Unexpected names are not assigned to a guessed species.
    return (
        "Enterobacter cloacae complex - unresolved",
        "Enterobacter",
        "",
        "Unresolved",
        "Unexpected or ambiguous original taxonomy"
    )

# ============================================================
# FILE VALIDATION
# ============================================================

print("=" * 100)
print("ECC FINAL TAXONOMIC NORMALIZATION")
print("=" * 100)

if not MASTER.exists():
    raise FileNotFoundError(
        "2000-2026 master not found:\n" + str(MASTER)
    )

if not REBUILD.exists():
    raise FileNotFoundError(
        "Figure rebuild script not found:\n" + str(REBUILD)
    )

df = pd.read_excel(
    MASTER,
    sheet_name="FINAL_MASTER_2000_2026",
    engine="openpyxl",
    dtype={"Assembly_Accession": str}
)

required = [
    "Assembly_Accession",
    "Species",
    "Country_Harmonized",
    "Final_Source_Category",
    "Verified_Collection_Year"
]

for column in required:
    if column not in df.columns:
        raise ValueError(
            f"Missing required column: {column}"
        )

if df.empty:
    raise ValueError("Input master contains no genomes")

if df["Assembly_Accession"].duplicated().any():
    raise ValueError("Duplicate assembly accession found")

years = pd.to_numeric(
    df["Verified_Collection_Year"],
    errors="coerce"
)

if years.isna().any():
    raise ValueError("Missing collection years in final master")

if not years.between(2000, 2026).all():
    raise ValueError("Out-of-range collection year found")

if not (years % 1 == 0).all():
    raise ValueError("Noninteger collection year found")

if not df["Final_Source_Category"].isin(
    VALID_SOURCES
).all():
    raise ValueError(
        "Unexpected source classification in final master"
    )

N = len(df)

print("\nFinal temporal master genomes:", N)
print("Years:", int(years.min()), "-", int(years.max()))

# ============================================================
# PRESERVE ORIGINAL ANNOTATIONS
# ============================================================

if "Original_Species" in df.columns:
    original = df["Original_Species"].copy()
else:
    original = df["Species"].copy()

df["Original_Species"] = original

# ============================================================
# NORMALIZE ALL TAXONOMIC LABELS
# ============================================================

print("\nNormalizing species names...")

normalized_rows = [
    normalize_species(x)
    for x in df["Original_Species"]
]

normalized = pd.DataFrame(
    normalized_rows,
    columns=[
        "Normalized_Species",
        "Normalized_Genus",
        "Original_Subspecies",
        "Taxonomic_Resolution",
        "Taxonomic_Normalization_Rule"
    ],
    index=df.index
)

for col in normalized.columns:
    df[col] = normalized[col]

# The existing figure-rebuild script reads "Species".
# Set that plotting field to normalized species while
# retaining the complete original name separately.
df["Species"] = df["Normalized_Species"]

df["Taxonomic_Normalization_Method"] = (
    "Rule-based normalization of original NCBI "
    "taxonomic metadata; no genomic reidentification"
)

# ============================================================
# VALIDATION OF GENOME COUNTS
# ============================================================

assert len(df) == N

assert df["Normalized_Species"].notna().all()

assert df["Species"].equals(df["Normalized_Species"])

assert not df["Assembly_Accession"].duplicated().any()

assert df["Normalized_Genus"].eq("Enterobacter").all()

# ============================================================
# TAXONOMIC SUMMARIES
# ============================================================

species_summary = (
    df["Normalized_Species"]
    .value_counts()
    .rename_axis("Normalized_Species")
    .reset_index(name="Genome_Count")
)

species_summary["Percentage"] = (
    species_summary["Genome_Count"] / N * 100
).round(3)

raw_summary = (
    df.groupby(
        ["Original_Species", "Normalized_Species"],
        dropna=False
    )
    .size()
    .reset_index(name="Genome_Count")
    .sort_values(
        ["Normalized_Species", "Genome_Count"],
        ascending=[True, False]
    )
)

taxonomy_rules = (
    df.groupby(
        [
            "Original_Species",
            "Normalized_Species",
            "Taxonomic_Resolution",
            "Taxonomic_Normalization_Rule"
        ],
        dropna=False
    )
    .size()
    .reset_index(name="Genome_Count")
)

# ============================================================
# YEAR/SPECIES SUMMARY
# ============================================================

df["Verified_Collection_Year"] = years.astype(int)

year_species = (
    df.groupby(
        [
            "Verified_Collection_Year",
            "Normalized_Species"
        ]
    )
    .size()
    .reset_index(name="Genome_Count")
)

year_counts = (
    df.groupby("Verified_Collection_Year")
    .size()
    .reindex(
        range(2000, 2027),
        fill_value=0
    )
    .rename("Genome_Count")
    .reset_index()
)

country_counts = (
    df.groupby("Country_Harmonized")
    .size()
    .reset_index(name="Genome_Count")
    .sort_values(
        "Genome_Count",
        ascending=False
    )
)

source_counts = (
    df.groupby("Final_Source_Category")
    .size()
    .reset_index(name="Genome_Count")
)

# ============================================================
# GENOME SELECTION AND EXCLUSION
# ============================================================

original_total = 10713
source_excluded = 653
country_excluded = 198
country_reported = 9862

if SOURCE_COUNTRY_MASTER.exists():
    original_country = pd.read_excel(
        SOURCE_COUNTRY_MASTER,
        sheet_name="FINAL_COUNTRY_REPORTED",
        usecols=[
            "Assembly_Accession",
            "Collection_Year"
        ],
        dtype={"Assembly_Accession": str}
    )

    if len(original_country) != country_reported:
        raise ValueError(
            "Country master size differs from 9862"
        )

    year_values = pd.to_numeric(
        original_country["Collection_Year"],
        errors="coerce"
    )

    strict_year = (
        original_country["Collection_Year"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.fullmatch(r"\d{4}(?:\.0+)?")
    )

    valid_exact = strict_year & year_values.notna()

    before_2000 = int(
        (valid_exact & (year_values < 2000)).sum()
    )

    after_2026 = int(
        (valid_exact & (year_values > 2026)).sum()
    )

    invalid_year = int((~valid_exact).sum())

    expected_final = (
        country_reported
        - before_2000
        - after_2026
        - invalid_year
    )

    if expected_final != N:
        raise ValueError(
            "Year exclusion accounting differs from "
            f"final master: {expected_final} versus {N}. "
            "Review date filtering before proceeding."
        )

else:
    raise FileNotFoundError(
        "Original 9862-genome country master is required "
        "to audit year exclusions."
    )

selection = pd.DataFrame([
    {
        "Stage": "Initial NCBI Enterobacter retrieval preview",
        "Genomes": 33805,
        "Decision": "Broad preliminary retrieval",
        "Caution": "Not a final QC denominator"
    },
    {
        "Stage": "ECC genome assembly collection",
        "Genomes": 31412,
        "Decision": "Selected ECC assemblies",
        "Caution": "Input to QC/deduplication"
    },
    {
        "Stage": "Quality-filtered nonredundant ECC master",
        "Genomes": 10713,
        "Decision": "Combined assembly QC and deduplication",
        "Caution": "Separate QC and dedup counts require logs"
    },
    {
        "Stage": "Source-unclassified excluded",
        "Genomes": source_excluded,
        "Decision": "Unclassified or unspecified source",
        "Caution": "Exclusion"
    },
    {
        "Stage": "Source-classified retained",
        "Genomes": 10060,
        "Decision": "Five One Health categories",
        "Caution": "Intermediate included dataset"
    },
    {
        "Stage": "Missing-country excluded",
        "Genomes": country_excluded,
        "Decision": "Country absent in original genome master",
        "Caution": "Exclusion"
    },
    {
        "Stage": "Country-reported retained",
        "Genomes": country_reported,
        "Decision": "Reported country and classified source",
        "Caution": "Intermediate included dataset"
    },
    {
        "Stage": "Collection year before 2000 excluded",
        "Genomes": before_2000,
        "Decision": "Outside temporal study period",
        "Caution": "Exclusion"
    },
    {
        "Stage": "Collection year after 2026 excluded",
        "Genomes": after_2026,
        "Decision": "Outside temporal study period",
        "Caution": "Exclusion"
    },
    {
        "Stage": "Missing/ambiguous/invalid year excluded",
        "Genomes": invalid_year,
        "Decision": "No valid four-digit collection year",
        "Caution": "Exclusion"
    },
    {
        "Stage": "FINAL ECC MASTER 2000-2026",
        "Genomes": N,
        "Decision": "Species-normalized final analysis dataset",
        "Caution": "Final included dataset"
    }
])

# ============================================================
# MANUSCRIPT QC METADATA
# ============================================================

metadata = pd.DataFrame([
    ("Genomes in final master", N),
    ("Collection year minimum", int(years.min())),
    ("Collection year maximum", int(years.max())),
    ("Distinct country labels", len(country_counts)),
    ("Distinct normalized species", len(species_summary)),
    ("Original NCBI taxonomy retained", "Yes"),
    ("Normalized genus", "Enterobacter"),
    ("Taxonomic assignment method",
     "Rule-based standardization of NCBI species labels"),
    ("Subspecies treatment",
     "Collapsed into parent species"),
    ("Unresolved complex labels",
     "Preserved in unresolved category"),
    ("Source categories", ", ".join(sorted(VALID_SOURCES))),
    ("Created", datetime.now().isoformat()),
], columns=["Parameter", "Value"])

# ============================================================
# PRESERVE OTHER EXISTING WORKSHEETS
# ============================================================

existing = pd.ExcelFile(MASTER)

preserved_sheets = {}

replaced = {
    "FINAL_MASTER_2000_2026",
    "DATASET_SUMMARY",
    "SELECTION_FLOW",
    "YEAR_COUNTS",
    "SOURCE_COUNTS",
    "COUNTRY_COUNTS",
    "SPECIES_COUNTS",
    "YEAR_SPECIES_COUNTS",
    "RAW_TAXONOMY_AUDIT",
    "TAXONOMY_RULES"
}

for sheet in existing.sheet_names:
    if sheet not in replaced:
        preserved_sheets[sheet] = pd.read_excel(
            MASTER, sheet_name=sheet
        )

# ============================================================
# BACKUP BEFORE REPLACEMENT
# ============================================================

if not BACKUP.exists():
    shutil.copy2(MASTER, BACKUP)

print("\nBackup created/preserved:")
print(BACKUP)

# ============================================================
# WRITE TEMPORARY NORMALIZED MASTER
# ============================================================

with tempfile.TemporaryDirectory(
    prefix="ECC_NORMALIZATION_",
    dir=str(MASTER_DIR)
) as tmp:

    temp_master = Path(tmp) / MASTER.name

    with pd.ExcelWriter(
        temp_master,
        engine="openpyxl"
    ) as writer:

        df.to_excel(
            writer,
            sheet_name="FINAL_MASTER_2000_2026",
            index=False
        )

        metadata.to_excel(
            writer,
            sheet_name="DATASET_SUMMARY",
            index=False
        )

        selection.to_excel(
            writer,
            sheet_name="SELECTION_FLOW",
            index=False
        )

        year_counts.to_excel(
            writer,
            sheet_name="YEAR_COUNTS",
            index=False
        )

        source_counts.to_excel(
            writer,
            sheet_name="SOURCE_COUNTS",
            index=False
        )

        country_counts.to_excel(
            writer,
            sheet_name="COUNTRY_COUNTS",
            index=False
        )

        species_summary.to_excel(
            writer,
            sheet_name="SPECIES_COUNTS",
            index=False
        )

        year_species.to_excel(
            writer,
            sheet_name="YEAR_SPECIES_COUNTS",
            index=False
        )

        raw_summary.to_excel(
            writer,
            sheet_name="RAW_TAXONOMY_AUDIT",
            index=False
        )

        taxonomy_rules.to_excel(
            writer,
            sheet_name="TAXONOMY_RULES",
            index=False
        )

        for sheet, data in preserved_sheets.items():
            data.to_excel(
                writer,
                sheet_name=sheet,
                index=False
            )

    check = pd.read_excel(
        temp_master,
        sheet_name="FINAL_MASTER_2000_2026",
        usecols=[
            "Assembly_Accession",
            "Original_Species",
            "Normalized_Species",
            "Verified_Collection_Year"
        ]
    )

    assert len(check) == N
    assert check["Normalized_Species"].notna().all()
    assert check["Assembly_Accession"].is_unique

    shutil.copy2(temp_master, MASTER)

print("\nFINAL SPECIES-NORMALIZED MASTER UPDATED")
print(MASTER)

# ============================================================
# SPECIES NORMALIZATION SUMMARY
# ============================================================

print("\n" + "=" * 100)
print("FINAL NORMALIZED SPECIES COUNTS")
print("=" * 100)

print(species_summary.to_string(index=False))

print("\nTAXONOMIC RESOLUTION")
print(
    df["Taxonomic_Resolution"]
    .value_counts()
    .to_string()
)

print("\nGENOME EXCLUSION FLOW")
print(selection.to_string(index=False))

# ============================================================
# REGENERATE THE TWO FIGURES
# ============================================================

print("\n" + "=" * 100)
print("REBUILDING FIGURES 1A AND 1B")
print("=" * 100)

result = subprocess.run(
    [sys.executable, str(REBUILD)],
    cwd=str(MASTER_DIR),
    check=False
)

if result.returncode != 0:
    raise RuntimeError(
        "Figure rebuilding failed. "
        "The normalized master and backup are available, "
        "but the supplementary folder was not changed."
    )

for path in EXPECTED_FIGURES:
    if not path.exists():
        raise FileNotFoundError(
            "Expected figure output missing: " + str(path)
        )

# ============================================================
# VERIFY FIGURE DATABASE COUNTS
# ============================================================

fig_a = pd.read_excel(
    EXPECTED_FIGURES[1],
    sheet_name="COUNTRY_INCOME_COUNTS"
)

fig_b = pd.read_excel(
    EXPECTED_FIGURES[3],
    sheet_name="ANNUAL_COUNTS"
)

if int(fig_a["Genome_Count"].sum()) != N:
    raise ValueError("Figure 1A sample count mismatch")

if int(fig_b["Genome_Count"].sum()) != N:
    raise ValueError("Figure 1B sample count mismatch")

# ============================================================
# UPDATE SUPPLEMENTARY DIRECTORY
#
# Only after taxonomy and both figures passed validation.
# ============================================================

print("\nUpdating supplementary master folder...")

SUPPLEMENT.mkdir(parents=True, exist_ok=True)

if (
    SUPPLEMENT.name != "ECC_SUPPLEMENTARY_AND_SOURCE_DATA"
    or SUPPLEMENT.parent.resolve() != PROJECT.resolve()
):
    raise RuntimeError("Unsafe supplementary folder path")

# Ensure Excel is closed before clearing the directory.
locked = list(SUPPLEMENT.rglob("~$*.xlsx"))

if locked:
    raise RuntimeError(
        "Close supplementary Excel workbooks first: "
        + str(locked)
    )

# Stage the new workbook outside the folder to be cleared.
with tempfile.TemporaryDirectory(
    prefix="ECC_SUPPLEMENT_STAGE_",
    dir=str(MASTER_DIR)
) as tmp:

    staged = Path(tmp) / MASTER.name
    shutil.copy2(MASTER, staged)

    for item in SUPPLEMENT.iterdir():
        if item.is_dir() and not item.is_symlink():
            shutil.rmtree(item)
        else:
            item.unlink()

    supplement_master = SUPPLEMENT / MASTER.name
    shutil.copy2(staged, supplement_master)

# ============================================================
# CHECK FILE IDENTITIES
# ============================================================

def file_hash(path):
    digest = hashlib.sha256()

    with open(path, "rb") as f:
        for block in iter(
            lambda: f.read(1048576), b""
        ):
            digest.update(block)

    return digest.hexdigest()

assert file_hash(MASTER) == file_hash(
    SUPPLEMENT / MASTER.name
)

assert len(list(SUPPLEMENT.iterdir())) == 1

# ============================================================
# FINAL REPORT
# ============================================================

print("\n" + "=" * 100)
print("FINAL ECC TAXONOMY AND FIGURES COMPLETE")
print("=" * 100)

print("\nFinal genomes:", N)
print("Collection period: 2000-2026")
print("Normalized species:", len(species_summary))

print("\nFinal primary master:")
print(MASTER)

print("\nSupplementary master:")
print(SUPPLEMENT / MASTER.name)

print("\nFigure 1A:")
print(EXPECTED_FIGURES[0])
print(EXPECTED_FIGURES[1])

print("\nFigure 1B:")
print(EXPECTED_FIGURES[2])
print(EXPECTED_FIGURES[3])

print("\nBoth figures generated and counts verified.")
print("Supplementary folder contains one master only.")
print("=" * 100)
