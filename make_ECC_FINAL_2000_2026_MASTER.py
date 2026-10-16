
from pathlib import Path
from datetime import datetime
import shutil
import re
import pandas as pd
import numpy as np

# ============================================================
# ENTEROBACTER CLOACAE COMPLEX
# FINAL MANUSCRIPT DATASET: COLLECTION YEARS 2000-2026
#
# INPUT: Final 9862 country-reported genomes
#
# INCLUDE:
# - Country reported
# - Source classified
# - Collection year 2000 through 2026, inclusive
#
# OUTPUT:
# - One final master workbook in primary master folder
# - One identical workbook in supplementary folder
#
# SUPPLEMENTARY FOLDER:
# Remove all previous content after successfully generating
# and validating the new master.
#
# The original input master is not modified.
# ============================================================

ROOT = Path(r"D:\Enterobacter")
PROJECT = ROOT / "ECC_NPJ_AMR_FIGURES"

MASTER_DIR = (
    PROJECT / "00_FINAL_COUNTRY_REPORTED_MASTER"
)

SUPPLEMENT_DIR = (
    PROJECT / "ECC_SUPPLEMENTARY_AND_SOURCE_DATA"
)

INPUT = (
    MASTER_DIR / "ECC_FINAL_COUNTRY_REPORTED_MASTER.xlsx"
)

OUTPUT_NAME = "ECC_FINAL_2000_2026_MASTER.xlsx"

FINAL_MASTER = MASTER_DIR / OUTPUT_NAME
SUPPLEMENT_MASTER = SUPPLEMENT_DIR / OUTPUT_NAME

MIN_YEAR = 2000
MAX_YEAR = 2026

VALID_SOURCES = {
    "Human",
    "Animal",
    "Environment",
    "Food",
    "Plant-associated",
}

# ============================================================
# 1. READ CURRENT MASTER
# ============================================================

print("=" * 100)
print("ENTEROBACTER CLOACAE COMPLEX")
print("FINAL MASTER DATASET: 2000-2026")
print("=" * 100)

if not INPUT.exists():
    raise FileNotFoundError(INPUT)

df = pd.read_excel(
    INPUT,
    sheet_name="FINAL_COUNTRY_REPORTED",
    engine="openpyxl",
    dtype={"Assembly_Accession": str},
)

print("\nOriginal country-reported genomes:", len(df))
print("Original columns:", len(df.columns))

if len(df) != 9862:
    raise ValueError(
        "Expected 9862 genomes in input master; "
        f"found {len(df)}"
    )

required = [
    "Assembly_Accession",
    "Country_Harmonized",
    "Collection_Year",
    "Final_Source_Category",
]

for col in required:
    if col not in df.columns:
        raise ValueError(
            f"Required column missing: {col}"
        )

df["Assembly_Accession"] = (
    df["Assembly_Accession"]
    .fillna("")
    .astype(str)
    .str.strip()
)

if df["Assembly_Accession"].eq("").any():
    raise ValueError("Empty assembly accession found")

if df["Assembly_Accession"].duplicated().any():
    raise ValueError("Duplicate assembly accession found")

# ============================================================
# 2. VERIFY EXISTING COUNTRY AND SOURCE FILTERS
# ============================================================

print("\n[1] Checking country and One Health source...")

df["Country_Harmonized"] = (
    df["Country_Harmonized"]
    .fillna("")
    .astype(str)
    .str.strip()
)

df["Final_Source_Category"] = (
    df["Final_Source_Category"]
    .fillna("")
    .astype(str)
    .str.strip()
)

missing_country = df["Country_Harmonized"].str.casefold().isin({
    "", "nan", "none", "unknown", "missing",
    "not reported", "not available", "n/a",
})

if missing_country.any():
    raise ValueError(
        f"{missing_country.sum()} missing-country records "
        "were unexpectedly found in the input master"
    )

bad_source = ~df[
    "Final_Source_Category"
].isin(VALID_SOURCES)

if bad_source.any():
    raise ValueError(
        f"{bad_source.sum()} unclassified-source records "
        "were unexpectedly found in the input master"
    )

print("Country verification: PASSED")
print("Source verification: PASSED")

# ============================================================
# 3. PARSE COLLECTION YEAR
# ============================================================

print("\n[2] Validating collection years...")

def parse_year(value):
    if pd.isna(value):
        return np.nan

    text = str(value).strip()

    if text.casefold() in {
        "", "nan", "none", "unknown",
        "missing", "not reported",
        "not available", "n/a",
    }:
        return np.nan

    # Exact four-digit year or Excel-like 2020.0.
    # Do not turn ambiguous year ranges into exact years.
    if re.fullmatch(r"\d{4}(?:\.0+)?", text):
        return int(float(text))

    return np.nan

df["Verified_Collection_Year"] = (
    df["Collection_Year"].map(parse_year)
)

# ============================================================
# 4. APPLY TEMPORAL SELECTION
# ============================================================

print("\n[3] Applying 2000-2026 selection criteria...")

year_valid = (
    df["Verified_Collection_Year"]
    .between(MIN_YEAR, MAX_YEAR)
)

older = df[
    df["Verified_Collection_Year"].notna()
    & (df["Verified_Collection_Year"] < MIN_YEAR)
].copy()

future = df[
    df["Verified_Collection_Year"].notna()
    & (df["Verified_Collection_Year"] > MAX_YEAR)
].copy()

missing_year = df[
    df["Verified_Collection_Year"].isna()
].copy()

final = df[year_valid].copy()

final["Verified_Collection_Year"] = (
    final["Verified_Collection_Year"].astype(int)
)

print("Before 2000:", len(older))
print("Missing/invalid year:", len(missing_year))
print("After 2026:", len(future))
print("Final included:", len(final))

assert (
    len(older)
    + len(missing_year)
    + len(future)
    + len(final)
    == len(df)
)

assert final[
    "Verified_Collection_Year"
].between(2000, 2026).all()

assert not final[
    "Assembly_Accession"
].duplicated().any()

# ============================================================
# 5. MAINTAIN GENOMIC METADATA
# ============================================================

final["Dataset_Status"] = "FINAL_INCLUDED_2000_2026"

final["Temporal_Selection"] = (
    "Reported collection year within 2000-2026"
)

# Preserve original fields, including:
# Species
# AMR genes and classes
# Plasmid replicons
# MLST
# One Health source
# Country and metadata
# Original NCBI accessions

# ============================================================
# 6. YEAR DISTRIBUTION
# ============================================================

year_counts = (
    final.groupby("Verified_Collection_Year")
    .size()
    .reindex(
        range(MIN_YEAR, MAX_YEAR + 1),
        fill_value=0
    )
    .rename("Genome_Count")
    .reset_index()
)

year_counts.columns = [
    "Collection_Year",
    "Genome_Count"
]

# ============================================================
# 7. ONE HEALTH SOURCE DISTRIBUTION
# ============================================================

source_order = [
    "Human",
    "Animal",
    "Environment",
    "Food",
    "Plant-associated"
]

source_counts = (
    final["Final_Source_Category"]
    .value_counts()
    .reindex(source_order, fill_value=0)
    .rename_axis("Source")
    .reset_index(name="Genome_Count")
)

source_counts["Percentage"] = (
    100 * source_counts["Genome_Count"] / len(final)
).round(3)

# ============================================================
# 8. COUNTRY DISTRIBUTION
# ============================================================

country_counts = (
    final.groupby("Country_Harmonized")
    .size()
    .reset_index(name="Genome_Count")
    .sort_values(
        "Genome_Count",
        ascending=False
    )
)

country_counts["Percentage"] = (
    100 * country_counts["Genome_Count"] / len(final)
).round(3)

# ============================================================
# 9. SPECIES DISTRIBUTION
# ============================================================

if "Species" in final.columns:
    species_counts = (
        final["Species"]
        .fillna("Not reported")
        .astype(str)
        .value_counts()
        .rename_axis("Species")
        .reset_index(name="Genome_Count")
    )
else:
    species_counts = pd.DataFrame(
        columns=["Species", "Genome_Count"]
    )

# ============================================================
# 10. YEAR BY SOURCE
# ============================================================

year_source = (
    final.groupby([
        "Verified_Collection_Year",
        "Final_Source_Category"
    ])
    .size()
    .reset_index(name="Genome_Count")
)

# ============================================================
# 11. YEAR BY SPECIES
# ============================================================

if "Species" in final.columns:

    year_species = (
        final.groupby([
            "Verified_Collection_Year",
            "Species"
        ], dropna=False)
        .size()
        .reset_index(name="Genome_Count")
    )

else:
    year_species = pd.DataFrame(
        columns=[
            "Verified_Collection_Year",
            "Species",
            "Genome_Count"
        ]
    )

# ============================================================
# 12. MASTER DATASET SUMMARY
# ============================================================

selection = pd.DataFrame([
    {
        "Stage": "Original high-quality nonredundant ECC",
        "Genome_Count": 10713
    },
    {
        "Stage": "Excluded unresolved source",
        "Genome_Count": 653
    },
    {
        "Stage": "Excluded missing country",
        "Genome_Count": 198
    },
    {
        "Stage": "Country-reported master",
        "Genome_Count": 9862
    },
    {
        "Stage": "Excluded collection year before 2000",
        "Genome_Count": len(older)
    },
    {
        "Stage": "Excluded missing/invalid year",
        "Genome_Count": len(missing_year)
    },
    {
        "Stage": "Excluded collection year after 2026",
        "Genome_Count": len(future)
    },
    {
        "Stage": "FINAL 2000-2026 MASTER",
        "Genome_Count": len(final)
    },
])

metadata = pd.DataFrame([
    ("Dataset", "ECC final 2000-2026"),
    ("Original assembly count", 10713),
    ("Country-reported input", 9862),
    ("Collection year minimum", MIN_YEAR),
    ("Collection year maximum", MAX_YEAR),
    ("Excluded before 2000", len(older)),
    ("Excluded missing/invalid year", len(missing_year)),
    ("Excluded after 2026", len(future)),
    ("Final included genomes", len(final)),
    ("Number of countries", len(country_counts)),
    ("Number of reported species labels", len(species_counts)),
    ("Source classification", "Final manually curated"),
    ("Created", datetime.now().isoformat()),
], columns=["Parameter", "Value"])

# ============================================================
# 13. SAVE PRIMARY FINAL MASTER
# ============================================================

print("\n[4] Writing new final master workbook...")

MASTER_DIR.mkdir(parents=True, exist_ok=True)

with pd.ExcelWriter(
    FINAL_MASTER,
    engine="openpyxl"
) as writer:

    final.to_excel(
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

    species_counts.to_excel(
        writer,
        sheet_name="SPECIES_COUNTS",
        index=False
    )

    year_source.to_excel(
        writer,
        sheet_name="YEAR_SOURCE_COUNTS",
        index=False
    )

    year_species.to_excel(
        writer,
        sheet_name="YEAR_SPECIES_COUNTS",
        index=False
    )

# ============================================================
# 14. VERIFY WRITTEN MASTER
# ============================================================

print("\n[5] Verifying saved final master...")

check = pd.read_excel(
    FINAL_MASTER,
    sheet_name="FINAL_MASTER_2000_2026",
    usecols=[
        "Assembly_Accession",
        "Verified_Collection_Year"
    ],
    engine="openpyxl"
)

assert len(check) == len(final)

assert check[
    "Verified_Collection_Year"
].between(2000, 2026).all()

assert not check[
    "Assembly_Accession"
].duplicated().any()

print("Saved master validation: PASSED")

# ============================================================
# 15. CLEAN SUPPLEMENTARY DIRECTORY
#
# This step is intentionally performed LAST,
# after the primary final master is successfully saved.
#
# The folder will contain ONLY the current master.
# ============================================================

print("\n[6] Updating supplementary directory...")

if not SUPPLEMENT_DIR.exists():
    SUPPLEMENT_DIR.mkdir(parents=True)

# Safety checks to ensure the correct directory is cleaned.
if SUPPLEMENT_DIR.name != "ECC_SUPPLEMENTARY_AND_SOURCE_DATA":
    raise RuntimeError("Unexpected supplementary folder name")

if SUPPLEMENT_DIR.parent.resolve() != PROJECT.resolve():
    raise RuntimeError("Unexpected supplementary folder location")

# Check for Windows temporary/locked Excel files.
locked = list(SUPPLEMENT_DIR.rglob("~$*.xlsx"))

if locked:
    raise RuntimeError(
        "Close open Excel files inside supplementary folder "
        "before cleaning:\n"
        + "\n".join(str(p) for p in locked)
    )

# Copy the verified new master first to a staging location
# outside the directory that will be cleaned.
staging = MASTER_DIR / "_ECC_SUPPLEMENT_STAGING.xlsx"

shutil.copy2(FINAL_MASTER, staging)

if not staging.exists():
    raise RuntimeError("Could not stage final master")

# Delete existing contents of supplementary folder.
for item in SUPPLEMENT_DIR.iterdir():

    if item.is_dir():
        shutil.rmtree(item)

    elif item.is_file() or item.is_symlink():
        item.unlink()

# Put back ONLY the verified master.
shutil.copy2(staging, SUPPLEMENT_MASTER)

staging.unlink()

# ============================================================
# 16. FINAL SUPPLEMENTARY VALIDATION
# ============================================================

contents = list(SUPPLEMENT_DIR.iterdir())

if len(contents) != 1:
    raise RuntimeError(
        "Supplementary folder has unexpected extra files"
    )

if contents[0].name != OUTPUT_NAME:
    raise RuntimeError(
        "Supplementary file name mismatch"
    )

if not SUPPLEMENT_MASTER.exists():
    raise RuntimeError("Supplementary master missing")

if (
    SUPPLEMENT_MASTER.stat().st_size
    != FINAL_MASTER.stat().st_size
):
    raise RuntimeError(
        "Primary and supplementary master sizes differ"
    )

# ============================================================
# 17. PRINT RESULTS
# ============================================================

print("\n" + "=" * 100)
print("FINAL ECC DATASET - 2000 TO 2026")
print("=" * 100)

print("\nGENOME SELECTION")
print(selection.to_string(index=False))

print("\nYEAR COVERAGE")
print("Minimum year:", int(final["Verified_Collection_Year"].min()))
print("Maximum year:", int(final["Verified_Collection_Year"].max()))
print("Genomes retained:", len(final))

print("\nSOURCE COUNTS")
print(source_counts.to_string(index=False))

print("\nTOP 20 COUNTRIES")
print(country_counts.head(20).to_string(index=False))

print("\nSPECIES COUNTS")
print(species_counts.head(20).to_string(index=False))

print("\nFINAL PRIMARY MASTER:")
print(FINAL_MASTER)

print("\nSUPPLEMENTARY MASTER:")
print(SUPPLEMENT_MASTER)

print("\nSUPPLEMENTARY FOLDER CONTENTS:")
for p in SUPPLEMENT_DIR.iterdir():
    print(" -", p.name)

print("\nALL VALIDATION CHECKS PASSED")
print("=" * 100)
