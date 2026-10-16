
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd

# ============================================================
# FIGURE 2C - COUNTRY-LEVEL ECC SPECIES DIVERSITY AUDIT
#
# INPUT:
# Final normalized 9157-genome master, 2000-2026
# Frozen Figure 1A country-income classifications
#
# OUTPUT:
# ONE Excel workbook only.
#
# No arbitrary country exclusion threshold.
# No changes to master or supplementary folders.
#
# Diversity metrics calculated from resolved species only.
# Unresolved ECC retained separately for transparency.
# ============================================================

ROOT = Path(r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES")

MASTER = (
    ROOT / "00_FINAL_COUNTRY_REPORTED_MASTER"
    / "ECC_FINAL_2000_2026_MASTER.xlsx"
)

INCOME_FILE = (
    ROOT / "Fig_01A_Global_Geographical_Distribution"
    / "Fig_01A_Country_Income_Database.xlsx"
)

OUT = ROOT / "Fig_02C_Country_Species_Diversity"
OUT.mkdir(parents=True, exist_ok=True)

OUTPUT = OUT / "Fig_02C_Country_Diversity_Audit.xlsx"

EXPECTED_N = 9157

UNRESOLVED = (
    "Enterobacter cloacae complex - unresolved"
)

INCOME_ORDER = [
    "H", "M", "L+LM", "Unclassified"
]

THRESHOLDS = [1, 5, 10, 20, 30, 50, 100, 200, 500]

# ============================================================
# 1. LOAD GENOMIC MASTER
# ============================================================

print("=" * 100)
print("FIGURE 2C - COUNTRY-LEVEL SPECIES DIVERSITY AUDIT")
print("=" * 100)

if not MASTER.exists():
    raise FileNotFoundError(MASTER)

if not INCOME_FILE.exists():
    raise FileNotFoundError(INCOME_FILE)

df = pd.read_excel(
    MASTER,
    sheet_name="FINAL_MASTER_2000_2026",
    engine="openpyxl",
    dtype={"Assembly_Accession": str}
)

required = [
    "Assembly_Accession",
    "Country_Harmonized",
    "Normalized_Species",
    "Verified_Collection_Year"
]

for col in required:
    if col not in df.columns:
        raise ValueError("Missing master column: " + col)

if len(df) != EXPECTED_N:
    raise ValueError(
        f"Expected {EXPECTED_N} genomes; found {len(df)}"
    )

if df["Assembly_Accession"].isna().any():
    raise ValueError("Missing assembly accession")

if df["Assembly_Accession"].duplicated().any():
    raise ValueError("Duplicate assembly accession")

years = pd.to_numeric(
    df["Verified_Collection_Year"],
    errors="coerce"
)

if (
    years.isna().any()
    or not years.between(2000, 2026).all()
    or not (years % 1 == 0).all()
):
    raise ValueError("Invalid collection year")

df["Year"] = years.astype(int)

df["Country"] = (
    df["Country_Harmonized"]
    .fillna("")
    .astype(str)
    .str.strip()
)

df["Species"] = (
    df["Normalized_Species"]
    .fillna("")
    .astype(str)
    .str.strip()
)

if df["Country"].eq("").any():
    raise ValueError("Missing country")

if df["Species"].eq("").any():
    raise ValueError("Missing normalized species")

N = len(df)

print("Final genomes:", N)
print("Countries:", df["Country"].nunique())
print("Normalized taxonomy categories:",
      df["Species"].nunique())

# ============================================================
# 2. USE EXACT FIGURE 1A INCOME CLASSIFICATIONS
# ============================================================

income = pd.read_excel(
    INCOME_FILE,
    sheet_name="COUNTRY_INCOME_COUNTS",
    engine="openpyxl"
)

for col in [
    "Country_Harmonized",
    "Income_Group",
    "Genome_Count"
]:
    if col not in income.columns:
        raise ValueError(
            "Missing Figure 1A income column: " + col
        )

income["Country"] = (
    income["Country_Harmonized"]
    .fillna("")
    .astype(str)
    .str.strip()
)

if income["Country"].duplicated().any():
    raise ValueError(
        "Duplicated country in Figure 1A income table"
    )

income["Income_Group"] = (
    income["Income_Group"]
    .fillna("Unclassified")
    .astype(str)
    .str.strip()
)

if not income["Income_Group"].isin(INCOME_ORDER).all():
    raise ValueError("Unexpected income group")

master_counts = df["Country"].value_counts()
income_counts = income.set_index("Country")["Genome_Count"]

check = pd.concat(
    [
        master_counts.rename("Master_Count"),
        income_counts.rename("Figure1A_Count")
    ],
    axis=1
)

if check.isna().any().any():
    raise ValueError(
        "Country mismatch between master and Figure 1A"
    )

if not (
    check["Master_Count"].astype(int)
    == check["Figure1A_Count"].astype(int)
).all():
    raise ValueError(
        "Country counts differ between master and Figure 1A"
    )

df = df.merge(
    income[["Country", "Income_Group"]],
    on="Country",
    how="left",
    validate="many_to_one"
)

if df["Income_Group"].isna().any():
    raise ValueError("Income assignment failed")

assert len(df) == N

# ============================================================
# 3. SPECIES-BY-COUNTRY MATRIX
# ============================================================

matrix = pd.crosstab(
    df["Country"],
    df["Species"]
)

country_order = (
    matrix.sum(axis=1)
    .sort_values(ascending=False)
    .index
)

matrix = matrix.reindex(country_order, fill_value=0)

if UNRESOLVED not in matrix.columns:
    matrix[UNRESOLVED] = 0

resolved_species = [
    col for col in matrix.columns
    if col != UNRESOLVED
]

resolved_matrix = matrix[resolved_species]

if int(matrix.to_numpy().sum()) != N:
    raise RuntimeError("Genome count mismatch")

# ============================================================
# 4. COUNTRY DIVERSITY METRICS
# ============================================================

print("\nCalculating diversity metrics...")

results = []

income_lookup = (
    income.set_index("Country")["Income_Group"]
    .to_dict()
)

for country in matrix.index:

    row = matrix.loc[country]

    all_count = int(row.sum())

    unresolved_count = int(row[UNRESOLVED])

    resolved = (
        resolved_matrix.loc[country]
        .to_numpy(dtype=int)
    )

    resolved_count = int(resolved.sum())

    richness = int((resolved > 0).sum())

    if resolved_count > 0:

        p = resolved[resolved > 0] / resolved_count

        shannon = float(
            -(p * np.log(p)).sum()
        )

        simpson = float(
            1.0 - np.sum(p ** 2)
        )

        inverse_simpson = float(
            1.0 / np.sum(p ** 2)
        )

        pielou = (
            float(shannon / np.log(richness))
            if richness > 1
            else np.nan
        )

        dominant_index = int(np.argmax(resolved))
        dominant_species = resolved_species[dominant_index]

        dominant_percent = (
            100 * resolved[dominant_index]
            / resolved_count
        )

    else:
        shannon = np.nan
        simpson = np.nan
        inverse_simpson = np.nan
        pielou = np.nan
        dominant_species = ""
        dominant_percent = np.nan

    results.append({
        "Country": country,
        "Income_Group": income_lookup[country],
        "Total_Genomes": all_count,
        "Resolved_Genomes": resolved_count,
        "Unresolved_Genomes": unresolved_count,
        "Unresolved_Percent": (
            100 * unresolved_count / all_count
        ),
        "Resolved_Species_Richness": richness,
        "Shannon_H": shannon,
        "Simpson_1_minus_D": simpson,
        "Inverse_Simpson": inverse_simpson,
        "Pielou_Evenness": pielou,
        "Dominant_Species": dominant_species,
        "Dominant_Species_Percent": dominant_percent
    })

country_metrics = pd.DataFrame(results)

country_metrics = country_metrics.sort_values(
    ["Total_Genomes", "Country"],
    ascending=[False, True]
).reset_index(drop=True)

assert (
    country_metrics["Total_Genomes"].sum()
    == N
)

assert (
    country_metrics["Resolved_Genomes"].sum()
    + country_metrics["Unresolved_Genomes"].sum()
    == N
)

# ============================================================
# 5. GENOME-DEPTH ELIGIBILITY
#
# Two denominators are reported:
# A) All available genomes
# B) Resolved genomes only
#
# Only B is eligible for resolved-species rarefaction.
# ============================================================

print("\nAuditing genome-depth thresholds...")

threshold_rows = []

for depth in THRESHOLDS:

    for group in (
        ["ALL_COUNTRIES"] + INCOME_ORDER
    ):

        subset = (
            country_metrics
            if group == "ALL_COUNTRIES"
            else country_metrics[
                country_metrics["Income_Group"] == group
            ]
        )

        total_eligible = subset[
            subset["Total_Genomes"] >= depth
        ]

        resolved_eligible = subset[
            subset["Resolved_Genomes"] >= depth
        ]

        threshold_rows.append({
            "Income_Group": group,
            "Minimum_Genomes": depth,
            "All_Countries": len(subset),
            "Countries_Total_N_Eligible": len(total_eligible),
            "Countries_Resolved_N_Eligible": len(resolved_eligible),
            "Genomes_In_Resolved_Eligible_Countries": int(
                resolved_eligible["Total_Genomes"].sum()
            ),
            "Resolved_Genomes_In_Eligible_Countries": int(
                resolved_eligible["Resolved_Genomes"].sum()
            ),
            "Countries_Excluded_At_Resolved_Depth": (
                len(subset) - len(resolved_eligible)
            )
        })

threshold_audit = pd.DataFrame(threshold_rows)

# ============================================================
# 6. INCOME-GROUP SUMMARY
# ============================================================

income_rows = []

for group in INCOME_ORDER:

    subset = country_metrics[
        country_metrics["Income_Group"] == group
    ]

    if subset.empty:
        continue

    diversity_valid = subset[
        subset["Resolved_Genomes"] > 0
    ]

    income_rows.append({
        "Income_Group": group,
        "Number_of_Countries": len(subset),
        "Total_Genomes": int(subset["Total_Genomes"].sum()),
        "Resolved_Genomes": int(
            subset["Resolved_Genomes"].sum()
        ),
        "Unresolved_Genomes": int(
            subset["Unresolved_Genomes"].sum()
        ),
        "Median_Genomes_Per_Country": float(
            subset["Total_Genomes"].median()
        ),
        "Median_Resolved_Genomes": float(
            subset["Resolved_Genomes"].median()
        ),
        "Median_Observed_Richness": (
            float(diversity_valid[
                "Resolved_Species_Richness"
            ].median())
            if len(diversity_valid) else np.nan
        ),
        "Median_Shannon_H": (
            float(diversity_valid["Shannon_H"].median())
            if len(diversity_valid) else np.nan
        ),
        "Countries_With_No_Resolved_Genomes": int(
            (subset["Resolved_Genomes"] == 0).sum()
        )
    })

income_summary = pd.DataFrame(income_rows)

# ============================================================
# 7. RANKINGS - DESCRIPTIVE ONLY
# ============================================================

richness_rank = country_metrics.sort_values(
    [
        "Resolved_Species_Richness",
        "Resolved_Genomes"
    ],
    ascending=[False, False]
).copy()

shannon_rank = country_metrics[
    country_metrics["Shannon_H"].notna()
].sort_values(
    "Shannon_H",
    ascending=False
).copy()

# These rankings are unadjusted for sampling depth.
# Never interpret them as final diversity comparisons.

# ============================================================
# 8. COUNTRY-SPECIES MATRIX FOR RAREFACTION
# ============================================================

matrix_export = matrix.copy()
matrix_export.insert(
    0, "Total_Genomes", matrix.sum(axis=1)
)
matrix_export.insert(
    1, "Resolved_Genomes", resolved_matrix.sum(axis=1)
)

matrix_export = (
    matrix_export
    .rename_axis("Country")
    .reset_index()
)

matrix_export = matrix_export.merge(
    income[["Country", "Income_Group"]],
    on="Country",
    how="left",
    validate="one_to_one"
)

# ============================================================
# 9. METHODS AND DEFINITIONS
# ============================================================

methods = pd.DataFrame([
    (
        "Dataset",
        "Final 9157 ECC genomes collected during 2000-2026"
    ),
    (
        "Income assignments",
        "Unchanged classifications from final Figure 1A"
    ),
    (
        "Country sampling",
        "All country-reported records retained in audit"
    ),
    (
        "Species richness",
        "Number of distinct resolved species per country"
    ),
    (
        "Shannon H",
        "-sum(p_i * natural_log(p_i))"
    ),
    (
        "Simpson 1-D",
        "1 - sum(p_i squared)"
    ),
    (
        "Inverse Simpson",
        "1 / sum(p_i squared)"
    ),
    (
        "Pielou evenness",
        "Shannon H / log(species richness), when richness > 1"
    ),
    (
        "Diversity denominator",
        "Resolved genomes within each country"
    ),
    (
        "Unresolved ECC genomes",
        "Excluded from species richness and diversity "
        "calculations, but retained in sample counts"
    ),
    (
        "Depth eligibility",
        "Computed separately for all genomes and "
        "resolved genomes"
    ),
    (
        "Rarefaction depth",
        "Not selected - determined after sampling audit"
    ),
    (
        "Observed diversity rankings",
        "Exploratory; not adjusted for sequencing effort"
    ),
    (
        "Created",
        datetime.now().isoformat()
    )
], columns=["Parameter", "Definition"])

# ============================================================
# 10. WRITE ONE AUDIT EXCEL WORKBOOK
# ============================================================

print("\nWriting final country diversity audit workbook...")

with pd.ExcelWriter(
    OUTPUT,
    engine="openpyxl"
) as writer:

    country_metrics.to_excel(
        writer,
        sheet_name="COUNTRY_DIVERSITY",
        index=False
    )

    matrix_export.to_excel(
        writer,
        sheet_name="COUNTRY_SPECIES_MATRIX",
        index=False
    )

    threshold_audit.to_excel(
        writer,
        sheet_name="SAMPLING_DEPTH_AUDIT",
        index=False
    )

    income_summary.to_excel(
        writer,
        sheet_name="INCOME_GROUP_SUMMARY",
        index=False
    )

    richness_rank.to_excel(
        writer,
        sheet_name="OBSERVED_RICHNESS_RANK",
        index=False
    )

    shannon_rank.to_excel(
        writer,
        sheet_name="OBSERVED_SHANNON_RANK",
        index=False
    )

    df[[
        "Assembly_Accession",
        "Country",
        "Income_Group",
        "Year",
        "Species"
    ]].to_excel(
        writer,
        sheet_name="GENOME_SOURCE_DATA",
        index=False
    )

    methods.to_excel(
        writer,
        sheet_name="METHODS_DEFINITIONS",
        index=False
    )

# ============================================================
# 11. VERIFY WORKBOOK
# ============================================================

check_matrix = pd.read_excel(
    OUTPUT,
    sheet_name="COUNTRY_SPECIES_MATRIX"
)

if int(check_matrix["Total_Genomes"].sum()) != N:
    raise RuntimeError("Excel country count validation failed")

# ============================================================
# 12. PRINT SCIENTIFIC AUDIT RESULTS
# ============================================================

print("\n" + "=" * 100)
print("COUNTRY SAMPLING DEPTH - ALL COUNTRIES")
print("=" * 100)

print(
    threshold_audit[
        threshold_audit["Income_Group"] == "ALL_COUNTRIES"
    ][[
        "Minimum_Genomes",
        "Countries_Total_N_Eligible",
        "Countries_Resolved_N_Eligible",
        "Countries_Excluded_At_Resolved_Depth"
    ]].to_string(index=False)
)

print("\n" + "=" * 100)
print("COUNTRY SAMPLING DEPTH BY INCOME")
print("=" * 100)

print(
    threshold_audit[
        threshold_audit["Income_Group"] != "ALL_COUNTRIES"
    ][[
        "Income_Group",
        "Minimum_Genomes",
        "Countries_Resolved_N_Eligible"
    ]].to_string(index=False)
)

print("\n" + "=" * 100)
print("INCOME-GROUP DIVERSITY SUMMARY")
print("=" * 100)

print(income_summary.to_string(index=False))

print("\n" + "=" * 100)
print("TOP 20 COUNTRIES BY AVAILABLE GENOMES")
print("=" * 100)

print(
    country_metrics.head(20)[[
        "Country",
        "Income_Group",
        "Total_Genomes",
        "Resolved_Genomes",
        "Resolved_Species_Richness",
        "Shannon_H",
        "Unresolved_Percent"
    ]].to_string(index=False)
)

print("\n" + "=" * 100)
print("AUDIT COMPLETED")
print("=" * 100)

print("Genomes:", N)
print("Countries:", len(country_metrics))
print("Resolved species:", len(resolved_species))
print("Output:", OUTPUT)
print("No figure generated yet")
print("Final master unchanged")
print("=" * 100)
