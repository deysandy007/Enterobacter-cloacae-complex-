
from pathlib import Path
from datetime import datetime
import pandas as pd

ROOT = Path(r"D:\Enterobacter")
OUT = ROOT / "ECC_NPJ_AMR_FIGURES" / "00_FINAL_COUNTRY_REPORTED_MASTER"
OUT.mkdir(parents=True, exist_ok=True)

MASTER = ROOT / "ECC_FINAL_HIGH_QUALITY_NONREDUNDANT_MASTER.xlsx"

SOURCE = (
    ROOT / "ECC_NPJ_AMR_FIGURES"
    / "00_Source_Metadata_Audit"
    / "10_FINAL_SOURCE_CLASSIFICATION"
    / "ECC_SOURCE_CLASSIFICATION_FINAL_CLEANED.xlsx"
)

VALID_SOURCES = [
    "Human", "Animal", "Environment",
    "Food", "Plant-associated"
]

MISSING_COUNTRY = {
    "", "na", "n/a", "nan", "none", "null",
    "unknown", "missing", "not available",
    "not reported", "not provided", "unspecified",
    "unavailable", "restricted access",
    "not determined", "not collected", "-"
}

ALIASES = {
    "USA": "United States",
    "United States of America": "United States",
    "Czech Republic": "Czechia",
    "Viet Nam": "Vietnam",
    "Russian Federation": "Russia",
    "Türkiye": "Turkey",
    "Republic of Korea": "South Korea"
}

def clean_country(value):
    if pd.isna(value):
        return ""
    text = " ".join(str(value).strip().split())
    if text.casefold() in MISSING_COUNTRY:
        return ""
    return text

def accession(series):
    return series.fillna("").astype(str).str.strip()

print("=" * 90)
print("ECC FINAL COUNTRY-REPORTED MASTER - CORRECTED BUILD")
print("=" * 90)

if not MASTER.exists():
    raise FileNotFoundError(MASTER)

if not SOURCE.exists():
    raise FileNotFoundError(SOURCE)

genomes = pd.read_excel(
    MASTER,
    sheet_name="FINAL_SELECTED_ECC",
    dtype={"Assembly_Accession": str}
)

sources = pd.read_excel(
    SOURCE,
    sheet_name="ALL_10713_FINAL",
    dtype={"Assembly_Accession": str}
)

for df in [genomes, sources]:
    df["Assembly_Accession"] = accession(
        df["Assembly_Accession"]
    )

    if df["Assembly_Accession"].eq("").any():
        raise ValueError("Missing accession detected")

    if df["Assembly_Accession"].duplicated().any():
        raise ValueError("Duplicate accession detected")

if len(genomes) != 10713 or len(sources) != 10713:
    raise ValueError("Input dataset size is incorrect")

if set(genomes["Assembly_Accession"]) != set(
    sources["Assembly_Accession"]
):
    raise ValueError("Accession sets do not match")

print("\nOriginal genomes:", len(genomes))
print("Source records:", len(sources))
print("Accession matching: PASSED")

# ============================================================
# MERGE ALL SOURCE ANNOTATIONS
# ============================================================

source_cols = [
    "Assembly_Accession",
    "Country",
    "Collection_Year",
    "Raw_Host",
    "Raw_Isolation_Source",
    "Normalized_Host",
    "Normalized_Isolation_Source",
    "Final_Source_Category",
    "Original_Manual_Category",
    "Classification_Reason",
    "Classification_Confidence",
    "Manual_Review_Required"
]

source_cols = [c for c in source_cols if c in sources.columns]

source_part = sources[source_cols].copy()

source_part = source_part.rename(columns={
    "Country": "Source_Audit_Country",
    "Collection_Year": "Source_Audit_Year"
})

for c in source_part.columns:
    if c == "Assembly_Accession":
        continue

    if c in genomes.columns:
        source_part = source_part.rename(
            columns={c: "Curated_" + c}
        )

df = genomes.merge(
    source_part,
    on="Assembly_Accession",
    how="left",
    validate="one_to_one",
    indicator=True
)

if not df["_merge"].eq("both").all():
    raise ValueError("Merge incomplete")

df = df.drop(columns="_merge")

df["Final_Source_Category"] = (
    df["Final_Source_Category"].fillna("")
    .astype(str).str.strip()
)

# ============================================================
# COMPARE ORIGINAL COUNTRY AGAINST SOURCE-AUDIT COUNTRY
# ============================================================

df["Country_Original_Clean"] = df["Country"].map(
    clean_country
)

df["Country_Audit_Clean"] = df[
    "Source_Audit_Country"
].map(clean_country)

df["Original_Country_Available"] = (
    df["Country_Original_Clean"].ne("")
)

df["Audit_Country_Available"] = (
    df["Country_Audit_Clean"].ne("")
)

df["Country_Comparison_Status"] = "MATCH"

different = (
    df["Country_Original_Clean"].str.casefold()
    != df["Country_Audit_Clean"].str.casefold()
)

df.loc[different, "Country_Comparison_Status"] = (
    "COUNTRY_VALUES_DIFFER"
)

df.loc[
    df["Original_Country_Available"]
    & ~df["Audit_Country_Available"],
    "Country_Comparison_Status"
] = "ORIGINAL_REPORTED_AUDIT_MISSING"

df.loc[
    ~df["Original_Country_Available"]
    & df["Audit_Country_Available"],
    "Country_Comparison_Status"
] = "ORIGINAL_MISSING_AUDIT_REPORTED"

df["Source_Eligible"] = df[
    "Final_Source_Category"
].isin(VALID_SOURCES)

country_discrepancies = df[
    df["Source_Eligible"]
    & (
        df["Country_Comparison_Status"] != "MATCH"
    )
].copy()

country_discrepancies.to_excel(
    OUT / "ECC_COUNTRY_METADATA_DISCREPANCIES.xlsx",
    index=False
)

# ============================================================
# FINAL EXCLUSION CRITERIA
# ============================================================

df["Exclusion_Reason"] = ""

df.loc[
    ~df["Source_Eligible"],
    "Exclusion_Reason"
] = "Unclassified or unspecified isolation source"

df.loc[
    df["Source_Eligible"]
    & ~df["Original_Country_Available"],
    "Exclusion_Reason"
] = "Country not reported in original genomic master"

included = df[
    df["Source_Eligible"]
    & df["Original_Country_Available"]
].copy()

excluded = df[
    df["Exclusion_Reason"].ne("")
].copy()

# Preserve original country label
included["Country_Harmonized"] = (
    included["Country_Original_Clean"]
    .replace(ALIASES)
)

included["Dataset_Status"] = "INCLUDED"

# ============================================================
# RECONCILIATION
# ============================================================

excluded_source = int(
    (~df["Source_Eligible"]).sum()
)

excluded_country = int(
    (
        df["Source_Eligible"]
        & ~df["Original_Country_Available"]
    ).sum()
)

original_reported_audit_missing = int(
    (
        df["Source_Eligible"]
        & df["Original_Country_Available"]
        & ~df["Audit_Country_Available"]
    ).sum()
)

original_missing_audit_reported = int(
    (
        df["Source_Eligible"]
        & ~df["Original_Country_Available"]
        & df["Audit_Country_Available"]
    ).sum()
)

# Geographic-location data are left unchanged.
# Antarctica and other special geographic labels remain
# reported values, pending income-classification review.

flow = pd.DataFrame([
    {
        "Stage": "Original ECC master",
        "Remaining": len(df),
        "Excluded": 0
    },
    {
        "Stage": "After source classification",
        "Remaining": int(df["Source_Eligible"].sum()),
        "Excluded": excluded_source
    },
    {
        "Stage": "After missing-country exclusion",
        "Remaining": len(included),
        "Excluded": excluded_country
    }
])

source_summary = (
    included["Final_Source_Category"]
    .value_counts()
    .reindex(VALID_SOURCES, fill_value=0)
    .rename_axis("Source")
    .reset_index(name="Genomes")
)

source_summary["Percent"] = (
    source_summary["Genomes"]
    / len(included) * 100
).round(3)

country_counts = (
    included.groupby("Country_Harmonized")
    .size()
    .reset_index(name="Genomes")
    .sort_values("Genomes", ascending=False)
)

country_source = (
    included.groupby([
        "Country_Harmonized",
        "Final_Source_Category"
    ])
    .size()
    .reset_index(name="Genomes")
)

included["Collection_Year_Numeric"] = pd.to_numeric(
    included["Collection_Year"],
    errors="coerce"
)

year_source = (
    included[
        included["Collection_Year_Numeric"]
        .between(1900, 2026)
    ]
    .groupby([
        "Collection_Year_Numeric",
        "Final_Source_Category"
    ])
    .size()
    .reset_index(name="Genomes")
)

# ============================================================
# VALIDATION
# ============================================================

assert len(df) == 10713
assert len(included) + len(excluded) == 10713
assert not included["Assembly_Accession"].duplicated().any()
assert included["Country_Original_Clean"].ne("").all()
assert included["Final_Source_Category"].isin(VALID_SOURCES).all()
assert source_summary["Genomes"].sum() == len(included)
assert country_counts["Genomes"].sum() == len(included)

# ============================================================
# SAVE FINAL MASTER
# ============================================================

final_xlsx = OUT / "ECC_FINAL_COUNTRY_REPORTED_MASTER.xlsx"
final_tsv = OUT / "ECC_FINAL_COUNTRY_REPORTED_MASTER.tsv"
audit_xlsx = OUT / "ECC_COUNTRY_EXCLUSION_AUDIT.xlsx"

included.to_csv(
    final_tsv,
    sep="\t",
    index=False,
    encoding="utf-8-sig"
)

with pd.ExcelWriter(
    final_xlsx, engine="openpyxl"
) as writer:

    included.to_excel(
        writer,
        sheet_name="FINAL_COUNTRY_REPORTED",
        index=False
    )

    source_summary.to_excel(
        writer,
        sheet_name="SOURCE_COUNTS",
        index=False
    )

    country_counts.to_excel(
        writer,
        sheet_name="COUNTRY_COUNTS",
        index=False
    )

    country_source.to_excel(
        writer,
        sheet_name="COUNTRY_SOURCE_COUNTS",
        index=False
    )

    year_source.to_excel(
        writer,
        sheet_name="YEAR_SOURCE_COUNTS",
        index=False
    )

    flow.to_excel(
        writer,
        sheet_name="EXCLUSION_FLOW",
        index=False
    )

with pd.ExcelWriter(
    audit_xlsx, engine="openpyxl"
) as writer:

    excluded.to_excel(
        writer,
        sheet_name="ALL_EXCLUDED",
        index=False
    )

    country_discrepancies.to_excel(
        writer,
        sheet_name="COUNTRY_DISCREPANCIES",
        index=False
    )

    flow.to_excel(
        writer,
        sheet_name="EXCLUSION_FLOW",
        index=False
    )

# ============================================================
# REPORT
# ============================================================

report = [
    "ECC COUNTRY-REPORTED DATASET",
    "=" * 75,
    f"Original genomes: {len(df)}",
    f"Source excluded: {excluded_source}",
    f"Country excluded after source filtering: {excluded_country}",
    f"Final included: {len(included)}",
    f"Overall retention: {len(included) / len(df) * 100:.3f}%",
    "",
    "COUNTRY AUDIT DISCREPANCIES",
    f"Original reported; audit missing: "
    f"{original_reported_audit_missing}",
    f"Original missing; audit reported: "
    f"{original_missing_audit_reported}",
    f"All differing country entries: "
    f"{len(country_discrepancies)}",
    "",
    "SOURCE DISTRIBUTION",
    source_summary.to_string(index=False),
    "",
    "TOP COUNTRIES",
    country_counts.head(20).to_string(index=False),
    "",
    "NOTE:",
    "Original master country column is used for filtering.",
    "Differences with source-audit metadata are retained",
    "in ECC_COUNTRY_METADATA_DISCREPANCIES.xlsx.",
    "An income group has NOT been assigned at this stage.",
    "",
    f"Generated: {datetime.now().isoformat()}"
]

report_path = OUT / "ECC_COUNTRY_SELECTION_REPORT.txt"

report_path.write_text(
    "\n".join(report),
    encoding="utf-8"
)

print("\n" + "=" * 90)
print("FINAL DATASET VALIDATION")
print("=" * 90)

print("\n".join(report))

print("\nOUTPUTS")
print(final_xlsx)
print(final_tsv)
print(audit_xlsx)
print(
    OUT / "ECC_COUNTRY_METADATA_DISCREPANCIES.xlsx"
)
print(report_path)

print("\nALL VALIDATION CHECKS PASSED")
