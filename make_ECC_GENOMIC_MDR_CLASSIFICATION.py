
#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import re
import pandas as pd
import numpy as np

ROOT = Path(r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES")

MASTER = (
    ROOT / "00_FINAL_COUNTRY_REPORTED_MASTER"
    / "ECC_FINAL_2000_2026_MASTER.xlsx"
)

INCOME_FILE = (
    ROOT / "Fig_01A_Global_Geographical_Distribution"
    / "Fig_01A_Country_Income_Database.xlsx"
)

OUTDIR = ROOT / "Fig_06_ECC_GENOMIC_MDR"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUTPUT = OUTDIR / "ECC_ALL_9157_GENOMIC_MDR_CLASSIFICATION.xlsx"

SPECIES = [
    "Enterobacter hormaechei",
    "Enterobacter cloacae",
    "Enterobacter roggenkampii",
    "Enterobacter asburiae",
    "Enterobacter kobei",
    "Enterobacter ludwigii",
    "Enterobacter bugandensis",
]

# The AMR_Classes column must contain distinct drug-class
# names separated by semicolons.
#
# IMPORTANT:
# This is an annotation-based genomic proxy, not
# phenotypically confirmed MDR.

def parse_classes(value):
    if pd.isna(value):
        return set()

    text = str(value).strip()

    if not text or text.lower() in {
        "nan", "none", "na", "n/a",
        "not detected", "no amr",
        "no amr classes", "-"
    }:
        return set()

    # Preserve class names as supplied by the master.
    # Separate only on semicolons, which avoids splitting
    # composite class names on hyphens or slashes.
    return {
        token.strip()
        for token in text.split(";")
        if token.strip()
    }


def classify(n):
    if n == 0:
        return "0 classes"
    if n == 1:
        return "1 class"
    if n == 2:
        return "2 classes"
    return "Genomic MDR >=3 classes"


def normalize_income(value):
    s = str(value).strip().upper().replace(" ", "")

    if s in {"H", "HIGH", "HIGHINCOME"}:
        return "H"

    if s in {
        "M", "UM", "UPPERMIDDLE",
        "UPPER-MIDDLE", "UPPERMIDDLEINCOME",
        "UPPER-MIDDLEINCOME"
    }:
        return "M"

    if s in {
        "L+LM", "L_LM", "LM+L",
        "LOW+LOWERMIDDLE",
        "LOWERMIDDLE+LOW",
        "LOWERMIDDLEANDLOW"
    }:
        return "L+LM"

    return str(value).strip()


print("=" * 100)
print("ECC — GENOME-LEVEL GENOMIC MDR CLASSIFICATION")
print("=" * 100)

if not MASTER.is_file():
    raise FileNotFoundError(MASTER)

df = pd.read_excel(
    MASTER,
    sheet_name="FINAL_MASTER_2000_2026",
    engine="openpyxl",
    dtype={"Assembly_Accession": str}
)

required = [
    "Assembly_Accession",
    "Normalized_Species",
    "Country_Harmonized",
    "AMR_Genes",
    "AMR_Classes"
]

missing = [c for c in required if c not in df.columns]

if missing:
    raise ValueError(
        "Required columns missing: " + str(missing)
    )

if len(df) != 9157:
    raise ValueError(
        f"Expected 9157 genomes, found {len(df)}"
    )

if (
    df["Assembly_Accession"].isna().any()
    or df["Assembly_Accession"].duplicated().any()
):
    raise ValueError(
        "Missing or duplicate assembly accession detected"
    )

print("Master genomes:", len(df))

# ============================================================
# FROZEN COUNTRY-INCOME MAPPING
# ============================================================

if INCOME_FILE.is_file():

    income = pd.read_excel(
        INCOME_FILE,
        sheet_name="COUNTRY_INCOME_COUNTS",
        engine="openpyxl"
    )

    needed = ["Country_Harmonized", "Income_Group"]

    if not set(needed).issubset(income.columns):
        raise ValueError(
            "Missing income mapping columns"
        )

    income["Country_Key"] = (
        income["Country_Harmonized"]
        .fillna("").astype(str).str.strip()
    )

    income["Income_Group"] = (
        income["Income_Group"]
        .map(normalize_income)
    )

    if income["Country_Key"].duplicated().any():
        raise ValueError(
            "Duplicate countries in income mapping"
        )

    df["Country_Key"] = (
        df["Country_Harmonized"]
        .fillna("").astype(str).str.strip()
    )

    mapping = income.set_index(
        "Country_Key"
    )["Income_Group"]

    df["MDR_Income_Group"] = df["Country_Key"].map(
        mapping
    )

else:
    print("WARNING: Frozen income mapping not found")

    if "Income_Group" in df.columns:
        df["MDR_Income_Group"] = (
            df["Income_Group"].map(normalize_income)
        )
    else:
        df["MDR_Income_Group"] = "Unknown"

df["MDR_Income_Group"] = (
    df["MDR_Income_Group"]
    .fillna("Unknown")
)

# ============================================================
# CALCULATE CLASS COUNTS
# ============================================================

df["MDR_Class_Set"] = (
    df["AMR_Classes"].apply(parse_classes)
)

df["MDR_Class_Count"] = (
    df["MDR_Class_Set"].apply(len)
)

df["MDR_Class_List"] = (
    df["MDR_Class_Set"].apply(
        lambda s: "; ".join(sorted(s))
    )
)

df["Genomic_MDR_Category"] = (
    df["MDR_Class_Count"].apply(classify)
)

df["Genomic_MDR_3plus"] = (
    df["MDR_Class_Count"] >= 3
)

# Audit possible formatting problems in AMR_Classes.
# Classes should be standardized before final inference.

class_tokens = Counter()

for class_set in df["MDR_Class_Set"]:
    class_tokens.update(class_set)

audit_rows = []

for token, count in sorted(
    class_tokens.items(),
    key=lambda x: (-x[1], x[0])
):
    audit_rows.append({
        "Class_Label": token,
        "Genome_Count": count,
        "Review": (
            "CHECK LABEL"
            if (
                "," in token
                or "|" in token
                or len(token) > 80
            )
            else ""
        )
    })

class_audit = pd.DataFrame(
    audit_rows,
    columns=["Class_Label", "Genome_Count", "Review"]
)

# ============================================================
# SUMMARIES
# ============================================================

category_order = [
    "0 classes",
    "1 class",
    "2 classes",
    "Genomic MDR >=3 classes"
]

overall = (
    df["Genomic_MDR_Category"]
    .value_counts()
    .reindex(category_order, fill_value=0)
    .rename_axis("Classification")
    .reset_index(name="Genome_Count")
)

overall["Percent"] = (
    100 * overall["Genome_Count"] / len(df)
)

species_df = df[
    df["Normalized_Species"].isin(SPECIES)
].copy()

def make_summary(frame, keys):
    result = (
        frame.groupby(keys, dropna=False)
        .agg(
            Total_Genomes=(
                "Assembly_Accession", "size"
            ),
            Genomic_MDR_Genomes=(
                "Genomic_MDR_3plus", "sum"
            ),
            Mean_AMR_Class_Count=(
                "MDR_Class_Count", "mean"
            ),
            Median_AMR_Class_Count=(
                "MDR_Class_Count", "median"
            )
        )
        .reset_index()
    )

    result["Genomic_MDR_Prevalence_Percent"] = (
        100
        * result["Genomic_MDR_Genomes"]
        / result["Total_Genomes"]
    )

    return result.sort_values(
        "Genomic_MDR_Prevalence_Percent",
        ascending=False
    )

species_summary = make_summary(
    species_df,
    ["Normalized_Species"]
)

income_summary = make_summary(
    species_df,
    ["MDR_Income_Group"]
)

species_income_summary = make_summary(
    species_df,
    [
        "Normalized_Species",
        "MDR_Income_Group"
    ]
)

country_summary = make_summary(
    species_df,
    ["Country_Harmonized"]
)

class_count_summary = (
    df["MDR_Class_Count"]
    .value_counts()
    .sort_index()
    .rename_axis("Number_Of_AMR_Classes")
    .reset_index(name="Genomes")
)

# ============================================================
# GENOME-LEVEL FINAL TABLE
# ============================================================

remove_internal = [
    "MDR_Class_Set",
    "Country_Key"
]

genome_table = df.drop(
    columns=remove_internal,
    errors="ignore"
)

# ============================================================
# DATA QUALITY AUDIT
# ============================================================

audit = pd.DataFrame([
    {
        "Check": "Total master genomes",
        "Value": len(df)
    },
    {
        "Check": "Genomes in seven resolved species",
        "Value": len(species_df)
    },
    {
        "Check": "Other/unresolved species",
        "Value": len(df) - len(species_df)
    },
    {
        "Check": "Genome-level MDR >=3 classes",
        "Value": int(df["Genomic_MDR_3plus"].sum())
    },
    {
        "Check": "Genomes with no annotated AMR class",
        "Value": int(
            (df["MDR_Class_Count"] == 0).sum()
        )
    },
    {
        "Check": "Missing or unrecognized income group",
        "Value": int(
            (~df["MDR_Income_Group"].isin(
                ["H", "M", "L+LM"]
            )).sum()
        )
    },
    {
        "Check": "Distinct AMR class labels",
        "Value": len(class_tokens)
    }
])

methods = pd.DataFrame([
    {
        "Parameter": "Classification",
        "Definition": (
            "Genomic MDR proxy: >=3 distinct AMR "
            "classes annotated in one genome"
        )
    },
    {
        "Parameter": "Unit",
        "Definition": "Individual genome / assembly"
    },
    {
        "Parameter": "Primary source",
        "Definition": (
            "AMR_Classes column of final ECC master"
        )
    },
    {
        "Parameter": "Class deduplication",
        "Definition": (
            "Each identical class label counted once "
            "per genome"
        )
    },
    {
        "Parameter": "Interpretation",
        "Definition": (
            "Genomic prediction only; not phenotypically "
            "confirmed MDR"
        )
    },
    {
        "Parameter": "Important caveat",
        "Definition": (
            "Verify whether annotations include intrinsic "
            "resistance, uncertain genes, or synonymous "
            "class labels before manuscript reporting"
        )
    },
    {
        "Parameter": "MDR standard",
        "Definition": (
            "Phenotypic MDR requires acquired "
            "nonsusceptibility in >=3 antimicrobial "
            "categories; gene counts alone are insufficient"
        )
    }
])

# ============================================================
# WRITE EXCEL
# ============================================================

with pd.ExcelWriter(
    OUTPUT,
    engine="openpyxl"
) as writer:

    genome_table.to_excel(
        writer,
        sheet_name="ALL_9157_GENOMES",
        index=False
    )

    genome_table[
        genome_table["Genomic_MDR_3plus"]
    ].to_excel(
        writer,
        sheet_name="GENOMIC_MDR_3PLUS",
        index=False
    )

    overall.to_excel(
        writer,
        sheet_name="OVERALL_CLASSIFICATION",
        index=False
    )

    class_count_summary.to_excel(
        writer,
        sheet_name="NUMBER_OF_CLASSES",
        index=False
    )

    species_summary.to_excel(
        writer,
        sheet_name="SPECIES_SUMMARY",
        index=False
    )

    income_summary.to_excel(
        writer,
        sheet_name="INCOME_SUMMARY",
        index=False
    )

    species_income_summary.to_excel(
        writer,
        sheet_name="SPECIES_INCOME_SUMMARY",
        index=False
    )

    country_summary.to_excel(
        writer,
        sheet_name="COUNTRY_SUMMARY",
        index=False
    )

    class_audit.to_excel(
        writer,
        sheet_name="CLASS_LABEL_AUDIT",
        index=False
    )

    audit.to_excel(
        writer,
        sheet_name="QC_SUMMARY",
        index=False
    )

    methods.to_excel(
        writer,
        sheet_name="METHODS",
        index=False
    )

# ============================================================
# FINAL REPORT
# ============================================================

n_mdr = int(df["Genomic_MDR_3plus"].sum())

print()
print("=" * 100)
print("GENOMIC MDR CLASSIFICATION FINISHED")
print("=" * 100)

print("Total genomes:", len(df))
print("Genomic MDR >=3 classes:", n_mdr)
print(
    "Genomic MDR percentage:",
    f"{100*n_mdr/len(df):.2f}%"
)

print("\nOverall classification:")
print(overall.to_string(index=False))

print("\nSpecies summary:")
print(species_summary.to_string(index=False))

print("\nIncome summary:")
print(income_summary.to_string(index=False))

print("\nClass-label audit:")
print(class_audit.to_string(index=False))

print("\nOutput:", OUTPUT)
print("=" * 100)
