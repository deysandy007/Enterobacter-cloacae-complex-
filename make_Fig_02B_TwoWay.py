
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import font_manager
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MaxNLocator
from scipy.stats import chi2_contingency

# ============================================================
# FIGURE 2B
# ENTEROBACTER CLOACAE COMPLEX:
# NORMALIZED SPECIES COMPOSITION ACROSS INCOME GROUPS
#
# TWO PANELS:
# 1. Absolute genome counts
# 2. Relative species composition (%)
#
# FINAL MASTER: 9157 GENOMES, 2000-2026
#
# INCOME GROUPS: READ FROM FINAL FIGURE 1A EXCEL
# (NO RECLASSIFICATION OR LIVE WORLD BANK DOWNLOAD)
#
# OUTPUTS:
# - ONE 600-DPI JPEG
# - ONE EXCEL WORKBOOK
#
# NO MODIFICATION OF MASTER OR FIGURE 1A
# ============================================================

ROOT = Path(
    r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES"
)

MASTER = (
    ROOT
    / "00_FINAL_COUNTRY_REPORTED_MASTER"
    / "ECC_FINAL_2000_2026_MASTER.xlsx"
)

INCOME_FILE = (
    ROOT
    / "Fig_01A_Global_Geographical_Distribution"
    / "Fig_01A_Country_Income_Database.xlsx"
)

OUT = (
    ROOT / "Fig_02B_Species_Income_Distribution"
)
OUT.mkdir(parents=True, exist_ok=True)

JPG = OUT / "Fig_02B_Species_Income_TwoWay.jpg"
XLSX = OUT / "Fig_02B_Species_Income_Data.xlsx"

EXPECTED_N = 9157

INCOME_ORDER = [
    "H",
    "M",
    "L+LM",
    "Unclassified"
]

INCOME_LABELS = {
    "H": "High income",
    "M": "Upper-middle income",
    "L+LM": "Low + lower-middle income",
    "Unclassified": "Unclassified income"
}

INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82",
    "Unclassified": "#999999"
}

SPECIES_ORDER = [
    "Enterobacter hormaechei",
    "Enterobacter cloacae",
    "Enterobacter roggenkampii",
    "Enterobacter asburiae",
    "Enterobacter kobei",
    "Enterobacter ludwigii",
    "Enterobacter bugandensis",
    "Enterobacter cloacae complex - unresolved"
]

SPECIES_COLORS = {
    "Enterobacter hormaechei": "#0072B2",
    "Enterobacter cloacae": "#D55E00",
    "Enterobacter roggenkampii": "#009E73",
    "Enterobacter asburiae": "#CC79A7",
    "Enterobacter kobei": "#7851A9",
    "Enterobacter ludwigii": "#A6761D",
    "Enterobacter bugandensis": "#56B4E9",
    "Enterobacter cloacae complex - unresolved":
        "#808080"
}

SPECIES_SHORT = {
    "Enterobacter hormaechei": "E. hormaechei",
    "Enterobacter cloacae": "E. cloacae",
    "Enterobacter roggenkampii": "E. roggenkampii",
    "Enterobacter asburiae": "E. asburiae",
    "Enterobacter kobei": "E. kobei",
    "Enterobacter ludwigii": "E. ludwigii",
    "Enterobacter bugandensis": "E. bugandensis",
    "Enterobacter cloacae complex - unresolved":
        "Unresolved ECC"
}

# ============================================================
# FONT
# ============================================================

available = {
    f.name for f in font_manager.fontManager.ttflist
}

FONT = next(
    (
        f for f in (
            "Arial",
            "Helvetica",
            "Liberation Sans",
            "DejaVu Sans"
        )
        if f in available
    ),
    "DejaVu Sans"
)

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 11,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "axes.linewidth": 0.75,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})

# ============================================================
# 1. LOAD FINAL MASTER
# ============================================================

print("=" * 100)
print("FIGURE 2B - SPECIES DISTRIBUTION BY INCOME GROUP")
print("=" * 100)

if not MASTER.exists():
    raise FileNotFoundError(MASTER)

if not INCOME_FILE.exists():
    raise FileNotFoundError(
        "Final Figure 1A income workbook not found:\n"
        + str(INCOME_FILE)
    )

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

for column in required:
    if column not in df.columns:
        raise ValueError("Missing master column: " + column)

N = len(df)

if N != EXPECTED_N:
    raise ValueError(
        f"Expected {EXPECTED_N} genomes; found {N}"
    )

if df["Assembly_Accession"].isna().any():
    raise ValueError("Missing accession")

if df["Assembly_Accession"].duplicated().any():
    raise ValueError("Duplicate accession")

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

df["Country_Harmonized"] = (
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

unexpected_species = (
    set(df["Species"]) - set(SPECIES_ORDER)
)

if unexpected_species:
    raise ValueError(
        "Unexpected normalized species: "
        + str(sorted(unexpected_species))
    )

print("Final genomes:", N)
print("Collection period:", df.Year.min(), "-", df.Year.max())
print("Normalized species:", df.Species.nunique())

# ============================================================
# 2. LOAD FROZEN FIGURE 1A INCOME DATABASE
# ============================================================

print("\nLoading saved Figure 1A income classifications...")

income = pd.read_excel(
    INCOME_FILE,
    sheet_name="COUNTRY_INCOME_COUNTS",
    engine="openpyxl"
)

income_required = [
    "Country_Harmonized",
    "Income_Group"
]

for column in income_required:
    if column not in income.columns:
        raise ValueError(
            "Missing Figure 1A income column: " + column
        )

income["Country_Harmonized"] = (
    income["Country_Harmonized"]
    .fillna("")
    .astype(str)
    .str.strip()
)

income["Income_Group"] = (
    income["Income_Group"]
    .fillna("Unclassified")
    .astype(str)
    .str.strip()
)

if income["Country_Harmonized"].duplicated().any():
    raise ValueError(
        "Duplicate country in Figure 1A income table"
    )

unknown_income = (
    set(income["Income_Group"]) - set(INCOME_ORDER)
)

if unknown_income:
    raise ValueError(
        "Unexpected income group: "
        + str(sorted(unknown_income))
    )

master_countries = set(df["Country_Harmonized"])
income_countries = set(income["Country_Harmonized"])

missing_countries = sorted(
    master_countries - income_countries
)

if missing_countries:
    raise ValueError(
        "Countries missing from saved Figure 1A:\n"
        + "\n".join(missing_countries)
    )

# Check that Figure 1A contains the same genome totals.
fig1a_counts = (
    income.set_index("Country_Harmonized")[
        "Genome_Count"
    ]
    if "Genome_Count" in income.columns
    else None
)

if fig1a_counts is not None:
    current_counts = (
        df["Country_Harmonized"].value_counts()
    )

    comparison = (
        pd.DataFrame({
            "Current_Master": current_counts,
            "Figure_1A": fig1a_counts
        })
        .fillna(0)
    )

    if not np.array_equal(
        comparison["Current_Master"].to_numpy(),
        comparison["Figure_1A"].to_numpy()
    ):
        raise ValueError(
            "Figure 1A country counts do not match "
            "the final master. Recheck Figure 1A."
        )

use_columns = [
    "Country_Harmonized",
    "Income_Group"
]

for optional in [
    "Income_Source",
    "Income_Status",
    "WB_Code",
    "WB_Income_Code",
    "ISO3"
]:
    if optional in income.columns:
        use_columns.append(optional)

income_mapping = income[use_columns].copy()

df = df.merge(
    income_mapping,
    on="Country_Harmonized",
    how="left",
    validate="many_to_one"
)

if df["Income_Group"].isna().any():
    raise ValueError("Income mapping failed")

assert len(df) == N

print("Country-income mapping validated")

# ============================================================
# 3. BUILD CONTINGENCY TABLE
# ============================================================

groups = [
    g for g in INCOME_ORDER
    if df["Income_Group"].eq(g).any()
]

counts = pd.crosstab(
    df["Income_Group"],
    df["Species"]
).reindex(
    index=groups,
    columns=SPECIES_ORDER,
    fill_value=0
)

assert int(counts.to_numpy().sum()) == N

income_totals = counts.sum(axis=1)
species_totals = counts.sum(axis=0)

percent = (
    counts.div(income_totals, axis=0) * 100
)

assert np.allclose(
    percent.sum(axis=1).to_numpy(),
    100.0
)

print("\nINCOME GROUP GENOME COUNTS")
print(income_totals.to_string())

print("\nSPECIES BY INCOME COUNTS")
print(counts.to_string())

print("\nSPECIES BY INCOME PERCENTAGES")
print(percent.round(2).to_string())

# ============================================================
# 4. STATISTICAL ANALYSIS
#
# A) Four-group descriptive table retains unclassified.
# B) Inferential chi-square uses classified groups only,
#    because unclassified is not a socioeconomic level.
# ============================================================

classified = [
    g for g in ["H", "M", "L+LM"]
    if g in counts.index
]

stats_rows = []
expected_excel = pd.DataFrame()
residual_excel = pd.DataFrame()

if len(classified) >= 2:

    test_counts = counts.loc[classified].copy()

    # Exclude columns absent entirely from classified data.
    test_counts = test_counts.loc[
        :, test_counts.sum(axis=0) > 0
    ]

    n_test = int(test_counts.to_numpy().sum())

    chi2, p_value, dof, expected = chi2_contingency(
        test_counts.to_numpy(),
        correction=False
    )

    effect = float(np.sqrt(
        chi2 / (
            n_test * min(
                test_counts.shape[0] - 1,
                test_counts.shape[1] - 1
            )
        )
    ))

    below5 = int((expected < 5).sum())
    below1 = int((expected < 1).sum())
    total_cells = expected.size

    expected_excel = pd.DataFrame(
        expected,
        index=test_counts.index,
        columns=test_counts.columns
    )

    residual_excel = pd.DataFrame(
        (
            test_counts.to_numpy() - expected
        ) / np.sqrt(expected),
        index=test_counts.index,
        columns=test_counts.columns
    )

    stats_rows = [
        ("Test", "Pearson chi-square"),
        ("Analyzed income groups", ", ".join(classified)),
        ("Genomes in statistical test", n_test),
        ("Excluded unclassified genomes",
         int(N - n_test)),
        ("Chi-square statistic", float(chi2)),
        ("Degrees of freedom", int(dof)),
        ("P-value", float(p_value)),
        ("Cramers V", effect),
        ("Expected cells <5", below5),
        ("Expected cells <1", below1),
        ("Expected cells <5 (%)",
         100 * below5 / total_cells),
        ("Test limitation",
         "Descriptive association in a nonrandom, "
         "country-clustered public genomic dataset")
    ]

    print("\nSTATISTICAL ANALYSIS")
    print("Chi-square:", round(chi2, 4))
    print("Degrees of freedom:", dof)
    print("P-value:", p_value)
    print("Cramers V:", round(effect, 4))
    print("Expected cells <5:", below5)

else:
    stats_rows = [
        ("Test", "Not performed"),
        ("Reason",
         "Fewer than two classified income groups")
    ]

stats_df = pd.DataFrame(
    stats_rows,
    columns=["Metric", "Value"]
)

# ============================================================
# 5. COUNTRY-LEVEL SAMPLING AUDIT
# ============================================================

country_counts = (
    df.groupby([
        "Income_Group",
        "Country_Harmonized"
    ])
    .size()
    .reset_index(name="Genome_Count")
)

country_counts["Income_Total"] = (
    country_counts["Income_Group"]
    .map(income_totals)
)

country_counts["Within_Income_Percent"] = (
    100
    * country_counts["Genome_Count"]
    / country_counts["Income_Total"]
).round(3)

country_counts = country_counts.sort_values(
    ["Income_Group", "Genome_Count"],
    ascending=[True, False]
)

country_species = (
    df.groupby([
        "Income_Group",
        "Country_Harmonized",
        "Species"
    ])
    .size()
    .reset_index(name="Genome_Count")
)

country_species = country_species.merge(
    country_counts[[
        "Income_Group",
        "Country_Harmonized",
        "Genome_Count"
    ]].rename(
        columns={"Genome_Count": "Country_Total"}
    ),
    on=["Income_Group", "Country_Harmonized"],
    validate="many_to_one"
)

country_species["Within_Country_Percent"] = (
    100
    * country_species["Genome_Count"]
    / country_species["Country_Total"]
).round(3)

# Sensitivity display: give each represented country
# equal weight, instead of each genome equal weight.
#
# This is descriptive, not a replacement for formal
# country-clustered statistical modeling.

country_matrix = pd.crosstab(
    [
        df["Income_Group"],
        df["Country_Harmonized"]
    ],
    df["Species"]
).reindex(
    columns=SPECIES_ORDER,
    fill_value=0
)

country_row_percent = country_matrix.div(
    country_matrix.sum(axis=1),
    axis=0
) * 100

equal_country = (
    country_row_percent
    .groupby(level=0)
    .mean()
    .reindex(groups)
)

equal_country = equal_country.fillna(0)

income_country_summary = []

for group in groups:

    country_subset = country_counts[
        country_counts["Income_Group"] == group
    ].sort_values(
        "Genome_Count",
        ascending=False
    )

    income_country_summary.append({
        "Income_Group": group,
        "Total_Genomes": int(income_totals[group]),
        "Number_of_Countries": len(country_subset),
        "Largest_Country": (
            country_subset.iloc[0][
                "Country_Harmonized"
            ]
        ),
        "Largest_Country_Genomes": int(
            country_subset.iloc[0]["Genome_Count"]
        ),
        "Largest_Country_Percent": round(
            100
            * country_subset.iloc[0]["Genome_Count"]
            / income_totals[group],
            3
        ),
        "Top_3_Country_Percent": round(
            100
            * country_subset.head(3)[
                "Genome_Count"
            ].sum()
            / income_totals[group],
            3
        )
    })

income_country_summary = pd.DataFrame(
    income_country_summary
)

# ============================================================
# 6. FIGURE 2B - TWO PANEL DESIGN
# ============================================================

print("\nGenerating Nature-style two-way figure...")

positions = np.arange(len(groups))
labels = [
    f"{INCOME_LABELS[g]}\n(n={income_totals[g]:,})"
    for g in groups
]

fig, (ax1, ax2) = plt.subplots(
    1, 2,
    figsize=(17.2, 9.5),
    gridspec_kw={
        "width_ratios": [1.09, 1.00],
        "wspace": 0.29
    }
)

fig.suptitle(
    "Enterobacter cloacae complex species "
    "composition across country-income groups",
    fontsize=18,
    fontweight="bold",
    y=0.975
)

fig.text(
    0.5, 0.932,
    f"{N:,} genomes | Collection period 2000–2026",
    ha="center",
    fontsize=12,
    color="#555555"
)

BAR_HEIGHT = 0.53

# ============================================================
# PANEL 1: ABSOLUTE COUNTS
# ============================================================

left = np.zeros(len(groups))

for species in SPECIES_ORDER:

    values = counts[species].to_numpy()

    ax1.barh(
        positions,
        values,
        left=left,
        height=BAR_HEIGHT,
        color=SPECIES_COLORS[species],
        edgecolor="white",
        linewidth=0.35
    )

    left += values

ax1.set_title(
    "Absolute genome counts",
    fontsize=14,
    fontweight="bold",
    loc="left",
    pad=13
)

ax1.set_xlabel("Number of genomes")

ax1.set_yticks(positions)
ax1.set_yticklabels(labels, fontsize=11)

ax1.invert_yaxis()

ax1.set_xlim(
    0,
    max(income_totals) * 1.13
)

ax1.xaxis.set_major_formatter(
    FuncFormatter(
        lambda x, pos: f"{x:,.0f}"
    )
)

for i, group in enumerate(groups):
    n = int(income_totals[group])

    ax1.text(
        n + max(income_totals) * 0.016,
        positions[i],
        f"{n:,}",
        va="center",
        ha="left",
        fontsize=11,
        fontweight="bold",
        color="#333333"
    )

# ============================================================
# PANEL 2: RELATIVE COMPOSITION
# ============================================================

left = np.zeros(len(groups))

for species in SPECIES_ORDER:

    values = percent[species].to_numpy()

    ax2.barh(
        positions,
        values,
        left=left,
        height=BAR_HEIGHT,
        color=SPECIES_COLORS[species],
        edgecolor="white",
        linewidth=0.35
    )

    for i, val in enumerate(values):

        if val >= 9:
            ax2.text(
                left[i] + val / 2,
                positions[i],
                f"{val:.1f}%",
                ha="center",
                va="center",
                fontsize=9.7,
                fontweight="bold",
                color="white",
                clip_on=True
            )

    left += values

ax2.set_title(
    "Relative species composition",
    fontsize=14,
    fontweight="bold",
    loc="left",
    pad=13
)

ax2.set_xlim(0, 100)

ax2.set_xlabel(
    "Species composition (%)"
)

ax2.set_xticks(np.arange(0, 101, 20))

ax2.xaxis.set_major_formatter(
    FuncFormatter(
        lambda x, pos: f"{x:.0f}%"
    )
)

ax2.set_yticks(positions)

ax2.set_yticklabels(
    [INCOME_LABELS[g] for g in groups],
    fontsize=11
)

ax2.invert_yaxis()

# ============================================================
# PANEL STYLING
# ============================================================

for ax in [ax1, ax2]:

    ax.set_axisbelow(True)

    ax.grid(
        axis="x",
        color="#E5E5E5",
        linewidth=0.55,
        alpha=0.75
    )

    for name in [
        "top", "right", "left"
    ]:
        ax.spines[name].set_visible(False)

    ax.spines["bottom"].set_color("#777777")

    ax.tick_params(
        axis="y",
        length=0,
        pad=7
    )

    ax.tick_params(
        axis="x",
        length=3,
        width=0.7,
        colors="#333333"
    )

# ============================================================
# SHARED LEGEND
# ============================================================

handles = [
    Patch(
        facecolor=SPECIES_COLORS[species],
        edgecolor="none",
        label=SPECIES_SHORT[species]
    )
    for species in SPECIES_ORDER
]

fig.legend(
    handles=handles,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.115),
    ncol=4,
    frameon=False,
    fontsize=11,
    columnspacing=1.7,
    handlelength=1.3,
    labelspacing=1.2
)

fig.subplots_adjust(
    left=0.18,
    right=0.965,
    top=0.85,
    bottom=0.245
)

fig.text(
    0.5,
    0.052,
    "Income classifications match Figure 1A. "
    "Percentages are calculated within each income group; "
    "they represent sampled genome composition, not "
    "population prevalence.",
    ha="center",
    fontsize=10,
    color="#555555"
)

# ============================================================
# 7. SAVE FIGURE
# ============================================================

fig.savefig(
    JPG,
    format="jpeg",
    dpi=600,
    facecolor="white",
    pil_kwargs={
        "quality": 95,
        "subsampling": 0,
        "optimize": True
    }
)

plt.close(fig)

# ============================================================
# 8. PREPARE EXCEL SOURCE DATA
# ============================================================

count_wide = (
    counts.rename(
        columns=SPECIES_SHORT
    ).copy()
)

count_wide.insert(
    0,
    "Total_Genomes",
    income_totals
)

count_wide = (
    count_wide
    .rename_axis("Income_Group")
    .reset_index()
)

percent_wide = (
    percent.rename(
        columns=SPECIES_SHORT
    )
    .round(4)
    .rename_axis("Income_Group")
    .reset_index()
)

count_long = (
    counts.rename_axis("Income_Group")
    .reset_index()
    .melt(
        id_vars="Income_Group",
        var_name="Normalized_Species",
        value_name="Genome_Count"
    )
)

percent_long = (
    percent.rename_axis("Income_Group")
    .reset_index()
    .melt(
        id_vars="Income_Group",
        var_name="Normalized_Species",
        value_name="Percentage"
    )
)

combined = count_long.merge(
    percent_long,
    on=["Income_Group", "Normalized_Species"],
    validate="one_to_one"
)

combined["Percentage"] = (
    combined["Percentage"].round(4)
)

species_summary = (
    species_totals
    .rename_axis("Normalized_Species")
    .reset_index(name="Genome_Count")
)

species_summary["Percentage"] = (
    100
    * species_summary["Genome_Count"]
    / N
).round(4)

income_summary = (
    income_totals
    .rename_axis("Income_Group")
    .reset_index(name="Genome_Count")
)

income_summary["Percentage"] = (
    100
    * income_summary["Genome_Count"]
    / N
).round(4)

equal_country_excel = (
    equal_country.rename_axis("Income_Group")
    .reset_index()
)

expected_excel = (
    expected_excel.rename_axis("Income_Group")
    .reset_index()
)

residual_excel = (
    residual_excel.rename_axis("Income_Group")
    .reset_index()
)

metadata = pd.DataFrame([
    ("Input master", str(MASTER)),
    ("Master genomes", N),
    ("Study years", "2000-2026"),
    ("Species variable", "Normalized_Species"),
    ("Income reference", str(INCOME_FILE)),
    ("Income assignment",
     "Exact classifications from final Figure 1A"),
    ("Percentage denominator",
     "Total genomes within income group"),
    ("Unclassified income",
     "Shown in descriptive plots, excluded from chi-square"),
    ("Statistical test", "Pearson chi-square"),
    ("Effect size", "Cramers V"),
    ("Sampling sensitivity",
     "Equal-country-weighted species composition"),
    ("Sampling interpretation",
     "Nonrandom NCBI genomic representation"),
    ("Figure format", "JPEG"),
    ("Figure resolution", "600 DPI"),
    ("Figure font", FONT),
    ("Created", datetime.now().isoformat())
], columns=["Parameter", "Value"])

# ============================================================
# 9. WRITE EXCEL WORKBOOK
# ============================================================

print("\nSaving Excel source database...")

with pd.ExcelWriter(
    XLSX,
    engine="openpyxl"
) as writer:

    count_wide.to_excel(
        writer,
        sheet_name="ABSOLUTE_COUNTS",
        index=False
    )

    percent_wide.to_excel(
        writer,
        sheet_name="RELATIVE_PERCENTAGES",
        index=False
    )

    combined.to_excel(
        writer,
        sheet_name="COUNTS_AND_PERCENTAGES",
        index=False
    )

    species_summary.to_excel(
        writer,
        sheet_name="SPECIES_TOTALS",
        index=False
    )

    income_summary.to_excel(
        writer,
        sheet_name="INCOME_TOTALS",
        index=False
    )

    stats_df.to_excel(
        writer,
        sheet_name="STATISTICAL_SUMMARY",
        index=False
    )

    expected_excel.to_excel(
        writer,
        sheet_name="EXPECTED_CELL_COUNTS",
        index=False
    )

    residual_excel.to_excel(
        writer,
        sheet_name="PEARSON_RESIDUALS",
        index=False
    )

    country_counts.to_excel(
        writer,
        sheet_name="COUNTRY_SAMPLE_COUNTS",
        index=False
    )

    income_country_summary.to_excel(
        writer,
        sheet_name="COUNTRY_IMBALANCE",
        index=False
    )

    country_species.to_excel(
        writer,
        sheet_name="COUNTRY_SPECIES",
        index=False
    )

    equal_country_excel.to_excel(
        writer,
        sheet_name="EQUAL_COUNTRY_SENSITIVITY",
        index=False
    )

    df[[
        "Assembly_Accession",
        "Country_Harmonized",
        "Year",
        "Normalized_Species",
        "Income_Group"
    ]].to_excel(
        writer,
        sheet_name="GENOME_SOURCE_DATA",
        index=False
    )

    metadata.to_excel(
        writer,
        sheet_name="FIGURE_METADATA",
        index=False
    )

# ============================================================
# 10. VERIFY OUTPUT
# ============================================================

check = pd.read_excel(
    XLSX,
    sheet_name="ABSOLUTE_COUNTS",
    engine="openpyxl"
)

if int(check["Total_Genomes"].sum()) != N:
    raise RuntimeError(
        "Excel genome count mismatch"
    )

if not JPG.exists() or JPG.stat().st_size < 100000:
    raise RuntimeError(
        "JPEG missing or unexpectedly small"
    )

print("\n" + "=" * 100)
print("FIGURE 2B COMPLETED SUCCESSFULLY")
print("=" * 100)

print("\nINCOME COUNTS")
print(income_summary.to_string(index=False))

print("\nSPECIES COUNTS")
print(species_summary.to_string(index=False))

print("\nSTATISTICAL SUMMARY")
print(stats_df.to_string(index=False))

print("\nCOUNTRY IMBALANCE")
print(income_country_summary.to_string(index=False))

print("\nFINAL OUTPUTS")
print(JPG)
print(XLSX)

print("\nFigure: two horizontal stacked-bar panels")
print("Resolution: 600 DPI")
print("Excel workbook: 14 worksheets")
print("Final master: unchanged")
print("Figure 1A income assignments: unchanged")
print("=" * 100)
