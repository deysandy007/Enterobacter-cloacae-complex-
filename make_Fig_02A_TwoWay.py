
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import font_manager
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
from scipy.stats import chi2_contingency

# ============================================================
# FIGURE 2A
# ECC SPECIES DISTRIBUTION ACROSS ONE HEALTH SOURCES
#
# PANEL 1: ABSOLUTE GENOME COUNTS
# PANEL 2: RELATIVE SPECIES COMPOSITION (%)
#
# FINAL DATASET: 9157 GENOMES, 2000-2026
# NORMALIZED SPECIES: 8 CATEGORIES
#
# OUTPUT:
# 1 JPEG (600 DPI)
# 1 EXCEL SOURCE DATA
#
# THE FINAL MASTER IS READ ONLY.
# ============================================================

ROOT = Path(
    r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES"
)

MASTER = (
    ROOT / "00_FINAL_COUNTRY_REPORTED_MASTER"
    / "ECC_FINAL_2000_2026_MASTER.xlsx"
)

OUT = (
    ROOT / "Fig_02A_Species_OneHealth_Distribution"
)

OUT.mkdir(parents=True, exist_ok=True)

JPG = OUT / "Fig_02A_Species_OneHealth_TwoWay.jpg"
XLSX = OUT / "Fig_02A_Species_OneHealth_Data.xlsx"

EXPECTED_N = 9157

# ============================================================
# GROUP DEFINITIONS
# ============================================================

SOURCES = [
    "Human",
    "Animal",
    "Environment",
    "Food",
    "Plant-associated"
]

SPECIES = [
    "Enterobacter hormaechei",
    "Enterobacter cloacae",
    "Enterobacter roggenkampii",
    "Enterobacter asburiae",
    "Enterobacter kobei",
    "Enterobacter ludwigii",
    "Enterobacter bugandensis",
    "Enterobacter cloacae complex - unresolved"
]

SHORT_NAMES = {
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

# Distinct colors consistent across both panels.

COLORS = {
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

# ============================================================
# PUBLICATION FONT
# ============================================================

installed = {
    f.name for f in font_manager.fontManager.ttflist
}

FONT = next(
    (
        f for f in [
            "Arial",
            "Helvetica",
            "Liberation Sans",
            "DejaVu Sans"
        ]
        if f in installed
    ),
    "DejaVu Sans"
)

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 11,
    "axes.titlesize": 15,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 12,
    "axes.linewidth": 0.75,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})

# ============================================================
# 1. READ FINAL MASTER
# ============================================================

print("=" * 95)
print("FIGURE 2A - ONE HEALTH SPECIES DISTRIBUTION")
print("=" * 95)

if not MASTER.exists():
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
    "Final_Source_Category",
    "Verified_Collection_Year"
]

for column in required:
    if column not in df.columns:
        raise ValueError(
            f"Missing master column: {column}"
        )

N = len(df)

if N != EXPECTED_N:
    raise ValueError(
        f"Expected {EXPECTED_N} genomes, found {N}"
    )

if df["Assembly_Accession"].duplicated().any():
    raise ValueError(
        "Duplicate assembly accession detected"
    )

years = pd.to_numeric(
    df["Verified_Collection_Year"],
    errors="coerce"
)

if (
    years.isna().any()
    or not years.between(2000, 2026).all()
    or not (years % 1 == 0).all()
):
    raise ValueError("Invalid collection years")

df["Species"] = (
    df["Normalized_Species"]
    .fillna("")
    .astype(str)
    .str.strip()
)

df["Source"] = (
    df["Final_Source_Category"]
    .fillna("")
    .astype(str)
    .str.strip()
)

unexpected_species = sorted(
    set(df["Species"]) - set(SPECIES)
)

unexpected_sources = sorted(
    set(df["Source"]) - set(SOURCES)
)

if unexpected_species:
    raise ValueError(
        "Unrecognized normalized species:\n"
        + str(unexpected_species)
    )

if unexpected_sources:
    raise ValueError(
        "Unrecognized One Health sources:\n"
        + str(unexpected_sources)
    )

print("Final genomes:", N)
print("Valid study period: 2000-2026")
print("Normalized species:", df["Species"].nunique())
print("One Health compartments:", df["Source"].nunique())

# ============================================================
# 2. CONTINGENCY TABLE
# ============================================================

counts = pd.crosstab(
    df["Source"],
    df["Species"]
).reindex(
    index=SOURCES,
    columns=SPECIES,
    fill_value=0
)

if (counts.sum(axis=0) == 0).any():
    raise ValueError(
        "An expected species category has zero genomes"
    )

if (counts.sum(axis=1) == 0).any():
    raise ValueError(
        "An expected source category has zero genomes"
    )

assert int(counts.to_numpy().sum()) == N

source_totals = counts.sum(axis=1)
species_totals = counts.sum(axis=0)

percent = counts.div(
    source_totals,
    axis=0
) * 100.0

assert np.allclose(
    percent.sum(axis=1).to_numpy(),
    100.0
)

print("\nSOURCE TOTALS")
print(source_totals.to_string())

print("\nSPECIES TOTALS")
print(species_totals.to_string())

print("\nABSOLUTE CONTINGENCY TABLE")
print(counts.to_string())

print("\nPERCENTAGE COMPOSITION")
print(percent.round(2).to_string())

# ============================================================
# 3. STATISTICS
# ============================================================

chi2, p, dof, expected = chi2_contingency(
    counts.values,
    correction=False
)

rows, cols = counts.shape

cramers_v = np.sqrt(
    chi2 / (
        N * min(rows - 1, cols - 1)
    )
)

expected_below_5 = int(
    (expected < 5).sum()
)

expected_below_1 = int(
    (expected < 1).sum()
)

fraction_below_5 = (
    100 * expected_below_5 / expected.size
)

stats = pd.DataFrame([
    ("Total genomes", N),
    ("Number of source categories", rows),
    ("Number of normalized species", cols),
    ("Chi-square statistic", chi2),
    ("Degrees of freedom", dof),
    ("P-value", p),
    ("Cramers V", cramers_v),
    ("Expected cells below 5", expected_below_5),
    ("Expected cells below 1", expected_below_1),
    ("Expected cells below 5 (%)", fraction_below_5),
    ("Statistical interpretation",
     "Association within sampled genomes, "
     "not a population prevalence comparison"),
    ("Assumption note",
     "Review expected counts and sample independence")
], columns=[
    "Metric", "Value"
])

# ============================================================
# 4. PREPARE PLOTTING ARRAYS
# ============================================================

species_colors = [
    COLORS[s] for s in SPECIES
]

source_y = np.arange(len(SOURCES))

source_labels = [
    f"{source} (n={source_totals[source]:,})"
    for source in SOURCES
]

# ============================================================
# 5. PUBLICATION FIGURE
# ============================================================

print("\nRendering Figure 2A...")

fig, (ax1, ax2) = plt.subplots(
    1, 2,
    figsize=(16.5, 9.2),
    gridspec_kw={
        "width_ratios": [1.10, 1.0],
        "wspace": 0.31
    }
)

fig.suptitle(
    "Species composition of the Enterobacter "
    "cloacae complex across One Health sources",
    fontsize=18,
    fontweight="bold",
    y=0.978
)

fig.text(
    0.5, 0.932,
    f"{N:,} species-normalized genomes "
    f"| Collection period 2000–2026",
    ha="center",
    fontsize=12,
    color="#555555"
)

BAR_HEIGHT = 0.60

# ============================================================
# 6. ABSOLUTE COUNT STACKED BARS
# ============================================================

left = np.zeros(len(SOURCES), dtype=float)

for species in SPECIES:

    values = counts[species].to_numpy(dtype=float)

    ax1.barh(
        source_y,
        values,
        left=left,
        height=BAR_HEIGHT,
        color=COLORS[species],
        edgecolor="white",
        linewidth=0.35
    )

    left += values

ax1.set_title(
    "Absolute genome counts",
    loc="left",
    fontsize=15,
    fontweight="bold",
    pad=13
)

ax1.set_xlabel(
    "Number of genomes",
    fontsize=12
)

ax1.set_yticks(source_y)
ax1.set_yticklabels(
    source_labels,
    fontsize=12
)

ax1.invert_yaxis()

ax1.set_xlim(
    0,
    max(source_totals) * 1.12
)

ax1.xaxis.set_major_formatter(
    FuncFormatter(
        lambda x, pos: f"{x:,.0f}"
    )
)

# Number of genomes at the end of each bar
for i, source in enumerate(SOURCES):
    ax1.text(
        source_totals[source] + max(source_totals) * 0.012,
        i,
        f"{source_totals[source]:,}",
        fontsize=11,
        va="center",
        ha="left",
        fontweight="bold",
        color="#303030"
    )

# ============================================================
# 7. RELATIVE COMPOSITION
# ============================================================

left = np.zeros(len(SOURCES), dtype=float)

for species in SPECIES:

    values = percent[species].to_numpy(dtype=float)

    ax2.barh(
        source_y,
        values,
        left=left,
        height=BAR_HEIGHT,
        color=COLORS[species],
        edgecolor="white",
        linewidth=0.35
    )

    # Display percentages for sufficiently large segments.
    for i, val in enumerate(values):

        if val >= 8.0:

            midpoint = left[i] + val / 2

            ax2.text(
                midpoint,
                source_y[i],
                f"{val:.1f}%",
                ha="center",
                va="center",
                fontsize=10,
                fontweight="bold",
                color="white",
                clip_on=True
            )

    left += values

ax2.set_title(
    "Relative species composition",
    loc="left",
    fontsize=15,
    fontweight="bold",
    pad=13
)

ax2.set_xlabel(
    "Species composition (%)",
    fontsize=12
)

ax2.set_xlim(0, 100)

ax2.set_xticks(
    np.arange(0, 101, 20)
)

ax2.xaxis.set_major_formatter(
    FuncFormatter(
        lambda x, pos: f"{x:.0f}%"
    )
)

ax2.set_yticks(source_y)
ax2.set_yticklabels(SOURCES, fontsize=12)
ax2.invert_yaxis()

# ============================================================
# 8. COMMON AXIS DESIGN
# ============================================================

for ax in [ax1, ax2]:

    ax.set_axisbelow(True)

    ax.grid(
        axis="x",
        color="#E6E6E6",
        linewidth=0.6,
        alpha=0.8
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    ax.spines["bottom"].set_color("#666666")

    ax.tick_params(
        axis="both",
        which="major",
        length=3,
        width=0.7,
        colors="#333333"
    )

    ax.tick_params(
        axis="y",
        length=0,
        pad=8
    )

# ============================================================
# 9. SHARED SPECIES LEGEND
# ============================================================

legend_handles = [
    Patch(
        facecolor=COLORS[s],
        edgecolor="none",
        label=SHORT_NAMES[s]
    )
    for s in SPECIES
]

fig.legend(
    handles=legend_handles,
    loc="lower center",
    bbox_to_anchor=(0.50, 0.085),
    ncol=4,
    frameon=False,
    fontsize=11,
    columnspacing=1.8,
    handlelength=1.3,
    labelspacing=1.1
)

# Allow space for a shared legend.
fig.subplots_adjust(
    left=0.17,
    right=0.965,
    top=0.85,
    bottom=0.23
)

fig.text(
    0.5, 0.031,
    "Species were standardized from NCBI metadata. "
    "Relative abundances describe the composition "
    "of sampled genomes, not population prevalence.",
    ha="center",
    fontsize=10,
    color="#555555"
)

# ============================================================
# 10. EXPORT 600-DPI JPEG
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
# 11. FORMAT EXCEL TABLES
# ============================================================

count_excel = counts.rename(
    columns=SHORT_NAMES
).copy()

count_excel.insert(
    0,
    "Total_Genomes",
    source_totals
)

count_excel = (
    count_excel
    .rename_axis("Source")
    .reset_index()
)

percentage_excel = (
    percent
    .rename(columns=SHORT_NAMES)
    .round(4)
    .rename_axis("Source")
    .reset_index()
)

count_long = (
    counts.rename_axis("Source")
    .reset_index()
    .melt(
        id_vars="Source",
        var_name="Normalized_Species",
        value_name="Genome_Count"
    )
)

percent_long = (
    percent.rename_axis("Source")
    .reset_index()
    .melt(
        id_vars="Source",
        var_name="Normalized_Species",
        value_name="Percentage"
    )
)

combined_long = count_long.merge(
    percent_long,
    on=["Source", "Normalized_Species"],
    validate="one_to_one"
)

combined_long["Percentage"] = (
    combined_long["Percentage"].round(4)
)

species_summary = (
    species_totals
    .rename_axis("Normalized_Species")
    .reset_index(name="Genome_Count")
)

species_summary["Percentage"] = (
    100 * species_summary["Genome_Count"] / N
).round(4)

source_summary = (
    source_totals
    .rename_axis("Source")
    .reset_index(name="Genome_Count")
)

source_summary["Percentage"] = (
    100 * source_summary["Genome_Count"] / N
).round(4)

expected_df = pd.DataFrame(
    expected,
    index=SOURCES,
    columns=[SHORT_NAMES[s] for s in SPECIES]
).rename_axis("Source").reset_index()

metadata = pd.DataFrame([
    ("Input master", str(MASTER)),
    ("Final genomes", N),
    ("Collection years", "2000-2026"),
    ("Genome uniqueness", "Assembly accession"),
    ("Species classification", "Normalized_Species"),
    ("One Health classification", "Final_Source_Category"),
    ("Absolute abundance", "Genome counts"),
    ("Relative composition", "Percent within each source"),
    ("Normalization denominator",
     "Total genomes in that One Health source"),
    ("Statistic", "Pearson chi-square"),
    ("Effect size", "Cramers V"),
    ("Image format", "JPEG"),
    ("Resolution", "600 DPI"),
    ("Font", FONT),
    ("Created", datetime.now().isoformat())
], columns=["Parameter", "Value"])

# ============================================================
# 12. SAVE EXCEL SOURCE DATA
# ============================================================

with pd.ExcelWriter(
    XLSX,
    engine="openpyxl"
) as writer:

    count_excel.to_excel(
        writer,
        sheet_name="ABSOLUTE_COUNTS",
        index=False
    )

    percentage_excel.to_excel(
        writer,
        sheet_name="RELATIVE_PERCENTAGES",
        index=False
    )

    combined_long.to_excel(
        writer,
        sheet_name="COUNTS_AND_PERCENTAGES",
        index=False
    )

    species_summary.to_excel(
        writer,
        sheet_name="SPECIES_TOTALS",
        index=False
    )

    source_summary.to_excel(
        writer,
        sheet_name="SOURCE_TOTALS",
        index=False
    )

    stats.to_excel(
        writer,
        sheet_name="STATISTICAL_SUMMARY",
        index=False
    )

    expected_df.to_excel(
        writer,
        sheet_name="EXPECTED_CELL_COUNTS",
        index=False
    )

    df[[
        "Assembly_Accession",
        "Normalized_Species",
        "Final_Source_Category",
        "Verified_Collection_Year"
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
# 13. FINAL VALIDATION
# ============================================================

check_counts = pd.read_excel(
    XLSX,
    sheet_name="ABSOLUTE_COUNTS"
)

if int(check_counts["Total_Genomes"].sum()) != N:
    raise RuntimeError("Excel genome count mismatch")

if not JPG.exists() or not XLSX.exists():
    raise RuntimeError("Output file missing")

print("\n" + "=" * 95)
print("FIGURE 2A COMPLETED")
print("=" * 95)

print("\nGENOMES:", N)

print("\nSOURCE DISTRIBUTION")
print(source_summary.to_string(index=False))

print("\nSPECIES DISTRIBUTION")
print(species_summary.to_string(index=False))

print("\nSTATISTICAL RESULTS")
print(stats.to_string(index=False))

print("\nFINAL OUTPUTS")
print(JPG)
print(XLSX)

print("\nFigure layout: 2 horizontal stacked-bar panels")
print("Resolution: 600 DPI JPEG")
print("Excel source data: 9 worksheets")
print("Original master unchanged")
print("=" * 95)
