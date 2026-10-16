
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.colors import (
    LinearSegmentedColormap,
    ListedColormap,
    BoundaryNorm
)
from matplotlib.patches import Rectangle, Patch
from matplotlib import font_manager

# ============================================================
# FIGURE 3D
# COMPLEX ONE HEALTH AMR-CLASS HEATMAP
#
# INPUT: FOUR EXISTING FIGURE 3C EXCEL WORKBOOKS
#
# UPPER PANEL:
# Country-equal AMR class prevalence (%)
#
# LOWER PANEL:
# Mann-Whitney pairwise significance, from
# the SAME previously calculated FDR results
#
# 4 species x 5 compartments = 20 rows
# 10 antibiotic classes x 3 incomes = 30 columns
#
# Missing prevalence = grey
# No statistical test = NT
# Tested, non-significant = ns
#
# FIXED INCOME COLORS:
# H    = #C98A00
# M    = #673AB7
# L+LM = #C94F82
#
# OUTPUT:
# ONE 600-DPI JPEG
# ONE EXCEL WORKBOOK
#
# NO NEW STATISTICAL TESTS PERFORMED
# ============================================================

ROOT = Path(
    r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES"
)

INPUT_ROOT = ROOT / "Fig_03C_AMR_Classes"

OUT = ROOT / "Fig_03D_Complex_AMR_Heatmap"
OUT.mkdir(parents=True, exist_ok=True)

JPG = OUT / "Fig_03D_Complex_AMR_Class_Heatmap.jpg"
XLSX = OUT / "Fig_03D_Complex_AMR_Class_Heatmap_Data.xlsx"

SPECIES = [
    ("Enterobacter hormaechei", "E_hormaechei"),
    ("Enterobacter cloacae", "E_cloacae"),
    ("Enterobacter roggenkampii", "E_roggenkampii"),
    ("Enterobacter asburiae", "E_asburiae")
]

SOURCES = [
    "Human",
    "Animal",
    "Environment",
    "Food",
    "Plant-associated"
]

CLASSES = [
    "Aminoglycosides",
    "Quinolones",
    "Sulfonamides",
    "Trimethoprim",
    "Tetracyclines",
    "Phenicols",
    "Macrolides",
    "Colistin",
    "Beta-lactams",
    "Fosfomycin"
]

INCOMES = ["H", "M", "L+LM"]

PAIRS = [
    ("H", "M", "H–M"),
    ("M", "L+LM", "M–L"),
    ("H", "L+LM", "H–L")
]

INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82"
}

SPECIES_COLORS = [
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7"
]

SOURCE_COLORS = {
    "Human": "#E69F00",
    "Animal": "#56B4E9",
    "Environment": "#009E73",
    "Food": "#B07AA1",
    "Plant-associated": "#A6761D"
}

# ============================================================
# FONT
# ============================================================

available = {
    f.name
    for f in font_manager.fontManager.ttflist
}

font_name = next(
    (
        x for x in [
            "Arial",
            "Helvetica",
            "Liberation Sans",
            "DejaVu Sans"
        ]
        if x in available
    ),
    "DejaVu Sans"
)

plt.rcParams.update({
    "font.family": font_name,
    "font.size": 11,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})

# ============================================================
# UTILITIES
# ============================================================

def significance(q, status):

    if str(status).strip().upper() != "TESTED":
        return "NT"

    if pd.isna(q):
        raise ValueError(
            "A TESTED comparison has no FDR value"
        )

    q = float(q)

    if not 0 <= q <= 1:
        raise ValueError("FDR outside [0,1]")

    if q < 0.0001:
        return "****"
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"

    return "ns"


def get_one(table, filters):

    mask = pd.Series(
        True,
        index=table.index
    )

    for column, value in filters.items():
        mask &= table[column].eq(value)

    sub = table.loc[mask]

    if len(sub) != 1:
        raise ValueError(
            "Expected exactly one matching row: "
            + str(filters)
            + f"; found {len(sub)}"
        )

    return sub.iloc[0]


def add_group_grid(ax, nrows, ncols):

    # Vertical lines between antibiotic classes.
    for class_index in range(1, len(CLASSES)):
        ax.axvline(
            class_index * 3 - 0.5,
            color="white",
            linewidth=2.3,
            zorder=5
        )

    # Horizontal separators between species.
    for species_index in range(1, len(SPECIES)):
        ax.axhline(
            species_index * len(SOURCES) - 0.5,
            color="#383838",
            linewidth=1.8,
            zorder=5
        )

    ax.set_xlim(-0.5, ncols - 0.5)
    ax.set_ylim(nrows - 0.5, -0.5)

    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_top_headers(ax, labels, colors):

    # Individual three-column sublabels.
    for i, label in enumerate(labels):
        ax.text(
            i,
            -0.86,
            label,
            ha="center",
            va="bottom",
            fontsize=8.3,
            fontweight="bold",
            color=colors[i],
            clip_on=False
        )

    # Class names over each block of 3 columns.
    for class_index, class_name in enumerate(CLASSES):

        x = class_index * 3 + 1

        ax.text(
            x,
            -1.63,
            class_name,
            rotation=34,
            ha="left",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            color="#222222",
            clip_on=False
        )


def draw_side_labels(ax):

    for s_idx, (species, _) in enumerate(SPECIES):

        start = s_idx * len(SOURCES)

        # Species identity strip.
        ax.add_patch(
            Rectangle(
                (-2.18, start - 0.48),
                0.20,
                len(SOURCES) - 0.04,
                facecolor=SPECIES_COLORS[s_idx],
                edgecolor="none",
                clip_on=False
            )
        )

        short = species.replace(
            "Enterobacter ", "E. "
        )

        ax.text(
            -2.36,
            start + 2,
            short,
            ha="right",
            va="center",
            fontsize=11.5,
            fontstyle="italic",
            fontweight="bold",
            color=SPECIES_COLORS[s_idx],
            clip_on=False
        )

        for j, source in enumerate(SOURCES):

            y = start + j

            # Thin source-colour strip.
            ax.add_patch(
                Rectangle(
                    (-0.78, y - 0.43),
                    0.14,
                    0.86,
                    facecolor=SOURCE_COLORS[source],
                    edgecolor="none",
                    clip_on=False
                )
            )

            ax.text(
                -0.88,
                y,
                source,
                ha="right",
                va="center",
                fontsize=10,
                color="#333333",
                clip_on=False
            )


# ============================================================
# 1. LOAD ALL FOUR WORKBOOKS
# ============================================================

print("=" * 100)
print("FIGURE 3D - COMPLEX AMR CLASS HEATMAP")
print("=" * 100)

prevalence_parts = []
test_parts = []

for species, tag in SPECIES:

    workbook = (
        INPUT_ROOT
        / tag
        / f"Fig03C_{tag}_AMR_Classes_Data.xlsx"
    )

    if not workbook.is_file():
        raise FileNotFoundError(
            "Missing Figure 3C workbook:\n"
            + str(workbook)
        )

    prevalence = pd.read_excel(
        workbook,
        sheet_name="CLASS_PREVALENCE",
        engine="openpyxl"
    )

    tests = pd.read_excel(
        workbook,
        sheet_name="MANN_WHITNEY_FDR",
        engine="openpyxl"
    )

    for col in [
        "Species",
        "Source",
        "Income_Group",
        "Antibiotic_Class",
        "Country_Equal_Percent",
        "Observed_Genomes",
        "Eligible_Countries"
    ]:
        if col not in prevalence.columns:
            raise ValueError(
                f"Missing prevalence column: {col}"
            )

    for col in [
        "Species",
        "Source",
        "Antibiotic_Class",
        "Reference",
        "Comparison",
        "Raw_P",
        "FDR",
        "Status"
    ]:
        if col not in tests.columns:
            raise ValueError(
                f"Missing statistical column: {col}"
            )

    if not prevalence["Species"].eq(species).all():
        raise ValueError("Wrong species in prevalence")

    if not tests["Species"].eq(species).all():
        raise ValueError("Wrong species in statistics")

    if len(prevalence) != 150:
        raise ValueError(
            f"Expected 150 prevalence records for {species}"
        )

    if len(tests) != 150:
        raise ValueError(
            f"Expected 150 statistical records for {species}"
        )

    prevalence_parts.append(prevalence)
    test_parts.append(tests)

    print("[LOADED]", species)

prevalence = pd.concat(
    prevalence_parts,
    ignore_index=True
)

tests = pd.concat(
    test_parts,
    ignore_index=True
)

# ============================================================
# 2. VALIDATION
# ============================================================

if len(prevalence) != 600:
    raise ValueError("Expected 600 prevalence records")

if len(tests) != 600:
    raise ValueError("Expected 600 test records")

prev_keys = [
    "Species",
    "Source",
    "Income_Group",
    "Antibiotic_Class"
]

test_keys = [
    "Species",
    "Source",
    "Antibiotic_Class",
    "Reference",
    "Comparison"
]

if prevalence.duplicated(prev_keys).any():
    raise ValueError("Duplicate prevalence records")

if tests.duplicated(test_keys).any():
    raise ValueError("Duplicate statistical records")

prevalence["Country_Equal_Percent"] = pd.to_numeric(
    prevalence["Country_Equal_Percent"],
    errors="coerce"
)

tests["Raw_P"] = pd.to_numeric(
    tests["Raw_P"],
    errors="coerce"
)

tests["FDR"] = pd.to_numeric(
    tests["FDR"],
    errors="coerce"
)

valid_prev = prevalence[
    "Country_Equal_Percent"
].dropna()

if not valid_prev.between(0, 100).all():
    raise ValueError("Invalid prevalence percentages")

tests["Significance"] = [
    significance(q, status)
    for q, status in zip(
        tests["FDR"],
        tests["Status"]
    )
]

# ============================================================
# 3. BUILD ALIGNED HEATMAP MATRICES
# ============================================================

NROWS = len(SPECIES) * len(SOURCES)
NCOLS = len(CLASSES) * len(INCOMES)

PREVALENCE_MATRIX = np.full(
    (NROWS, NCOLS),
    np.nan
)

SIGNIFICANCE_MATRIX = np.full(
    (NROWS, NCOLS),
    0,
    dtype=int
)

SIGNIFICANCE_TEXT = np.full(
    (NROWS, NCOLS),
    "",
    dtype=object
)

PREVALENCE_N = np.zeros(
    (NROWS, NCOLS),
    dtype=int
)

COUNTRY_N = np.zeros(
    (NROWS, NCOLS),
    dtype=int
)

matrix_rows = []
stat_matrix_rows = []

for species_idx, (species, _) in enumerate(SPECIES):

    for source_idx, source in enumerate(SOURCES):

        row_index = (
            species_idx * len(SOURCES)
            + source_idx
        )

        for class_index, class_name in enumerate(CLASSES):

            for income_index, income in enumerate(INCOMES):

                col_index = class_index * 3 + income_index

                item = get_one(
                    prevalence,
                    {
                        "Species": species,
                        "Source": source,
                        "Antibiotic_Class": class_name,
                        "Income_Group": income
                    }
                )

                value = item["Country_Equal_Percent"]

                genomes = int(item["Observed_Genomes"])
                countries = int(item["Eligible_Countries"])

                if pd.notna(value):
                    PREVALENCE_MATRIX[
                        row_index, col_index
                    ] = float(value)

                PREVALENCE_N[
                    row_index, col_index
                ] = genomes

                COUNTRY_N[
                    row_index, col_index
                ] = countries

                matrix_rows.append({
                    "Species": species,
                    "Source": source,
                    "Antibiotic_Class": class_name,
                    "Income_Group": income,
                    "Country_Equal_Prevalence": value,
                    "Genomes": genomes,
                    "Eligible_Countries": countries
                })

            for pair_index, (
                ref, comp, pair_label
            ) in enumerate(PAIRS):

                col_index = class_index * 3 + pair_index

                test = get_one(
                    tests,
                    {
                        "Species": species,
                        "Source": source,
                        "Antibiotic_Class": class_name,
                        "Reference": ref,
                        "Comparison": comp
                    }
                )

                label = significance(
                    test["FDR"],
                    test["Status"]
                )

                SIGNIFICANCE_TEXT[
                    row_index, col_index
                ] = label

                if label == "NT":
                    code = 0
                elif label == "ns":
                    code = 1
                elif label == "*":
                    code = 2
                elif label == "**":
                    code = 3
                elif label == "***":
                    code = 4
                else:
                    code = 5

                SIGNIFICANCE_MATRIX[
                    row_index, col_index
                ] = code

                stat_matrix_rows.append({
                    "Species": species,
                    "Source": source,
                    "Antibiotic_Class": class_name,
                    "Reference": ref,
                    "Comparison": comp,
                    "Raw_P": test["Raw_P"],
                    "FDR": test["FDR"],
                    "Significance": label,
                    "Status": test["Status"]
                })

# ============================================================
# 4. COLOR PALETTES
# ============================================================

prevalence_cmap = LinearSegmentedColormap.from_list(
    "AMR_Prevalence",
    [
        "#F8FBFD",
        "#D3EBF1",
        "#91CAD8",
        "#3A9FB6",
        "#20708D",
        "#173B61"
    ]
)

prevalence_cmap.set_bad("#E0E0E0")

significance_cmap = ListedColormap([
    "#F0F0F0",  # NT
    "#FFFFFF",  # ns
    "#F8D7D7",  # *
    "#EB9E9E",  # **
    "#CD5151",  # ***
    "#962828"   # ****
])

significance_norm = BoundaryNorm(
    np.arange(-0.5, 6.5, 1),
    significance_cmap.N
)

# ============================================================
# 5. BUILD FIGURE
# ============================================================

print("\nRendering integrated complex heatmap...")

fig = plt.figure(
    figsize=(26, 21.5),
    facecolor="white"
)

gs = fig.add_gridspec(
    2,
    1,
    left=0.18,
    right=0.93,
    top=0.855,
    bottom=0.115,
    hspace=0.19
)

ax_prevalence = fig.add_subplot(gs[0])
ax_significance = fig.add_subplot(gs[1])

fig.suptitle(
    "Antibiotic resistance class repertoire across "
    "the Enterobacter cloacae complex",
    fontsize=22,
    fontweight="bold",
    y=0.983
)

fig.text(
    0.5,
    0.954,
    "Species-specific and One Health-stratified "
    "comparison across country-income groups",
    ha="center",
    fontsize=13,
    color="#555555"
)

# ============================================================
# UPPER HEATMAP: PREVALENCE
# ============================================================

masked = np.ma.masked_invalid(
    PREVALENCE_MATRIX
)

im = ax_prevalence.imshow(
    masked,
    cmap=prevalence_cmap,
    vmin=0,
    vmax=100,
    aspect="auto",
    interpolation="nearest"
)

add_group_grid(
    ax_prevalence,
    NROWS,
    NCOLS
)

draw_side_labels(ax_prevalence)

upper_labels = INCOMES * len(CLASSES)

upper_colors = [
    INCOME_COLORS[x] for x in upper_labels
]

draw_top_headers(
    ax_prevalence,
    upper_labels,
    upper_colors
)

ax_prevalence.set_title(
    "Annotated AMR class-carrier prevalence (%)",
    loc="left",
    fontsize=15,
    fontweight="bold",
    pad=105
)

# Label individual heatmap cells with percentages.
for i in range(NROWS):
    for j in range(NCOLS):

        value = PREVALENCE_MATRIX[i, j]

        if not np.isfinite(value):
            txt = "—"
            color = "#777777"
        else:
            txt = f"{value:.0f}"
            color = (
                "white"
                if value >= 67
                else "#202020"
            )

        ax_prevalence.text(
            j,
            i,
            txt,
            ha="center",
            va="center",
            fontsize=8.1,
            color=color
        )

# Prevalence colorbar.
cax = fig.add_axes([
    0.947,
    0.535,
    0.013,
    0.235
])

cbar = fig.colorbar(
    im,
    cax=cax
)

cbar.set_label(
    "Genomes carrying class determinant (%)",
    fontsize=11
)

cbar.set_ticks(
    [0, 20, 40, 60, 80, 100]
)

# ============================================================
# LOWER HEATMAP: STATISTICAL SIGNIFICANCE
# ============================================================

ax_significance.imshow(
    SIGNIFICANCE_MATRIX,
    cmap=significance_cmap,
    norm=significance_norm,
    aspect="auto",
    interpolation="nearest"
)

add_group_grid(
    ax_significance,
    NROWS,
    NCOLS
)

draw_side_labels(ax_significance)

lower_labels = [
    x[2]
    for _ in CLASSES
    for x in PAIRS
]

draw_top_headers(
    ax_significance,
    lower_labels,
    ["#444444"] * len(lower_labels)
)

ax_significance.set_title(
    "Pairwise income comparisons: "
    "FDR-adjusted Mann–Whitney U significance",
    loc="left",
    fontsize=15,
    fontweight="bold",
    pad=105
)

for i in range(NROWS):
    for j in range(NCOLS):

        label = SIGNIFICANCE_TEXT[i, j]

        if label in ["***", "****"]:
            text_color = "white"
        elif label == "NT":
            text_color = "#888888"
        elif label in ["*", "**"]:
            text_color = "#8C1616"
        else:
            text_color = "#444444"

        ax_significance.text(
            j,
            i,
            label,
            ha="center",
            va="center",
            fontsize=8.6,
            color=text_color,
            fontweight=(
                "bold"
                if "*" in label
                else "normal"
            )
        )

# ============================================================
# LEGENDS
# ============================================================

income_handles = [
    Patch(
        facecolor=INCOME_COLORS[x],
        edgecolor="none",
        label={
            "H": "High income",
            "M": "Upper-middle income",
            "L+LM": "Low + lower-middle income"
        }[x]
    )
    for x in INCOMES
]

source_handles = [
    Patch(
        facecolor=SOURCE_COLORS[x],
        edgecolor="none",
        label=x
    )
    for x in SOURCES
]

fig.legend(
    handles=income_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.082),
    ncol=3,
    frameon=False,
    fontsize=11
)

fig.legend(
    handles=source_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.052),
    ncol=5,
    frameon=False,
    fontsize=10
)

fig.text(
    0.5,
    0.026,
    "Significance: ns = tested but not significant; "
    "NT = insufficient countries; "
    "* q<0.05; ** q<0.01; "
    "*** q<0.001; **** q<0.0001.",
    ha="center",
    fontsize=10,
    color="#444444"
)

fig.text(
    0.5,
    0.009,
    "Class-carrier prevalence is based on exact "
    "AMRFinder class labels, not phenotypic drug resistance. "
    "Composite annotations are not split into separate classes.",
    ha="center",
    fontsize=9.4,
    color="#666666"
)

# ============================================================
# 6. SAVE JPEG
# ============================================================

fig.savefig(
    JPG,
    dpi=600,
    format="jpeg",
    facecolor="white",
    pil_kwargs={
        "quality": 94,
        "subsampling": 0,
        "optimize": True
    }
)

plt.close(fig)

# ============================================================
# 7. EXCEL SOURCE DATA
# ============================================================

prevalence_table = pd.DataFrame(
    matrix_rows
)

statistics_table = pd.DataFrame(
    stat_matrix_rows
)

prevalence_wide = prevalence_table.pivot(
    index=["Species", "Source"],
    columns=["Antibiotic_Class", "Income_Group"],
    values="Country_Equal_Prevalence"
)

statistics_wide = statistics_table.pivot(
    index=["Species", "Source"],
    columns=[
        "Antibiotic_Class",
        "Reference",
        "Comparison"
    ],
    values="Significance"
)

prevalence_columns = pd.MultiIndex.from_product(
    [CLASSES, INCOMES]
)

statistics_columns = pd.MultiIndex.from_tuples([
    (cls, ref, comp)
    for cls in CLASSES
    for ref, comp, _ in PAIRS
])

prevalence_wide = prevalence_wide.reindex(
    columns=prevalence_columns
)

statistics_wide = statistics_wide.reindex(
    columns=statistics_columns
)

methods = pd.DataFrame([
    (
        "Dataset",
        "Four ECC species, five One Health sources, "
        "three income groups, ten antibiotic classes"
    ),
    (
        "Input",
        "Existing Figure 3C statistical Excel workbooks"
    ),
    (
        "Primary prevalence outcome",
        "Country-equal mean percentage of genomes "
        "carrying the annotated AMR class"
    ),
    (
        "Missing prevalence",
        "Grey cells; no eligible country-level estimate"
    ),
    (
        "Significance test",
        "Existing two-sided Mann-Whitney U on "
        "country-level prevalence"
    ),
    (
        "FDR",
        "Existing Benjamini-Hochberg values; "
        "not recalculated"
    ),
    (
        "Significance NT",
        "Not tested due to insufficient eligible countries"
    ),
    (
        "Income colors",
        "H=#C98A00; M=#673AB7; L+LM=#C94F82"
    ),
    (
        "Class caveat",
        "Direct class annotations only; composite labels "
        "not split; not equivalent to phenotype"
    ),
    (
        "Created",
        datetime.now().isoformat()
    )
], columns=["Parameter", "Description"])

with pd.ExcelWriter(
    XLSX,
    engine="openpyxl"
) as writer:

    prevalence_table.to_excel(
        writer,
        sheet_name="PREVALENCE_LONG",
        index=False
    )

    statistics_table.to_excel(
        writer,
        sheet_name="SIGNIFICANCE_LONG",
        index=False
    )

    prevalence_wide.to_excel(
        writer,
        sheet_name="PREVALENCE_MATRIX"
    )

    statistics_wide.to_excel(
        writer,
        sheet_name="SIGNIFICANCE_MATRIX"
    )

    prevalence.to_excel(
        writer,
        sheet_name="FULL_PREVALENCE_SOURCE",
        index=False
    )

    tests.to_excel(
        writer,
        sheet_name="FULL_STATISTICS_SOURCE",
        index=False
    )

    methods.to_excel(
        writer,
        sheet_name="METHODS",
        index=False
    )

# ============================================================
# 8. REPORT
# ============================================================

print("\n" + "=" * 100)
print("FIGURE 3D - COMPLEX HEATMAP COMPLETED")
print("=" * 100)

print("\nSpecies:", len(SPECIES))
print("One Health compartments:", len(SOURCES))
print("Antibiotic classes:", len(CLASSES))

print("Prevalence matrix:", PREVALENCE_MATRIX.shape)
print("Significance matrix:", SIGNIFICANCE_MATRIX.shape)

print("\nSIGNIFICANCE COUNTS")

counts = tests["Significance"].value_counts()

for label in [
    "****",
    "***",
    "**",
    "*",
    "ns",
    "NT"
]:
    print(
        f"{label:>4}: {int(counts.get(label, 0))}"
    )

print("\nFDR-SIGNIFICANT RESULTS")

significant = tests[
    tests["Significance"].isin([
        "*", "**", "***", "****"
    ])
]

if significant.empty:

    print(
        "No FDR-significant comparisons "
        "in existing statistical workbooks."
    )

else:

    print(
        significant[
            [
                "Species",
                "Source",
                "Antibiotic_Class",
                "Reference",
                "Comparison",
                "Raw_P",
                "FDR",
                "Significance"
            ]
        ].round(6).to_string(index=False)
    )

print("\nOUTPUT JPEG:")
print(JPG)

print("\nOUTPUT EXCEL:")
print(XLSX)

print("\nJPEG DPI: 600")
print("Excel worksheets: 7")
print("Existing Figure 3C data unchanged")

print("=" * 100)
