
from pathlib import Path
import gc

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import font_manager

# ============================================================
# FIGURE 3G - FINAL PUBLICATION-STYLE FIGURES
#
# USE EXISTING FINAL STATISTICAL RESULTS
#
# 7 species, 3 income panels per species
# Human, Animal, Environment
#
# Changes:
# - No country-count labels
# - Display genome sample sizes only
# - Exact BH-FDR q-values on each graph
# - Significance brackets above boxplot area
# - Soft income-specific colored backgrounds
# - Improved titles and spacing
#
# No statistical recalculation
# No previous files deleted
#
# Output: 7 new JPEG figures at 600 DPI
# ============================================================

ROOT = Path(
    r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES"
    r"\Fig_03G_OneHealth_Genome_Burden"
)

INPUT_EXCEL = (
    ROOT / "Fig03G_FINAL_ALL_GENOMES_STATISTICS.xlsx"
)

OUTPUT = ROOT / "FINAL_BEAUTIFIED_FIGURES"
OUTPUT.mkdir(parents=True, exist_ok=True)

SPECIES = [
    "Enterobacter hormaechei",
    "Enterobacter cloacae",
    "Enterobacter roggenkampii",
    "Enterobacter asburiae",
    "Enterobacter kobei",
    "Enterobacter ludwigii",
    "Enterobacter bugandensis"
]

INCOMES = ["H", "M", "L+LM"]

SOURCES = ["Human", "Animal", "Environment"]

PAIRS = [
    ("Human", "Animal"),
    ("Animal", "Environment"),
    ("Human", "Environment")
]

INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82"
}

# Light backgrounds derived from fixed income colors.
INCOME_BACKGROUND = {
    "H": "#FCF6E9",
    "M": "#F4F0FA",
    "L+LM": "#FBEFF3"
}

INCOME_TITLES = {
    "H": "High income",
    "M": "Upper-middle income",
    "L+LM": "Low + lower-middle income"
}

SOURCE_COLORS = {
    "Human": "#E69F00",
    "Animal": "#56B4E9",
    "Environment": "#009E73"
}

POSITIONS = {
    "Human": 1,
    "Animal": 2,
    "Environment": 3
}

MAX_DISPLAY_POINTS = 350
DPI = 600
SEED = 20261010

# ============================================================
# FONT AND STYLE
# ============================================================

available = {
    f.name
    for f in font_manager.fontManager.ttflist
}

FONT = next(
    (
        name
        for name in [
            "Arial",
            "Liberation Sans",
            "DejaVu Sans"
        ]
        if name in available
    ),
    "DejaVu Sans"
)

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 11,
    "axes.titlesize": 14,
    "axes.labelsize": 11,
    "axes.linewidth": 0.75,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})

# ============================================================
# FUNCTIONS
# ============================================================

def short_species(name):
    return name.replace(
        "Enterobacter ", "E. "
    )


def tag_species(name):
    return name.replace(
        "Enterobacter ", "E_"
    )


def fdr_text(value):

    if pd.isna(value):
        return ""

    value = float(value)

    if value < 0.0001:
        return f"q={value:.2e}"

    if value < 0.001:
        return f"q={value:.4f}"

    return f"q={value:.3f}"


def get_statistics(tests, species, income, a, b):

    result = tests.loc[
        tests["Species"].eq(species)
        & tests["Income_Group"].eq(income)
        & tests["Source_A"].eq(a)
        & tests["Source_B"].eq(b)
    ]

    if len(result) != 1:
        raise RuntimeError(
            f"Missing/duplicate result: "
            f"{species}, {income}, {a}, {b}"
        )

    return result.iloc[0]


def plot_box(ax, position, values, source, rng):

    bp = ax.boxplot(
        [values],
        positions=[position],
        widths=0.51,
        patch_artist=True,
        showfliers=False,
        whis=1.5,
        medianprops={
            "color": "#202020",
            "linewidth": 1.8
        },
        boxprops={
            "edgecolor": "#474747",
            "linewidth": 0.95
        },
        whiskerprops={
            "color": "#555555",
            "linewidth": 0.95
        },
        capprops={
            "color": "#555555",
            "linewidth": 0.95
        }
    )

    bp["boxes"][0].set_facecolor(
        SOURCE_COLORS[source]
    )

    bp["boxes"][0].set_alpha(0.65)

    if len(values) > MAX_DISPLAY_POINTS:

        displayed = rng.choice(
            values,
            size=MAX_DISPLAY_POINTS,
            replace=False
        )

    else:
        displayed = values

    jitter = rng.uniform(
        -0.17,
        0.17,
        size=len(displayed)
    )

    ax.scatter(
        position + jitter,
        displayed,
        s=8,
        color=SOURCE_COLORS[source],
        alpha=0.30,
        edgecolor="none",
        zorder=3
    )


def add_stat_bracket(
    ax,
    x1,
    x2,
    y_fraction,
    label,
    qvalue
):

    # Coordinates: x in data units, y in axes fraction.
    # Brackets occupy a dedicated area above the axes,
    # rather than sharing space with genome observations.

    transform = ax.get_xaxis_transform()

    ax.plot(
        [x1, x1, x2, x2],
        [
            y_fraction - 0.013,
            y_fraction,
            y_fraction,
            y_fraction - 0.013
        ],
        transform=transform,
        color="#414141",
        linewidth=0.85,
        clip_on=False,
        zorder=10
    )

    is_significant = "*" in str(label)

    color = (
        "#A72734"
        if is_significant
        else "#4B4B4B"
    )

    ax.text(
        (x1 + x2) / 2,
        y_fraction + 0.013,
        f"{label}  ({fdr_text(qvalue)})",
        transform=transform,
        ha="center",
        va="bottom",
        fontsize=9.1,
        fontweight=(
            "bold" if is_significant
            else "normal"
        ),
        color=color,
        clip_on=False,
        zorder=11
    )


# ============================================================
# LOAD EXISTING FINAL STATISTICS
# ============================================================

print("=" * 105)
print("FIGURE 3G - PUBLICATION STYLE REFINEMENT")
print("=" * 105)

if not INPUT_EXCEL.is_file():
    raise FileNotFoundError(
        "Run the previous final analysis first:\n"
        + str(INPUT_EXCEL)
    )

tests = pd.read_excel(
    INPUT_EXCEL,
    sheet_name="ALL_63_COMPARISONS",
    engine="openpyxl"
)

genomes = pd.read_excel(
    INPUT_EXCEL,
    sheet_name="GENOME_AMR_DATA",
    engine="openpyxl"
)

required_tests = [
    "Species",
    "Income_Group",
    "Source_A",
    "Source_B",
    "N_A",
    "N_B",
    "Raw_P",
    "FDR",
    "Significance",
    "Status"
]

required_genomes = [
    "Assembly_Accession",
    "Species",
    "Income_Group",
    "Source",
    "Country",
    "AMR"
]

for column in required_tests:
    if column not in tests.columns:
        raise ValueError(
            "Missing test column: " + column
        )

for column in required_genomes:
    if column not in genomes.columns:
        raise ValueError(
            "Missing genome column: " + column
        )

if len(tests) != 63:
    raise ValueError(
        f"Expected 63 comparisons, found {len(tests)}"
    )

if tests[[
    "Species",
    "Income_Group",
    "Source_A",
    "Source_B"
]].duplicated().any():
    raise ValueError(
        "Duplicate statistical comparisons"
    )

if genomes["Assembly_Accession"].duplicated().any():
    raise ValueError(
        "Duplicate genome accessions"
    )

genomes["AMR"] = pd.to_numeric(
    genomes["AMR"],
    errors="coerce"
)

tests["FDR"] = pd.to_numeric(
    tests["FDR"],
    errors="coerce"
)

tests["Raw_P"] = pd.to_numeric(
    tests["Raw_P"],
    errors="coerce"
)

if (
    genomes["AMR"].isna().any()
    or (genomes["AMR"] < 0).any()
):
    raise ValueError(
        "Invalid genome AMR values"
    )

if (
    tests.loc[
        tests["Status"].eq("TESTED"),
        "FDR"
    ].isna().any()
):
    raise ValueError(
        "Missing FDR in tested comparisons"
    )

# Confirm no inconsistency between workbook sample
# counts and the actual genome-level data.

count_lookup = genomes.groupby(
    ["Species", "Income_Group", "Source"]
).size().to_dict()

for record in tests.itertuples():

    actual_a = count_lookup.get(
        (
            record.Species,
            record.Income_Group,
            record.Source_A
        ),
        0
    )

    actual_b = count_lookup.get(
        (
            record.Species,
            record.Income_Group,
            record.Source_B
        ),
        0
    )

    if (
        actual_a != record.N_A
        or actual_b != record.N_B
    ):
        raise ValueError(
            "Genome sample sizes differ from "
            "the statistical workbook"
        )

# ============================================================
# GENERATE SEVEN REFINED FIGURES
# ============================================================

generated = []

for species_index, species in enumerate(SPECIES):

    print("\nProcessing:", species, flush=True)

    species_data = genomes.loc[
        genomes["Species"].eq(species)
    ]

    rng = np.random.default_rng(
        SEED + species_index
    )

    # Common y-scale across the three income panels
    # within this species.

    maximum = float(
        species_data["AMR"].max()
    )

    ymax = max(
        6.0,
        maximum * 1.10 + 1.0
    )

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(19.2, 8.9),
        sharey=True,
        facecolor="white"
    )

    # Leave a large dedicated area for panel titles
    # and all statistical comparison brackets.
    fig.subplots_adjust(
        left=0.075,
        right=0.975,
        bottom=0.205,
        top=0.705,
        wspace=0.20
    )

    fig.suptitle(
        short_species(species),
        fontsize=23,
        fontweight="bold",
        fontstyle="italic",
        y=0.985
    )

    fig.text(
        0.5,
        0.947,
        "One Health variation in annotated "
        "AMR gene burden per genome",
        ha="center",
        va="center",
        fontsize=12.2,
        color="#505050"
    )

    for panel_index, income in enumerate(INCOMES):

        ax = axes[panel_index]

        panel = species_data.loc[
            species_data["Income_Group"].eq(income)
        ]

        # Income-specific pale background,
        # deliberately low contrast.
        ax.set_facecolor(
            INCOME_BACKGROUND[income]
        )

        # Upper strip uses the established
        # original income color.
        ax.text(
            0.5,
            1.395,
            INCOME_TITLES[income],
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=14,
            fontweight="bold",
            color=INCOME_COLORS[income],
            clip_on=False
        )

        ax.plot(
            [0, 1],
            [1.345, 1.345],
            transform=ax.transAxes,
            color=INCOME_COLORS[income],
            linewidth=3.4,
            solid_capstyle="round",
            clip_on=False
        )

        for source in SOURCES:

            group = panel.loc[
                panel["Source"].eq(source)
            ]

            values = group[
                "AMR"
            ].to_numpy(dtype=float)

            x = POSITIONS[source]

            if len(values) == 0:

                ax.text(
                    x,
                    ymax / 2,
                    "No genomes",
                    ha="center",
                    va="center",
                    rotation=90,
                    fontsize=10,
                    color="#999999"
                )

                continue

            plot_box(
                ax,
                x,
                values,
                source,
                rng
            )

            # Only genome sample size.
            # No country count anywhere on figure.
            ax.text(
                x,
                0.965,
                f"n={len(values):,}",
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=9.8,
                fontweight="medium",
                color="#343434"
            )

        # Three significance brackets, all located
        # above the genomic distribution panel.

        bracket_positions = [
            ("Human", "Animal", 1.085),
            ("Animal", "Environment", 1.165),
            ("Human", "Environment", 1.265)
        ]

        for source_a, source_b, y_frac in bracket_positions:

            result = get_statistics(
                tests,
                species,
                income,
                source_a,
                source_b
            )

            # Empty comparisons are not annotated.
            if result["Status"] != "TESTED":
                continue

            add_stat_bracket(
                ax,
                POSITIONS[source_a],
                POSITIONS[source_b],
                y_frac,
                result["Significance"],
                result["FDR"]
            )

        ax.set_xlim(0.4, 3.6)
        ax.set_ylim(0, ymax)

        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(
            SOURCES,
            fontsize=11,
            fontweight="medium"
        )

        ax.tick_params(
            axis="x",
            length=0,
            pad=11
        )

        ax.tick_params(
            axis="y",
            labelsize=10
        )

        ax.grid(
            axis="y",
            color="#D8D8D8",
            alpha=0.6,
            linewidth=0.65
        )

        ax.set_axisbelow(True)

        for spine in [
            "top",
            "right"
        ]:
            ax.spines[spine].set_visible(False)

        ax.spines["left"].set_color(
            "#999999"
        )

        ax.spines["bottom"].set_color(
            "#999999"
        )

        if panel_index == 0:
            ax.set_ylabel(
                "Unique annotated AMR genes per genome",
                fontsize=12,
                labelpad=10
            )

    # ========================================================
    # LEGEND
    # ========================================================

    handles = [
        Patch(
            facecolor=SOURCE_COLORS[source],
            edgecolor="#555555",
            alpha=0.70,
            label=source
        )
        for source in SOURCES
    ]

    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.130),
        ncol=3,
        frameon=False,
        fontsize=11
    )

    # ========================================================
    # FIGURE NOTES
    # ========================================================

    fig.text(
        0.5,
        0.080,
        "Boxplots: median and interquartile range; "
        "whiskers: 1.5 × IQR. "
        "All genomes contribute to calculations.",
        ha="center",
        fontsize=10,
        color="#555555"
    )

    fig.text(
        0.5,
        0.055,
        "Pairwise two-sided Mann–Whitney U tests; "
        "q = BH-adjusted P-value. "
        "ns = not significant; "
        "* q<0.05; ** q<0.01; "
        "*** q<0.001; **** q<0.0001.",
        ha="center",
        fontsize=9.8,
        color="#444444"
    )

    fig.text(
        0.5,
        0.030,
        "Genome-level inference is exploratory "
        "and unadjusted for country, BioProject, "
        "collection year and genomic lineage.",
        ha="center",
        fontsize=9.2,
        color="#777777"
    )

    # ========================================================
    # SAVE
    # ========================================================

    filename = (
        OUTPUT
        / (
            f"Fig03G_{tag_species(species)}"
            f"_BEAUTIFIED_WITH_FDR.jpg"
        )
    )

    fig.savefig(
        filename,
        dpi=DPI,
        format="jpeg",
        facecolor="white",
        pil_kwargs={
            "quality": 94,
            "subsampling": 0,
            "optimize": True
        }
    )

    plt.close(fig)
    gc.collect()

    if not filename.is_file():
        raise RuntimeError(
            "Failed to save: " + str(filename)
        )

    generated.append(filename)

    significant = tests.loc[
        tests["Species"].eq(species)
        & (tests["FDR"] < 0.05)
    ]

    print(
        f"[SAVED] {filename.name}"
        f" | significant comparisons: {len(significant)}",
        flush=True
    )

# ============================================================
# FINAL REPORT
# ============================================================

if len(generated) != 7:
    raise RuntimeError(
        "Expected seven final species figures"
    )

print("\n" + "=" * 105)
print("FIGURE 3G - FINAL BEAUTIFIED FIGURES COMPLETED")
print("=" * 105)

print("Species:", len(SPECIES))
print("JPEG figures:", len(generated))
print("JPEG resolution:", DPI, "DPI")
print("Source genomes:", len(genomes))

print(
    "Statistical comparisons tested:",
    int(tests["Status"].eq("TESTED").sum())
)

print(
    "FDR-significant comparisons:",
    int((tests["FDR"] < 0.05).sum())
)

print("\nOUTPUT DIRECTORY:")
print(OUTPUT)

print("\nStatistics read from:")
print(INPUT_EXCEL)

print("\nNo P-values recalculated.")
print("No original figures or workbooks modified.")
print("=" * 105)

