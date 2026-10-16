from pathlib import Path
import gc
import shutil
import tempfile

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import font_manager

# ============================================================
# FIGURE 3G - FINAL VIOLIN PLOT VERSION
#
# Uses existing statistics from:
# Fig03G_FINAL_ALL_GENOMES_STATISTICS.xlsx
#
# Changes:
# - violin plot instead of boxplot
# - mean shown as black diamond + mean text
# - no country counts
# - only n and mean displayed
# - deep blue q-values on significance brackets
# - glossy/prominent source colors
# - soft income-specific panel backgrounds
# - safe replacement of previous Figure 3G images
#
# No statistical recalculation
# Excel workbook preserved
# ============================================================

ROOT = Path(
    r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES"
    r"\Fig_03G_OneHealth_Genome_Burden"
)

STAT_FILE = ROOT / "Fig03G_FINAL_ALL_GENOMES_STATISTICS.xlsx"
OUTPUT = ROOT / "FINAL_VIOLIN_MEAN_FDR_FIGURES"

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

INCOME_LIGHT = {
    "H": "#FCF6E9",
    "M": "#F4F0FA",
    "L+LM": "#FBEFF3"
}

INCOME_TITLES = {
    "H": "High income",
    "M": "Upper-middle income",
    "L+LM": "Low + lower-middle income"
}

# Fixed source colors across the project
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
DEEP_BLUE = "#0B3D91"

# ============================================================
# FONT / STYLE
# ============================================================

available_fonts = {
    f.name for f in font_manager.fontManager.ttflist
}

FONT = next(
    (
        f for f in [
            "Arial",
            "Liberation Sans",
            "DejaVu Sans"
        ]
        if f in available_fonts
    ),
    "DejaVu Sans"
)

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 11,
    "axes.titlesize": 14,
    "axes.labelsize": 11,
    "axes.linewidth": 0.8,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})

# ============================================================
# FUNCTIONS
# ============================================================

def short_species(name):
    return name.replace("Enterobacter ", "E. ")

def species_tag(name):
    return name.replace("Enterobacter ", "E_")

def fdr_text(value):
    if pd.isna(value):
        return ""
    value = float(value)
    if value < 0.0001:
        return f"q={value:.2e}"
    if value < 0.001:
        return f"q={value:.4f}"
    return f"q={value:.3f}"

def get_result(tests, species, income, a, b):
    result = tests.loc[
        tests["Species"].eq(species)
        & tests["Income_Group"].eq(income)
        & tests["Source_A"].eq(a)
        & tests["Source_B"].eq(b)
    ]
    if len(result) != 1:
        raise RuntimeError(
            f"Missing/duplicate result for "
            f"{species}, {income}, {a}, {b}"
        )
    return result.iloc[0]

def draw_violin(ax, position, values, source, rng):
    parts = ax.violinplot(
        [values],
        positions=[position],
        widths=0.82,
        showmeans=False,
        showmedians=False,
        showextrema=False
    )

    for body in parts["bodies"]:
        body.set_facecolor(SOURCE_COLORS[source])
        body.set_edgecolor("#333333")
        body.set_linewidth(1.05)
        body.set_alpha(0.92)

    # median line
    median = float(np.median(values))
    ax.hlines(
        median,
        position - 0.17,
        position + 0.17,
        color="#1E1E1E",
        linewidth=1.8,
        zorder=5
    )

    # mean point (black diamond)
    mean_value = float(np.mean(values))
    ax.scatter(
        position,
        mean_value,
        marker="D",
        s=34,
        color="#111111",
        edgecolor="white",
        linewidth=0.45,
        zorder=6
    )

    # scatter points (subsampled only for display)
    if len(values) > MAX_DISPLAY_POINTS:
        shown = rng.choice(
            values,
            size=MAX_DISPLAY_POINTS,
            replace=False
        )
    else:
        shown = values

    jitter = rng.uniform(
        -0.16, 0.16,
        size=len(shown)
    )

    ax.scatter(
        np.full(len(shown), position) + jitter,
        shown,
        s=8,
        color=SOURCE_COLORS[source],
        alpha=0.24,
        edgecolor="none",
        zorder=4
    )

    return mean_value, median

def add_stat_bracket(ax, x1, x2, y_fraction, label, qvalue):
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
        color="#3F3F3F",
        linewidth=0.9,
        clip_on=False,
        zorder=10
    )

    ax.text(
        (x1 + x2) / 2,
        y_fraction + 0.014,
        f"{label}  ({fdr_text(qvalue)})",
        transform=transform,
        ha="center",
        va="bottom",
        fontsize=9.15,
        fontweight="bold" if "*" in str(label) else "normal",
        color=DEEP_BLUE,
        clip_on=False,
        zorder=11
    )

# ============================================================
# LOAD EXISTING RESULTS
# ============================================================

print("=" * 105)
print("FIGURE 3G - FINAL VIOLIN REFINEMENT")
print("=" * 105)

if not STAT_FILE.is_file():
    raise FileNotFoundError(
        "Statistics workbook not found:\n" + str(STAT_FILE)
    )

tests = pd.read_excel(
    STAT_FILE,
    sheet_name="ALL_63_COMPARISONS",
    engine="openpyxl"
)

genomes = pd.read_excel(
    STAT_FILE,
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
        raise ValueError("Missing test column: " + column)

for column in required_genomes:
    if column not in genomes.columns:
        raise ValueError("Missing genome column: " + column)

if len(tests) != 63:
    raise ValueError(f"Expected 63 comparisons, found {len(tests)}")

if tests[[
    "Species", "Income_Group", "Source_A", "Source_B"
]].duplicated().any():
    raise ValueError("Duplicate test entries found")

if genomes["Assembly_Accession"].duplicated().any():
    raise ValueError("Duplicate genome accessions found")

genomes["AMR"] = pd.to_numeric(
    genomes["AMR"],
    errors="coerce"
)

tests["FDR"] = pd.to_numeric(
    tests["FDR"],
    errors="coerce"
)

if (
    genomes["AMR"].isna().any()
    or (genomes["AMR"] < 0).any()
):
    raise ValueError("Invalid genome AMR values")

# verify consistency between workbook counts and data
count_lookup = genomes.groupby(
    ["Species", "Income_Group", "Source"]
).size().to_dict()

for rec in tests.itertuples():
    actual_a = count_lookup.get(
        (rec.Species, rec.Income_Group, rec.Source_A),
        0
    )
    actual_b = count_lookup.get(
        (rec.Species, rec.Income_Group, rec.Source_B),
        0
    )
    if actual_a != rec.N_A or actual_b != rec.N_B:
        raise ValueError(
            "Genome counts in workbook do not match "
            "genome-level data"
        )

# ============================================================
# SAFE TEMPORARY GENERATION
# ============================================================

OUTPUT.mkdir(parents=True, exist_ok=True)

temp_dir = Path(
    tempfile.mkdtemp(
        prefix="Fig03G_violin_",
        dir=str(ROOT)
    )
)

print("Temporary directory:", temp_dir)

generated_temp = []

try:

    for species_index, species in enumerate(SPECIES):

        print("\nProcessing:", species, flush=True)

        rng = np.random.default_rng(SEED + species_index)

        species_data = genomes.loc[
            genomes["Species"].eq(species)
        ].copy()

        if species_data.empty:
            raise RuntimeError(
                f"No genome data for {species}"
            )

        ymax = max(
            6.0,
            float(species_data["AMR"].max()) * 1.10 + 1.2
        )

        fig, axes = plt.subplots(
            1, 3,
            figsize=(19.4, 9.1),
            sharey=True,
            facecolor="white"
        )

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
            "One Health variation in annotated AMR gene burden per genome",
            ha="center",
            va="center",
            fontsize=12.2,
            color="#4C4C4C"
        )

        for panel_index, income in enumerate(INCOMES):

            ax = axes[panel_index]

            panel = species_data.loc[
                species_data["Income_Group"].eq(income)
            ]

            ax.set_facecolor(INCOME_LIGHT[income])

            # panel title
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
                [0.08, 0.92],
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

                values = group["AMR"].to_numpy(dtype=float)
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
                        color="#9B9B9B"
                    )
                    continue

                mean_value, median_value = draw_violin(
                    ax, x, values, source, rng
                )

                # only n and mean
                ax.text(
                    x,
                    0.968,
                    f"n={len(values):,}\nMean={mean_value:.2f}",
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    va="top",
                    fontsize=9.6,
                    fontweight="medium",
                    color="#303030"
                )

            # significance brackets
            bracket_positions = [
                ("Human", "Animal", 1.085),
                ("Animal", "Environment", 1.170),
                ("Human", "Environment", 1.265)
            ]

            for source_a, source_b, y_frac in bracket_positions:

                result = get_result(
                    tests, species, income, source_a, source_b
                )

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
                color="#D9D9D9",
                alpha=0.60,
                linewidth=0.65
            )

            ax.set_axisbelow(True)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_color("#999999")
            ax.spines["bottom"].set_color("#999999")

            if panel_index == 0:
                ax.set_ylabel(
                    "Unique annotated AMR genes per genome",
                    fontsize=12,
                    labelpad=10
                )

        # legend
        handles = [
            Patch(
                facecolor=SOURCE_COLORS[source],
                edgecolor="#444444",
                alpha=0.85,
                label=source
            )
            for source in SOURCES
        ]

        fig.legend(
            handles=handles,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.128),
            ncol=3,
            frameon=False,
            fontsize=11
        )

        fig.text(
            0.5,
            0.080,
            "Violin plots show the distribution of genome-level AMR burden; "
            "black diamonds indicate the mean. All genomes contribute to statistics.",
            ha="center",
            fontsize=10,
            color="#555555"
        )

        fig.text(
            0.5,
            0.055,
            "Pairwise two-sided Mann–Whitney U tests; q = BH-adjusted P-value. "
            "ns = not significant; * q<0.05; ** q<0.01; *** q<0.001; **** q<0.0001.",
            ha="center",
            fontsize=9.8,
            color=DEEP_BLUE
        )

        fig.text(
            0.5,
            0.031,
            "Genome-level inference is exploratory and unadjusted for country, "
            "BioProject, collection year and genomic lineage.",
            ha="center",
            fontsize=9.2,
            color="#767676"
        )

        temp_file = temp_dir / (
            f"Fig03G_{species_tag(species)}_VIOLIN_MEAN_FDR.jpg"
        )

        fig.savefig(
            temp_file,
            dpi=DPI,
            format="jpeg",
            facecolor="white",
            pil_kwargs={
                "quality": 95,
                "subsampling": 0,
                "optimize": True
            }
        )

        plt.close(fig)
        gc.collect()

        if (not temp_file.is_file()) or temp_file.stat().st_size == 0:
            raise RuntimeError(
                "Failed to create figure: " + str(temp_file)
            )

        generated_temp.append(temp_file)

        significant_n = int(
            tests.loc[
                tests["Species"].eq(species)
                & (tests["Status"].eq("TESTED"))
                & (tests["FDR"] < 0.05)
            ].shape[0]
        )

        print(
            f"[TEMP SAVED] {temp_file.name} | "
            f"significant comparisons: {significant_n}",
            flush=True
        )

    # ========================================================
    # VALIDATE 7 FIGURES CREATED
    # ========================================================

    if len(generated_temp) != 7:
        raise RuntimeError(
            f"Expected 7 temporary figures, found {len(generated_temp)}"
        )

    # copy final files
    final_files = []

    for temp_file in generated_temp:
        final_file = OUTPUT / temp_file.name
        shutil.copy2(temp_file, final_file)

        if (not final_file.is_file()) or final_file.stat().st_size == 0:
            raise RuntimeError(
                "Failed to copy final figure: " + str(final_file)
            )

        final_files.append(final_file)

    print("\nAll 7 new figures copied successfully.")

    # ========================================================
    # DELETE OLD FIGURE 3G IMAGE FILES
    # ========================================================

    final_set = {f.resolve() for f in final_files}
    output_resolved = OUTPUT.resolve()
    temp_resolved = temp_dir.resolve()

    removed = 0

    for path in ROOT.rglob("*"):

        if not path.is_file():
            continue

        if path.suffix.lower() not in {".jpg", ".jpeg", ".png"}:
            continue

        resolved = path.resolve()

        # keep newly created final images
        if resolved in final_set:
            continue

        # keep temporary files until cleanup
        if temp_resolved in resolved.parents:
            continue

        # delete previous figure images anywhere
        # inside Figure 3G directory
        try:
            path.unlink()
            removed += 1
        except Exception as exc:
            raise RuntimeError(
                f"Could not delete old figure file:\n{path}\n{exc}"
            )

    print(f"Deleted old Figure 3G image files: {removed}")

finally:
    try:
        shutil.rmtree(temp_dir, ignore_errors=True)
    except Exception:
        pass

# ============================================================
# FINAL REPORT
# ============================================================

print("\n" + "=" * 105)
print("FIGURE 3G - FINAL VIOLIN FIGURES COMPLETED")
print("=" * 105)

print("Species figures created:", 7)
print("Resolution:", DPI, "DPI")
print("Input statistics:", STAT_FILE)
print("Final output folder:", OUTPUT)
print("No statistics recalculated.")
print("Excel workbook preserved.")
print("=" * 105)