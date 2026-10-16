
from pathlib import Path
import gc
import shutil
import tempfile

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from matplotlib.patches import (
    FancyBboxPatch,
    Patch
)
from matplotlib.colors import to_rgb
from PIL import Image

# ============================================================
# FIGURE 3G — FINAL PUBLICATION VIOLIN FIGURES
#
# Uses the EXISTING statistics and genome-level values.
# Does not recalculate Mann-Whitney tests or FDR.
#
# Seven species; three income groups per species.
# Human, Animal, Environment source comparisons.
#
# Full species name, glossy violins, smaller widths,
# rounded framed panels, mean diamond + value,
# sample size n, no country counts.
#
# Statistical annotations:
# nonsignificant = deep blue
# significant stars + q = magenta
#
# Output: 7 JPEG figures at 600 DPI.
#
# Old Figure 3G images deleted ONLY AFTER all
# replacement images have been verified.
# Excel workbooks and Python scripts preserved.
# ============================================================

ROOT = Path(
    r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES"
    r"\Fig_03G_OneHealth_Genome_Burden"
)

INPUT = ROOT / "Fig03G_FINAL_ALL_GENOMES_STATISTICS.xlsx"

OUT = ROOT / "FINAL_PUBLICATION_VIOLINS"

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

# Income-group identity colors: UNCHANGED
INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82"
}

INCOME_TINTS = {
    "H": "#FFF8E9",
    "M": "#F5F1FD",
    "L+LM": "#FFF0F5"
}

INCOME_TITLES = {
    "H": "High income",
    "M": "Upper-middle income",
    "L+LM": "Low + lower-middle income"
}

# Updated One Health palette for this final figure
SOURCE_COLORS = {
    "Human": "#E8792B",
    "Animal": "#2896CF",
    "Environment": "#139C83"
}

POSITIONS = {
    "Human": 1,
    "Animal": 2,
    "Environment": 3
}

DEEP_BLUE = "#153E82"
SIGNIFICANT_COLOR = "#AD1457"
DARK = "#263449"

DPI = 600
MAX_POINTS = 350

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10.5,
    "axes.labelsize": 11,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})

# ============================================================
# FUNCTIONS
# ============================================================

def display_q(q):
    if pd.isna(q):
        return ""

    q = float(q)

    if q < 0.0001:
        return f"q={q:.2e}"
    if q < 0.001:
        return f"q={q:.4f}"
    return f"q={q:.3f}"


def lighten(color, fraction):
    rgb = np.asarray(to_rgb(color))
    return tuple(
        rgb + (1.0 - rgb) * fraction
    )


def darken(color, fraction):
    rgb = np.asarray(to_rgb(color))
    return tuple(rgb * (1.0 - fraction))


def get_test(tests, species, income, a, b):
    matches = tests.loc[
        tests["Species"].eq(species)
        & tests["Income_Group"].eq(income)
        & tests["Source_A"].eq(a)
        & tests["Source_B"].eq(b)
    ]

    if len(matches) != 1:
        raise ValueError(
            "Missing or duplicated comparison: "
            f"{species} | {income} | {a} | {b}"
        )

    return matches.iloc[0]


def draw_violin(ax, values, x, source, rng):
    color = SOURCE_COLORS[source]

    # KDE is undefined with one observation or
    # no variation. Display those groups with points.
    if len(values) >= 2 and np.ptp(values) > 0:

        violin = ax.violinplot(
            [values],
            positions=[x],
            widths=0.60,
            showmeans=False,
            showmedians=False,
            showextrema=False,
            bw_method="scott",
            points=120
        )

        for body in violin["bodies"]:
            body.set_facecolor(color)
            body.set_alpha(0.88)
            body.set_edgecolor(
                darken(color, 0.35)
            )
            body.set_linewidth(1.05)
            body.set_zorder(3)

            # Subtle glossy central highlight.
            # This does not change the distribution.
            ax.plot(
                [x - 0.065, x - 0.065],
                [
                    np.quantile(values, .25),
                    np.quantile(values, .75)
                ],
                color=lighten(color, .70),
                alpha=.55,
                linewidth=2.2,
                solid_capstyle="round",
                zorder=4
            )

    # Points are subsampled only for display.
    if len(values) > MAX_POINTS:
        shown = rng.choice(
            values,
            MAX_POINTS,
            replace=False
        )
    else:
        shown = values

    offsets = rng.uniform(
        -0.115, .115,
        size=len(shown)
    )

    ax.scatter(
        x + offsets,
        shown,
        s=9,
        facecolor=lighten(color, .34),
        edgecolor="none",
        alpha=.42,
        zorder=5
    )

    median = float(np.median(values))
    mean = float(np.mean(values))

    # Median marker
    ax.hlines(
        median,
        x - .15,
        x + .15,
        color="#243042",
        linewidth=1.55,
        zorder=7
    )

    # Mean diamond
    ax.scatter(
        [x],
        [mean],
        marker="D",
        s=53,
        facecolor="#171D27",
        edgecolor="white",
        linewidth=.9,
        zorder=9
    )

    return mean


def draw_bracket(ax, a, b, y, label, q):
    if label is None:
        return

    significant = "*" in str(label)

    text_color = (
        SIGNIFICANT_COLOR
        if significant
        else DEEP_BLUE
    )

    transform = ax.get_xaxis_transform()

    x1 = POSITIONS[a]
    x2 = POSITIONS[b]

    ax.plot(
        [x1, x1, x2, x2],
        [y-.014, y, y, y-.014],
        transform=transform,
        color=text_color if significant else "#566071",
        linewidth=1.00,
        clip_on=False,
        zorder=10
    )

    ax.text(
        (x1 + x2) / 2,
        y + .012,
        f"{label} ({display_q(q)})",
        transform=transform,
        ha="center",
        va="bottom",
        fontsize=9.0,
        weight="bold" if significant else "medium",
        color=text_color,
        clip_on=False,
        zorder=11
    )


# ============================================================
# INPUT VALIDATION
# ============================================================

print("=" * 90)
print("FIGURE 3G — PUBLICATION VIOLIN PLOTS")
print("=" * 90)

if not INPUT.is_file():
    raise FileNotFoundError(INPUT)

tests = pd.read_excel(
    INPUT,
    sheet_name="ALL_63_COMPARISONS",
    engine="openpyxl"
)

genomes = pd.read_excel(
    INPUT,
    sheet_name="GENOME_AMR_DATA",
    engine="openpyxl"
)

test_columns = [
    "Species", "Income_Group",
    "Source_A", "Source_B",
    "N_A", "N_B", "FDR",
    "Significance", "Status"
]

genome_columns = [
    "Assembly_Accession", "Species",
    "Income_Group", "Source", "AMR"
]

for c in test_columns:
    if c not in tests:
        raise ValueError(f"Missing test column: {c}")

for c in genome_columns:
    if c not in genomes:
        raise ValueError(f"Missing genome column: {c}")

if len(tests) != 63:
    raise ValueError(
        f"Expected 63 comparisons; got {len(tests)}"
    )

if tests[
    ["Species", "Income_Group", "Source_A", "Source_B"]
].duplicated().any():
    raise ValueError("Duplicated statistical comparisons")

if genomes["Assembly_Accession"].duplicated().any():
    raise ValueError("Duplicated genome accessions")

genomes["AMR"] = pd.to_numeric(
    genomes["AMR"],
    errors="coerce"
)

tests["FDR"] = pd.to_numeric(
    tests["FDR"],
    errors="coerce"
)

if genomes["AMR"].isna().any():
    raise ValueError("Missing AMR burden values")

if not genomes["Species"].isin(SPECIES).all():
    raise ValueError("Unexpected species")

if not genomes["Income_Group"].isin(INCOMES).all():
    raise ValueError("Unexpected income category")

if not genomes["Source"].isin(SOURCES).all():
    raise ValueError("Unexpected source")

if tests.loc[
    tests["Status"].eq("TESTED"),
    "FDR"
].isna().any():
    raise ValueError("Tested comparisons missing FDR")

counts = genomes.groupby(
    ["Species", "Income_Group", "Source"]
).size().to_dict()

for record in tests.itertuples():
    n_a = counts.get(
        (record.Species, record.Income_Group, record.Source_A),
        0
    )
    n_b = counts.get(
        (record.Species, record.Income_Group, record.Source_B),
        0
    )

    if n_a != record.N_A or n_b != record.N_B:
        raise ValueError("Sample count mismatch")

OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# RENDER ALL SEVEN FIGURES IN TEMPORARY DIRECTORY
# ============================================================

with tempfile.TemporaryDirectory(
    prefix="Figure3G_staging_",
    dir=str(ROOT)
) as temporary:

    STAGE = Path(temporary)

    staged = []

    for index, species in enumerate(SPECIES):

        rng = np.random.default_rng(20261010 + index)

        species_data = genomes.loc[
            genomes["Species"].eq(species)
        ]

        if species_data.empty:
            raise ValueError(
                f"No genomes available for {species}"
            )

        species_max = float(species_data["AMR"].max())
        y_max = max(7.0, species_max + 4.0)

        fig, axes = plt.subplots(
            1, 3,
            figsize=(17.3, 9.5),
            sharey=True,
            facecolor="white"
        )

        # Reduced spacing between panels.
        fig.subplots_adjust(
            left=.073,
            right=.978,
            top=.690,
            bottom=.190,
            wspace=.095
        )

        # Frame around entire figure.
        outer = FancyBboxPatch(
            (.008, .010),
            .984,
            .978,
            transform=fig.transFigure,
            boxstyle="round,pad=0.003,rounding_size=0.012",
            facecolor="none",
            edgecolor="#52678E",
            linewidth=1.15,
            zorder=30
        )

        fig.add_artist(outer)

        fig.suptitle(
            species,
            fontsize=24,
            fontstyle="italic",
            fontweight="bold",
            color="#14274B",
            y=.974
        )

        fig.text(
            .5,
            .926,
            "One Health variation in annotated AMR gene burden per genome",
            ha="center",
            fontsize=12,
            color="#35445B"
        )

        for panel_index, income in enumerate(INCOMES):

            ax = axes[panel_index]
            position = ax.get_position()

            # Income-specific rounded panel background.
            # Frame extends to encompass annotations.
            left = position.x0 - .006
            right = position.x1 + .006

            panel_frame = FancyBboxPatch(
                (left, .159),
                right-left,
                .714,
                transform=fig.transFigure,
                boxstyle=(
                    "round,pad=0.005,"
                    "rounding_size=0.014"
                ),
                facecolor=INCOME_TINTS[income],
                edgecolor=lighten(
                    INCOME_COLORS[income], .25
                ),
                linewidth=1.45,
                zorder=-4
            )

            fig.add_artist(panel_frame)

            # Saturated colored header band.
            header = FancyBboxPatch(
                (left, .823),
                right-left,
                .050,
                transform=fig.transFigure,
                boxstyle=(
                    "round,pad=0.004,"
                    "rounding_size=0.012"
                ),
                facecolor=INCOME_COLORS[income],
                edgecolor="none",
                linewidth=0,
                zorder=-2
            )

            fig.add_artist(header)

            fig.text(
                (left+right)/2,
                .848,
                INCOME_TITLES[income],
                fontsize=13.2,
                fontweight="bold",
                color="white",
                ha="center",
                va="center",
                zorder=20
            )

            ax.set_facecolor(INCOME_TINTS[income])

            panel_data = species_data.loc[
                species_data["Income_Group"].eq(income)
            ]

            for source in SOURCES:

                subset = panel_data.loc[
                    panel_data["Source"].eq(source)
                ]

                values = subset["AMR"].to_numpy(dtype=float)
                x = POSITIONS[source]

                if not len(values):
                    ax.text(
                        x, y_max*.5,
                        "No genomes",
                        ha="center",
                        va="center",
                        rotation=90,
                        color="#7E8896"
                    )
                    continue

                mean_value = draw_violin(
                    ax,
                    values,
                    x,
                    source,
                    rng
                )

                ax.text(
                    x,
                    .970,
                    f"n={len(values):,}\nMean={mean_value:.2f}",
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    va="top",
                    fontsize=9.2,
                    color="#162C52",
                    weight="medium"
                )

            # Brackets are above the data region.
            bracket_spec = [
                ("Human", "Animal", 1.055),
                ("Animal", "Environment", 1.133),
                ("Human", "Environment", 1.222)
            ]

            for a, b, height in bracket_spec:

                result = get_test(
                    tests, species, income, a, b
                )

                if result["Status"] != "TESTED":
                    continue

                draw_bracket(
                    ax,
                    a,
                    b,
                    height,
                    result["Significance"],
                    result["FDR"]
                )

            ax.set_xlim(.48, 3.52)
            ax.set_ylim(0, y_max)

            ax.set_xticks([1, 2, 3])
            ax.set_xticklabels(
                SOURCES,
                fontsize=10.4,
                color="#172B4D"
            )

            ax.tick_params(
                axis="x",
                length=0,
                pad=10
            )

            ax.tick_params(
                axis="y",
                labelsize=9.5,
                colors="#24354F"
            )

            ax.grid(
                axis="y",
                color="#ABB5C3",
                alpha=.27,
                linewidth=.65
            )

            ax.set_axisbelow(True)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["bottom"].set_color("#657184")
            ax.spines["left"].set_color("#657184")

            if panel_index == 0:
                ax.set_ylabel(
                    "Unique annotated AMR genes per genome",
                    fontsize=11.6,
                    labelpad=8,
                    color="#172B4D"
                )
            else:
                ax.tick_params(axis="y", labelleft=False)

        # Unified legend.
        handles = [
            Patch(
                facecolor=SOURCE_COLORS[source],
                edgecolor=darken(
                    SOURCE_COLORS[source], .25
                ),
                label=source
            )
            for source in SOURCES
        ]

        fig.legend(
            handles=handles,
            loc="lower center",
            bbox_to_anchor=(.5, .105),
            ncol=3,
            frameon=True,
            fancybox=True,
            edgecolor="#C1CAD6",
            facecolor="#FFFFFF",
            fontsize=11
        )

        fig.text(
            .5,
            .067,
            "Black diamonds indicate means; horizontal lines indicate medians.",
            ha="center",
            fontsize=9.4,
            color="#45556E"
        )

        fig.text(
            .5,
            .046,
            "Two-sided Mann–Whitney U; q = BH-adjusted P-value. "
            "Deep blue: ns. Magenta: significant.",
            ha="center",
            fontsize=9.1,
            color=DEEP_BLUE
        )

        fig.text(
            .5,
            .026,
            "Exploratory genome-level inference; "
            "not adjusted for country, BioProject or lineage.",
            ha="center",
            fontsize=8.5,
            color="#667185"
        )

        tag = species.replace(
            "Enterobacter ", "E_"
        ).replace(" ", "_")

        file = STAGE / (
            f"Fig03G_{tag}_FINAL_PUBLICATION.jpg"
        )

        fig.savefig(
            file,
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

        # Verify image can actually be opened.
        with Image.open(file) as im:
            im.verify()

        if file.stat().st_size < 10000:
            raise RuntimeError(
                "Invalid output image: " + str(file)
            )

        staged.append(file)

        print(
            f"[RENDERED {len(staged)}/7] "
            f"{species}",
            flush=True
        )

    if len(staged) != 7:
        raise RuntimeError(
            "Seven images were not generated"
        )

    # Copy and verify every final image.
    final_files = []

    for source_file in staged:

        destination = OUT / source_file.name
        shutil.copy2(source_file, destination)

        with Image.open(destination) as im:
            im.verify()

        final_files.append(destination.resolve())

    if len(final_files) != 7:
        raise RuntimeError(
            "Seven final files were not verified"
        )

    # Delete ONLY older Figure 3G JPEG/PNG files.
    # No Excel, scripts, other figure directories.
    final_set = set(final_files)
    removed = 0

    for old in ROOT.rglob("*"):

        if not old.is_file():
            continue

        if old.suffix.lower() not in {
            ".jpg", ".jpeg", ".png"
        }:
            continue

        resolved = old.resolve()

        if resolved in final_set:
            continue

        if STAGE.resolve() in resolved.parents:
            continue

        old.unlink()
        removed += 1

    print("\nOld Figure 3G image files deleted:", removed)

print("\n" + "=" * 90)
print("FIGURE 3G — FINAL PUBLICATION FIGURES COMPLETED")
print("=" * 90)

print("Species figures:", 7)
print("Resolution:", DPI, "DPI")
print("Original statistical workbook:", INPUT)
print("New figure directory:", OUT)
print("P-values and FDR not recalculated.")
print("=" * 90)
