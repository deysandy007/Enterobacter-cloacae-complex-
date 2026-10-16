
from pathlib import Path
from datetime import datetime
import gc

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib import font_manager

# ============================================================
# FIGURE 3G - GENOME-LEVEL ONE HEALTH AMR BURDEN
#
# Seven resolved Enterobacter cloacae complex species
# Three income categories:
# H / M / L+LM
#
# Three One Health sources:
# Human / Animal / Environment
#
# Outcome:
# Unique_AMR_Gene_Count per genome
#
# Figure:
# Seven individual species figures
# Three income panels per figure
# Boxplot + individual genome observations
#
# Descriptive analysis only:
# NO misleading genome-level significance tests
#
# Outputs:
# 7 JPEG figures at 600 DPI
# 1 master Excel workbook
#
# Data:
# Final 9157-genome ECC master
# Frozen Figure 1A country income assignments
#
# Previous figures and files are not modified.
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

OUT = ROOT / "Fig_03G_OneHealth_Genome_Burden"
OUT.mkdir(parents=True, exist_ok=True)

EXCEL = OUT / "Fig03G_OneHealth_AMR_Burden_Master.xlsx"

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

INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82"
}

SOURCE_COLORS = {
    "Human": "#E69F00",
    "Animal": "#56B4E9",
    "Environment": "#009E73"
}

SPECIES_COLORS = {
    "Enterobacter hormaechei": "#0072B2",
    "Enterobacter cloacae": "#D55E00",
    "Enterobacter roggenkampii": "#009E73",
    "Enterobacter asburiae": "#CC79A7",
    "Enterobacter kobei": "#7851A9",
    "Enterobacter ludwigii": "#A6761D",
    "Enterobacter bugandensis": "#56B4E9"
}

RNG = np.random.default_rng(20261010)

DPI = 600
MAX_DISPLAY_POINTS_PER_GROUP = 350

# ============================================================
# FONTS AND STYLE
# ============================================================

installed = {
    f.name for f in font_manager.fontManager.ttflist
}

FONT = next(
    (
        name for name in [
            "Arial",
            "Helvetica",
            "Liberation Sans",
            "DejaVu Sans"
        ] if name in installed
    ),
    "DejaVu Sans"
)

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 11,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "axes.linewidth": 0.9,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": DPI
})

# ============================================================
# FUNCTIONS
# ============================================================

def short_species(name):
    return name.replace("Enterobacter ", "E. ")


def filename_species(name):
    return name.replace("Enterobacter ", "E_").replace(
        " ", "_"
    )


def draw_source_box(ax, position, values, color):

    bp = ax.boxplot(
        [values],
        positions=[position],
        widths=0.54,
        patch_artist=True,
        showfliers=False,
        whis=1.5,
        medianprops={
            "color": "#202020",
            "linewidth": 1.6
        },
        boxprops={
            "linewidth": 1.0,
            "edgecolor": "#444444"
        },
        whiskerprops={
            "color": "#444444",
            "linewidth": 1.0
        },
        capprops={
            "color": "#444444",
            "linewidth": 1.0
        }
    )

    bp["boxes"][0].set_facecolor(color)
    bp["boxes"][0].set_alpha(0.48)

    # All genomes are used in the boxplot.
    # Points are downsampled only for visual clarity.
    if len(values) > MAX_DISPLAY_POINTS_PER_GROUP:

        selected = RNG.choice(
            values,
            size=MAX_DISPLAY_POINTS_PER_GROUP,
            replace=False
        )

    else:
        selected = values

    jitter = RNG.uniform(
        -0.19,
        0.19,
        size=len(selected)
    )

    ax.scatter(
        np.full(len(selected), position) + jitter,
        selected,
        s=7,
        facecolor=color,
        edgecolor="none",
        alpha=0.24,
        zorder=3,
        rasterized=True
    )


# ============================================================
# LOAD FINAL MASTER
# ============================================================

print("=" * 110)
print("FIGURE 3G - GENOME-LEVEL ONE HEALTH AMR BURDEN")
print("=" * 110)

if not MASTER.is_file():
    raise FileNotFoundError(MASTER)

if not INCOME_FILE.is_file():
    raise FileNotFoundError(INCOME_FILE)

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
    "Country_Harmonized",
    "BioProject",
    "Verified_Collection_Year",
    "Unique_AMR_Gene_Count",
    "AMRFinder_Status"
]

missing = [
    col for col in required
    if col not in df.columns
]

if missing:
    raise ValueError(
        "Missing master columns: " + str(missing)
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
        "Missing or duplicate assembly accessions"
    )

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

df["Country"] = (
    df["Country_Harmonized"]
    .fillna("")
    .astype(str)
    .str.strip()
)

df["AMR"] = pd.to_numeric(
    df["Unique_AMR_Gene_Count"],
    errors="coerce"
)

df["Year"] = pd.to_numeric(
    df["Verified_Collection_Year"],
    errors="coerce"
)

status = (
    df["AMRFinder_Status"]
    .fillna("")
    .astype(str)
    .str.strip()
    .str.upper()
)

if not status.eq("COMPLETE").all():
    raise ValueError(
        "Some AMRFinder annotations are incomplete"
    )

if (
    df["AMR"].isna().any()
    or (df["AMR"] < 0).any()
    or (df["AMR"] % 1 != 0).any()
):
    raise ValueError(
        "Invalid unique AMR gene counts"
    )

if df["Country"].eq("").any():
    raise ValueError(
        "Missing harmonized country names"
    )

# ============================================================
# FROZEN INCOME CLASSIFICATION
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
            "Missing income column: " + col
        )

income["Country"] = (
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

if income["Country"].duplicated().any():
    raise ValueError(
        "Duplicate income country in reference"
    )

audit = pd.concat(
    [
        df["Country"].value_counts().rename("Master"),
        income.set_index("Country")[
            "Genome_Count"
        ].rename("Income_Reference")
    ],
    axis=1
)

if (
    audit.isna().any().any()
    or not (
        audit["Master"].astype(int)
        == audit["Income_Reference"].astype(int)
    ).all()
):
    raise ValueError(
        "Income reference does not match master"
    )

df = df.merge(
    income[["Country", "Income_Group"]],
    on="Country",
    how="left",
    validate="many_to_one"
)

if df["Income_Group"].isna().any():
    raise ValueError(
        "Missing frozen income assignments"
    )

# ============================================================
# FILTER THE THREE REQUESTED ONE HEALTH SOURCES
# ============================================================

data = df[
    df["Species"].isin(SPECIES)
    & df["Source"].isin(SOURCES)
    & df["Income_Group"].isin(INCOMES)
].copy()

print("\nFinal master genomes:", len(df))
print("Analyzed genomes:", len(data))

print("\nSOURCE DISTRIBUTION")
print(
    data["Source"]
    .value_counts()
    .to_string()
)

print("\nINCOME DISTRIBUTION")
print(
    data["Income_Group"]
    .value_counts()
    .to_string()
)

# ============================================================
# DESCRIPTIVE STATISTICS
# ============================================================

full_index = pd.MultiIndex.from_product(
    [SPECIES, INCOMES, SOURCES],
    names=["Species", "Income_Group", "Source"]
)

summary = (
    data.groupby(
        ["Species", "Income_Group", "Source"]
    )
    .agg(
        Genome_N=("AMR", "size"),
        Country_N=("Country", "nunique"),
        BioProject_N=("BioProject", "nunique"),
        Mean_AMR=("AMR", "mean"),
        Median_AMR=("AMR", "median"),
        SD_AMR=("AMR", "std"),
        Q1_AMR=("AMR", lambda x: x.quantile(.25)),
        Q3_AMR=("AMR", lambda x: x.quantile(.75)),
        Minimum_AMR=("AMR", "min"),
        Maximum_AMR=("AMR", "max")
    )
    .reindex(full_index)
    .reset_index()
)

for col in [
    "Genome_N",
    "Country_N",
    "BioProject_N"
]:
    summary[col] = (
        summary[col]
        .fillna(0)
        .astype(int)
    )

if summary["Genome_N"].sum() != len(data):
    raise RuntimeError(
        "Descriptive summary does not reconcile"
    )

country_summary = (
    data.groupby(
        [
            "Species",
            "Income_Group",
            "Source",
            "Country"
        ]
    )
    .agg(
        Genome_N=("AMR", "size"),
        Mean_AMR=("AMR", "mean"),
        Median_AMR=("AMR", "median"),
        BioProject_N=("BioProject", "nunique")
    )
    .reset_index()
)

country_summary["Country_Eligible_N5"] = (
    country_summary["Genome_N"] >= 5
)

# ============================================================
# GENOME-LEVEL DATA
# ============================================================

genome_table = data[[
    "Assembly_Accession",
    "Species",
    "Income_Group",
    "Source",
    "Country",
    "BioProject",
    "Year",
    "AMR"
]].copy()

genome_table = genome_table.rename(columns={
    "AMR": "Unique_AMR_Gene_Count"
})

# ============================================================
# WRITE EXCEL FIRST
# ============================================================

print("\nWriting Excel workbook...")

methods = pd.DataFrame([
    (
        "Primary outcome",
        "Unique annotated AMR genes per genome"
    ),
    (
        "Master",
        "Final country-reported ECC master, "
        "9157 genomes (2000-2026)"
    ),
    (
        "Included genomes",
        len(data)
    ),
    (
        "Species",
        "Seven resolved ECC species"
    ),
    (
        "Sources",
        "Human; Animal; Environment"
    ),
    (
        "Income categories",
        "H; M; L+LM"
    ),
    (
        "Plot",
        "Species-specific boxplots in three "
        "income panels"
    ),
    (
        "Boxplot",
        "Median, Q1, Q3 and 1.5 x IQR whiskers"
    ),
    (
        "All observations",
        "Included in summary and boxplot calculations"
    ),
    (
        "Scatter display",
        f"Up to {MAX_DISPLAY_POINTS_PER_GROUP} "
        "genome points per group for readability"
    ),
    (
        "Statistical tests",
        "None; descriptive analysis only"
    ),
    (
        "Country clustering",
        "Not ignored for inferential testing; "
        "no independent-genome P-values reported"
    ),
    (
        "AMR caveat",
        "Annotated gene repertoire, not acquired-only "
        "or phenotypic resistance"
    ),
    (
        "Income colors",
        str(INCOME_COLORS)
    ),
    (
        "Source colors",
        str(SOURCE_COLORS)
    ),
    (
        "Created",
        datetime.now().isoformat()
    )
], columns=["Parameter", "Description"])

with pd.ExcelWriter(
    EXCEL,
    engine="openpyxl"
) as writer:

    summary.to_excel(
        writer,
        sheet_name="ALL_63_GROUPS",
        index=False
    )

    genome_table.to_excel(
        writer,
        sheet_name="GENOME_LEVEL_DATA",
        index=False
    )

    country_summary.to_excel(
        writer,
        sheet_name="COUNTRY_LEVEL_SUMMARY",
        index=False
    )

    methods.to_excel(
        writer,
        sheet_name="METHODS",
        index=False
    )

print("Excel saved:", EXCEL)

# ============================================================
# FIGURE GENERATION
# ============================================================

print("\nCreating seven species-specific figures...")

global_ymax = max(
    8,
    float(data["AMR"].max()) + 3
)

results = []

for species in SPECIES:

    species_data = data[
        data["Species"] == species
    ]

    tag = filename_species(species)

    output = OUT / (
        f"Fig03G_{tag}_OneHealth_AMR_Burden.jpg"
    )

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(16.8, 7.3),
        facecolor="white",
        sharey=True
    )

    fig.subplots_adjust(
        left=0.085,
        right=0.985,
        top=0.80,
        bottom=0.21,
        wspace=0.18
    )

    for panel_index, income_group in enumerate(INCOMES):

        ax = axes[panel_index]

        for source_index, source in enumerate(SOURCES):

            group = species_data[
                species_data["Income_Group"].eq(
                    income_group
                )
                & species_data["Source"].eq(source)
            ]

            values = group[
                "AMR"
            ].to_numpy(dtype=float)

            position = source_index + 1

            if len(values) > 0:

                draw_source_box(
                    ax,
                    position,
                    values,
                    SOURCE_COLORS[source]
                )

                median = float(np.median(values))

                n = len(values)

                countries = group[
                    "Country"
                ].nunique()

                # Display sample size above each group.
                ax.text(
                    position,
                    global_ymax - 1.15,
                    f"n={n:,}\nC={countries}",
                    ha="center",
                    va="top",
                    fontsize=10,
                    color="#333333"
                )

                # Display exact median in the lower
                # annotation area.
                ax.text(
                    position,
                    0.025,
                    f"Med. {median:.1f}",
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    color="#444444"
                )

            else:

                ax.text(
                    position,
                    global_ymax * 0.5,
                    "No data",
                    ha="center",
                    va="center",
                    fontsize=10,
                    rotation=90,
                    color="#9A9A9A"
                )

                ax.text(
                    position,
                    global_ymax - 1.15,
                    "n=0\nC=0",
                    ha="center",
                    va="top",
                    fontsize=10,
                    color="#777777"
                )

        # Fixed income color strip above each panel.
        ax.add_patch(
            Rectangle(
                (0, 1.04),
                1,
                0.035,
                transform=ax.transAxes,
                facecolor=INCOME_COLORS[income_group],
                edgecolor="none",
                clip_on=False
            )
        )

        income_name = {
            "H": "High income",
            "M": "Upper-middle income",
            "L+LM": "Low + lower-middle income"
        }[income_group]

        ax.set_title(
            income_name,
            fontsize=14,
            fontweight="bold",
            pad=27
        )

        ax.set_xlim(0.4, 3.6)
        ax.set_ylim(0, global_ymax)

        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(
            SOURCES,
            fontsize=11
        )

        ax.grid(
            axis="y",
            color="#E5E5E5",
            linewidth=0.7
        )

        ax.set_axisbelow(True)

        ax.tick_params(
            axis="x",
            length=0,
            pad=10
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if panel_index == 0:

            ax.set_ylabel(
                "Unique annotated AMR genes per genome",
                fontsize=12,
                fontweight="medium"
            )

    fig.suptitle(
        short_species(species),
        fontsize=22,
        fontweight="bold",
        fontstyle="italic",
        y=0.974
    )

    fig.text(
        0.5,
        0.914,
        "One Health comparison of AMR gene burden "
        "within country-income categories",
        ha="center",
        fontsize=12,
        color="#444444"
    )

    legend_handles = [
        Patch(
            facecolor=SOURCE_COLORS[source],
            edgecolor="#555555",
            label=source,
            alpha=0.65
        )
        for source in SOURCES
    ]

    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.116),
        ncol=3,
        frameon=False,
        fontsize=11
    )

    fig.text(
        0.5,
        0.060,
        "Boxes: median and interquartile range; "
        "whiskers: 1.5 × IQR. "
        "n = genomes; C = represented countries.",
        ha="center",
        fontsize=10,
        color="#555555"
    )

    fig.text(
        0.5,
        0.037,
        "All genomes contribute to summary statistics; "
        "individual points may be subsampled for display. "
        "Descriptive visualization only.",
        ha="center",
        fontsize=9.6,
        color="#666666"
    )

    fig.savefig(
        output,
        dpi=DPI,
        format="jpeg",
        facecolor="white",
        pil_kwargs={
            "quality": 93,
            "subsampling": 0,
            "optimize": True
        }
    )

    plt.close(fig)
    gc.collect()

    if not output.is_file():
        raise RuntimeError(
            "Figure not created: " + str(output)
        )

    results.append({
        "Species": species,
        "Genomes": len(species_data),
        "File": str(output)
    })

    print(
        "[SAVED]",
        short_species(species),
        "| genomes:",
        len(species_data),
        "|",
        output.name,
        flush=True
    )

# ============================================================
# FINAL REPORT
# ============================================================

print("\n" + "=" * 110)
print("FIGURE 3G - ALL SPECIES COMPLETED")
print("=" * 110)

print("Analyzed genomes:", len(data))
print("Species:", len(SPECIES))
print("Income groups:", len(INCOMES))
print("One Health sources:", len(SOURCES))
print("Expected source-income-species groups:", len(summary))
print("JPEGs generated:", len(results))
print("Resolution: 600 DPI")
print("\nExcel:", EXCEL)

print("\nGENOME BURDEN SUMMARY")

print(
    summary[[
        "Species",
        "Income_Group",
        "Source",
        "Genome_N",
        "Country_N",
        "Mean_AMR",
        "Median_AMR"
    ]].round(3).to_string(index=False)
)

print("\nNo statistical significance tests were performed.")
print("No previous figures or workbooks were deleted.")
print("=" * 110)
