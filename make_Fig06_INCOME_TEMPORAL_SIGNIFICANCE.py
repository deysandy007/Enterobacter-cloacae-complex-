
#!/usr/bin/env python3

# ============================================================
# FIGURE 6B — ECC INCOME-STRATIFIED TEMPORAL AMR BURDEN
#
# Seven ECC species, 2000–2026
# Each panel: H, M and L+LM trajectories
# Mean ± SEM; complete five-year intervals
#
# Statistical analysis:
# - Kruskal–Wallis within each species × period
# - BH-FDR across all valid omnibus tests
# - Mann–Whitney U pairwise tests with Holm correction
# - Pairwise significance gated on omnibus q < 0.05
#
# Outputs:
# 600-DPI JPEG
# Complete statistical Excel workbook
#
# Input:
# Final 9,157-genome ECC master
# Frozen Figure 1A country-income database
#
# No dependency on previously generated Figure 6 files.
# ============================================================

from pathlib import Path
import re
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter

# ============================================================
# PATHS
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

OUT = (
    ROOT / "Fig_06_ECC_Temporal_AMR_Burden"
    / "INCOME_STRATIFIED_FINAL"
)

OUT.mkdir(parents=True, exist_ok=True)

JPG = (
    OUT /
    "Fig06B_ECC_Income_Temporal_AMR_Significance_600DPI.jpg"
)

XLS = (
    OUT /
    "Fig06B_ECC_Income_Temporal_AMR_Statistics.xlsx"
)

DPI = 600
MIN_N = 5

# To enforce a particular column, enter its name here.
# Leave None for automatic selection.
YEAR_COLUMN = None

YEAR_CANDIDATES = [
    "Isolation_Year",
    "Collection_Year",
    "Isolate_Collection_Year",
    "Year_of_Collection",
    "Collection_Date",
    "Isolation_Date",
    "Isolate_Collection_Date",
    "Sample_Collection_Date",
    "Year",
    "Release_Year",
    "Assembly_Release_Year",
    "Assembly_Release_Date",
    "Submission_Year",
    "Submission_Date",
    "Release_Date",
]

# ============================================================
# CATEGORIES
# ============================================================

SPECIES = [
    ("Enterobacter hormaechei", "#0072B2"),
    ("Enterobacter cloacae", "#D55E00"),
    ("Enterobacter roggenkampii", "#009E73"),
    ("Enterobacter asburiae", "#CC79A7"),
    ("Enterobacter kobei", "#7851A9"),
    ("Enterobacter ludwigii", "#A6761D"),
    ("Enterobacter bugandensis", "#56B4E9"),
]

INCOMES = ["H", "M", "L+LM"]

INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82",
}

INCOME_LABELS = {
    "H": "High income",
    "M": "Upper-middle income",
    "L+LM": "Lower-middle + low income",
}

PERIODS = [
    (2000, 2004),
    (2005, 2009),
    (2010, 2014),
    (2015, 2019),
    (2020, 2024),
    (2025, 2026),
]

PERIOD_LABELS = [
    f"{a}–{b}" for a, b in PERIODS
]

OFFSETS = {
    "H": -0.14,
    "M": 0.0,
    "L+LM": 0.14,
}

# ============================================================
# FUNCTIONS
# ============================================================

def get_year(value):
    if pd.isna(value):
        return np.nan

    if isinstance(value, pd.Timestamp):
        return value.year

    s = str(value).strip()

    if re.fullmatch(r"\d{4}(?:\.0+)?", s):
        return int(float(s))

    match = re.search(
        r"(?<!\d)(19\d{2}|20\d{2})(?!\d)",
        s
    )

    return int(match.group(1)) if match else np.nan


def get_period(year):
    if pd.isna(year):
        return None

    for i, (start, end) in enumerate(PERIODS):
        if start <= int(year) <= end:
            return i

    return None


def genes_from_cell(value):
    if pd.isna(value):
        return set()

    s = str(value).strip()

    if s.lower() in {
        "", "nan", "none", "na",
        "n/a", "no hits", "not detected", "-"
    }:
        return set()

    return {
        gene.strip()
        for gene in s.split(";")
        if gene.strip()
    }


def normalize_income(value):
    s = str(value).strip().upper()
    s = re.sub(r"[\s_\-]+", "", s)

    mapping = {
        "H": "H",
        "HIGH": "H",
        "HIGHINCOME": "H",
        "M": "M",
        "UM": "M",
        "UPPERMIDDLE": "M",
        "UPPERMIDDLEINCOME": "M",
        "L+LM": "L+LM",
        "LM+L": "L+LM",
        "LOWERMIDDLE+LOW": "L+LM",
        "LOW+LOWERMIDDLE": "L+LM",
        "LOWERMIDDLEANDLOW": "L+LM",
    }

    return mapping.get(s, "Unclassified")


def stars(p):
    if not np.isfinite(p):
        return "NA"
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def mean_sem(values):
    a = np.asarray(values, dtype=float)

    if len(a) == 0:
        return np.nan, np.nan

    mean = float(np.mean(a))

    if len(a) < 2:
        return mean, np.nan

    return (
        mean,
        float(
            np.std(a, ddof=1) / np.sqrt(len(a))
        )
    )


def style_workbook(writer):
    for ws in writer.book.worksheets:
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions

        for cell in ws[1]:
            cell.fill = PatternFill(
                "solid", fgColor="173352"
            )
            cell.font = Font(
                bold=True, color="FFFFFF"
            )

        for i, cell in enumerate(ws[1], 1):
            ws.column_dimensions[
                get_column_letter(i)
            ].width = min(
                45,
                max(15, len(str(cell.value)) + 3)
            )


# ============================================================
# LOAD
# ============================================================

print("=" * 100)
print("FIGURE 6B — INCOME-STRATIFIED ECC AMR BURDEN")
print("=" * 100)

for file in [MASTER, INCOME_FILE]:
    if not file.is_file():
        raise FileNotFoundError(file)

df = pd.read_excel(
    MASTER,
    sheet_name="FINAL_MASTER_2000_2026",
    engine="openpyxl",
    dtype={"Assembly_Accession": str}
)

if len(df) != 9157:
    raise ValueError(
        f"Expected 9157 genomes, found {len(df)}"
    )

required = [
    "Assembly_Accession",
    "Normalized_Species",
    "Country_Harmonized",
    "AMR_Genes",
]

missing = [c for c in required if c not in df.columns]

if missing:
    raise ValueError(
        f"Missing master columns: {missing}"
    )

if (
    df["Assembly_Accession"].isna().any()
    or df["Assembly_Accession"].duplicated().any()
):
    raise ValueError(
        "Missing or duplicated accession"
    )

income = pd.read_excel(
    INCOME_FILE,
    sheet_name="COUNTRY_INCOME_COUNTS",
    engine="openpyxl"
)

for col in ["Country_Harmonized", "Income_Group"]:
    if col not in income.columns:
        raise ValueError(
            f"Missing income column: {col}"
        )

income["Country_Key"] = (
    income["Country_Harmonized"]
    .fillna("")
    .astype(str)
    .str.strip()
)

income["Income_Group"] = (
    income["Income_Group"]
    .map(normalize_income)
)

if income["Country_Key"].duplicated().any():
    raise ValueError(
        "Duplicate countries in income database"
    )

income_map = income.set_index(
    "Country_Key"
)["Income_Group"]

df["Country_Key"] = (
    df["Country_Harmonized"]
    .fillna("")
    .astype(str)
    .str.strip()
)

df["Income_Used"] = (
    df["Country_Key"]
    .map(income_map)
    .fillna("Unclassified")
)

df["Species_Used"] = (
    df["Normalized_Species"]
    .astype(str)
    .str.strip()
)

# ============================================================
# YEAR FIELD
# ============================================================

if YEAR_COLUMN is not None:
    if YEAR_COLUMN not in df.columns:
        raise ValueError(
            f"YEAR_COLUMN missing: {YEAR_COLUMN}"
        )
    year_col = YEAR_COLUMN

else:
    year_col = next(
        (
            c for c in YEAR_CANDIDATES
            if c in df.columns
            and df[c].map(get_year).between(
                2000, 2026
            ).any()
        ),
        None
    )

if year_col is None:
    raise ValueError(
        "No usable year column found.\n"
        + "\n".join(df.columns)
    )

year_l = year_col.lower()

if "collection" in year_l or "isolation" in year_l:
    year_type = "Collection/isolation"
elif "release" in year_l or "submi" in year_l:
    year_type = "Submission/release"
else:
    year_type = "Unspecified"

df["Year_Used"] = df[year_col].map(get_year)
df["Period_Index"] = df["Year_Used"].map(
    get_period
)

df["AMR_Genes_Per_Genome"] = (
    df["AMR_Genes"].map(
        lambda x: len(genes_from_cell(x))
    )
)

if "AMRFinder_Status" in df.columns:
    df["Annotation_Failed"] = (
        df["AMRFinder_Status"]
        .fillna("")
        .astype(str)
        .str.lower()
        .str.contains(
            r"fail|error|not.run|not processed",
            regex=True
        )
    )
else:
    df["Annotation_Failed"] = False

valid_species = {
    species for species, _ in SPECIES
}

df["Figure_Eligible"] = (
    df["Species_Used"].isin(valid_species)
    & df["Income_Used"].isin(INCOMES)
    & df["Period_Index"].notna()
    & ~df["Annotation_Failed"]
)

work = df[
    df["Figure_Eligible"]
].copy()

if work.empty:
    raise ValueError(
        "No eligible genomes for Figure 6B"
    )

work["Period_Index"] = (
    work["Period_Index"].astype(int)
)

print("Master genomes:", len(df))
print("Income-eligible temporal genomes:", len(work))
print("Year column:", year_col)
print("Year interpretation:", year_type)

# ============================================================
# SUMMARY AND STATISTICAL ANALYSIS
# ============================================================

summary_rows = []
omnibus_rows = []
raw_pairwise = []

for species, _ in SPECIES:

    species_part = work[
        work["Species_Used"].eq(species)
    ]

    for period_i, period_name in enumerate(PERIOD_LABELS):

        period_part = species_part[
            species_part["Period_Index"].eq(period_i)
        ]

        arrays = {}

        for inc in INCOMES:

            a = period_part.loc[
                period_part["Income_Used"].eq(inc),
                "AMR_Genes_Per_Genome"
            ].to_numpy(float)

            arrays[inc] = a

            avg, sem = mean_sem(a)

            summary_rows.append({
                "Species": species,
                "Period": period_name,
                "Period_Index": period_i,
                "Income_Group": inc,
                "N": len(a),
                "Mean": avg,
                "SEM": sem,
                "Median": (
                    float(np.median(a))
                    if len(a) else np.nan
                ),
                "Q1": (
                    float(np.percentile(a, 25))
                    if len(a) else np.nan
                ),
                "Q3": (
                    float(np.percentile(a, 75))
                    if len(a) else np.nan
                )
            })

        sufficient = {
            key: value
            for key, value in arrays.items()
            if len(value) >= MIN_N
        }

        if len(sufficient) >= 2:
            try:
                stat, p = kruskal(
                    *sufficient.values()
                )
            except ValueError:
                stat, p = 0.0, 1.0
        else:
            stat, p = np.nan, np.nan

        omnibus_rows.append({
            "Species": species,
            "Period": period_name,
            "Period_Index": period_i,
            "Income_Groups_Tested": (
                "; ".join(sufficient)
            ),
            "Kruskal_Wallis_H": stat,
            "P_raw": p,
        })

        pairs = [
            ("H", "M"),
            ("H", "L+LM"),
            ("M", "L+LM"),
        ]

        for inc1, inc2 in pairs:

            a = arrays[inc1]
            b = arrays[inc2]

            if len(a) < MIN_N or len(b) < MIN_N:
                continue

            p_pair = mannwhitneyu(
                a, b,
                alternative="two-sided"
            ).pvalue

            raw_pairwise.append({
                "Species": species,
                "Period": period_name,
                "Period_Index": period_i,
                "Income_1": inc1,
                "Income_2": inc2,
                "N_1": len(a),
                "N_2": len(b),
                "Mean_Difference_2_minus_1": (
                    float(np.mean(b) - np.mean(a))
                ),
                "P_raw": p_pair
            })

summary = pd.DataFrame(summary_rows)
omnibus = pd.DataFrame(omnibus_rows)
pairwise = pd.DataFrame(raw_pairwise)

# BH across all testable species x period comparisons.
omnibus["Q_BH"] = np.nan

valid = omnibus["P_raw"].notna()

if valid.any():
    omnibus.loc[valid, "Q_BH"] = multipletests(
        omnibus.loc[valid, "P_raw"],
        method="fdr_bh"
    )[1]

# Holm correction among three income comparisons
# WITHIN each species-period comparison.
if not pairwise.empty:

    pairwise["P_Holm"] = np.nan

    for (_, _), group in pairwise.groupby(
        ["Species", "Period"],
        sort=False
    ):
        pairwise.loc[
            group.index, "P_Holm"
        ] = multipletests(
            group["P_raw"].to_numpy(float),
            method="holm"
        )[1]

    pairwise = pairwise.merge(
        omnibus[
            ["Species", "Period", "Q_BH"]
        ],
        on=["Species", "Period"],
        how="left",
        validate="many_to_one"
    )

    pairwise["Significant"] = (
        (pairwise["P_Holm"] < 0.05)
        & (pairwise["Q_BH"] < 0.05)
    )

    pairwise["Stars"] = [
        stars(p) if ok else "ns"
        for p, ok in zip(
            pairwise["P_Holm"],
            pairwise["Significant"]
        )
    ]

else:
    pairwise = pd.DataFrame(columns=[
        "Species", "Period", "Period_Index",
        "Income_1", "Income_2", "N_1",
        "N_2", "Mean_Difference_2_minus_1",
        "P_raw", "P_Holm", "Q_BH",
        "Significant", "Stars"
    ])

# ============================================================
# PLOT — 7 SPECIES WITH 3 INCOME TRAJECTORIES
# ============================================================

print("Creating seven-panel income temporal figure...")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "figure.facecolor": "white",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

fig, axes = plt.subplots(
    2, 4,
    figsize=(23.2, 14.4),
    facecolor="white"
)

axes = axes.ravel()

finite_means = summary["Mean"].dropna()
finite_sems = summary["SEM"].fillna(0)

global_max = float(
    (summary["Mean"] + finite_sems).max()
) if len(finite_means) else 5.0

global_max = max(4.0, global_max)

# Make room for significance annotations.
YMAX = global_max * 1.35

for panel_i, (species, species_color) in enumerate(SPECIES):

    ax = axes[panel_i]

    ss = summary[
        summary["Species"].eq(species)
    ]

    for inc in INCOMES:

        sub = ss[
            ss["Income_Group"].eq(inc)
        ].sort_values("Period_Index")

        x = (
            sub["Period_Index"]
            .to_numpy(float)
            + OFFSETS[inc]
        )

        y = sub["Mean"].to_numpy(float)
        sem = sub["SEM"].to_numpy(float)
        n = sub["N"].to_numpy(int)

        valid_points = (
            (n >= MIN_N)
            & np.isfinite(y)
        )

        sparse = (
            (n > 0)
            & ~valid_points
            & np.isfinite(y)
        )

        y_plot = np.where(
            valid_points, y, np.nan
        )

        sem_plot = np.where(
            valid_points,
            np.nan_to_num(sem, nan=0),
            np.nan
        )

        ax.errorbar(
            x,
            y_plot,
            yerr=sem_plot,
            color=INCOME_COLORS[inc],
            linewidth=2.1,
            marker="o",
            markersize=5.6,
            markeredgecolor="white",
            markeredgewidth=0.7,
            elinewidth=1.2,
            capsize=3,
            zorder=4
        )

        if sparse.any():
            ax.scatter(
                x[sparse],
                y[sparse],
                marker="x",
                color=INCOME_COLORS[inc],
                s=38,
                linewidth=1.5,
                zorder=6
            )

    # Print actual sample counts under each five-year period.
    # Each triplet is in order H/M/L+LM.
    counts_by_period = []

    for p in range(6):
        period_sub = ss[
            ss["Period_Index"].eq(p)
        ]

        count_map = dict(
            zip(
                period_sub["Income_Group"],
                period_sub["N"]
            )
        )

        counts_by_period.append(
            "/".join(
                str(int(count_map.get(inc, 0)))
                for inc in INCOMES
            )
        )

    ax.set_xlim(-0.48, 5.48)
    ax.set_ylim(0, YMAX)

    ax.set_xticks(range(6))
    ax.set_xticklabels(
        PERIOD_LABELS,
        rotation=40,
        ha="right",
        fontsize=8.1
    )

    ax.set_xlabel(
        "Five-year interval",
        fontsize=9.5,
        labelpad=7
    )

    ax.set_ylabel(
        "AMR genes per genome (mean ± SEM)",
        fontsize=9.0,
        labelpad=5
    )

    ax.yaxis.set_major_locator(
        MaxNLocator(nbins=5)
    )

    ax.grid(
        axis="y",
        linestyle=":",
        alpha=0.24
    )

    ax.set_title(
        species.replace("Enterobacter", "E."),
        fontsize=13.0,
        fontstyle="italic",
        fontweight="bold",
        color=species_color,
        loc="left",
        pad=15
    )

    # Significance follows the reference-paper approach.
    # Show only omnibus FDR-significant periods as stars.
    # Specific income comparisons remain in Excel.
    tests = omnibus[
        omnibus["Species"].eq(species)
    ]

    for r in tests.itertuples(index=False):

        if (
            not np.isfinite(r.Q_BH)
            or r.Q_BH >= 0.05
        ):
            continue

        p_idx = int(r.Period_Index)

        period_data = ss[
            ss["Period_Index"].eq(p_idx)
        ]

        tops = (
            period_data["Mean"]
            + period_data["SEM"].fillna(0)
        ).dropna()

        top = (
            float(tops.max())
            if len(tops) else 0.0
        )

        yy = min(
            YMAX * 0.88,
            top + YMAX * 0.11
        )

        # Bracket spans three income positions.
        ax.plot(
            [
                p_idx - 0.19,
                p_idx - 0.19,
                p_idx + 0.19,
                p_idx + 0.19
            ],
            [
                yy - YMAX * 0.015,
                yy,
                yy,
                yy - YMAX * 0.015
            ],
            color="#283B50",
            linewidth=1.0,
            zorder=7
        )

        ax.text(
            p_idx,
            yy + YMAX * 0.015,
            stars(r.Q_BH),
            ha="center",
            va="bottom",
            fontsize=11.0,
            fontweight="bold",
            color="#283B50",
            zorder=8
        )

    # Period-specific group counts are displayed in a
    # separate strip to avoid overlap with graph points.
    for p, counts in enumerate(counts_by_period):

        ax.text(
            p,
            -0.47,
            counts,
            ha="center",
            va="top",
            fontsize=6.9,
            color="#4F6174",
            transform=ax.get_xaxis_transform(),
            clip_on=False
        )

    ax.text(
        -0.02,
        -0.35,
        "n H/M/L+LM:",
        ha="right",
        va="top",
        fontsize=7.3,
        transform=ax.transAxes,
        color="#405468",
        clip_on=False
    )

# Eighth panel: clean legend and statistical key.
axes[7].axis("off")

handles = [
    Line2D(
        [0], [0],
        color=INCOME_COLORS[inc],
        marker="o",
        markersize=7,
        linewidth=2.5,
        label=INCOME_LABELS[inc]
    )
    for inc in INCOMES
]

axes[7].legend(
    handles=handles,
    loc="upper left",
    frameon=False,
    fontsize=11.0,
    title="Country-income categories",
    title_fontsize=12.0
)

axes[7].text(
    0.04,
    0.55,
    "Points: mean ± SEM\n"
    "Income-group n: shown below each period\n"
    "Bracket: income-group omnibus test\n"
    "Kruskal–Wallis with BH-FDR correction\n"
    "Pairwise Mann–Whitney with Holm correction\n"
    "×: fewer than five genomes\n"
    "2025–2026: partial two-year interval",
    transform=axes[7].transAxes,
    fontsize=10.1,
    color="#344B61",
    va="top",
    linespacing=1.65
)

axes[7].text(
    0.04,
    0.10,
    "* q < 0.05     ** q < 0.01\n"
    "*** q < 0.001     **** q < 0.0001",
    transform=axes[7].transAxes,
    fontsize=10.2,
    fontweight="medium",
    color="#344B61"
)

fig.suptitle(
    "Income-stratified temporal dynamics of AMR gene burden "
    "in the Enterobacter cloacae complex",
    fontsize=19.7,
    fontweight="bold",
    color="#173352",
    y=0.992
)

fig.text(
    0.5, 0.967,
    "Seven species analyzed separately across six time intervals",
    ha="center",
    fontsize=11,
    color="#52697E"
)

fig.text(
    0.5, 0.017,
    f"Time field: {year_col} ({year_type}). "
    "Burden represents unique annotated AMR genes per genome, "
    "including intrinsic-associated annotations.",
    ha="center",
    fontsize=8.8,
    color="#5B6F82"
)

fig.subplots_adjust(
    left=0.073,
    right=0.986,
    top=0.928,
    bottom=0.165,
    wspace=0.38,
    hspace=0.94
)

fig.savefig(
    JPG,
    dpi=DPI,
    format="jpeg",
    facecolor="white",
    pil_kwargs={
        "quality": 93,
        "subsampling": 1
    }
)

plt.close(fig)

print("[SAVED]", JPG)

# ============================================================
# EXCEL SOURCE AND STATISTICS
# ============================================================

export_cols = [
    "Assembly_Accession",
    "Species_Used",
    "Country_Harmonized",
    "Income_Used",
    "Year_Used",
    "Period_Index",
    "AMR_Genes_Per_Genome",
    "AMR_Genes",
    "Annotation_Failed",
    "Figure_Eligible",
]

genome_export = df[export_cols].copy()

genome_export["Period"] = genome_export[
    "Period_Index"
].map(
    lambda x: (
        PERIOD_LABELS[int(x)]
        if pd.notna(x) else ""
    )
)

qc = pd.DataFrame([
    ["Master genomes", len(df)],
    ["Temporal income-eligible genomes", len(work)],
    ["Year column", year_col],
    ["Year meaning", year_type],
    ["Minimum group size", MIN_N],
    ["Omnibus test", "Kruskal-Wallis by species and period"],
    ["Omnibus correction", "BH-FDR across all valid species-period tests"],
    ["Pairwise test", "Mann-Whitney U"],
    ["Pairwise correction", "Holm within each species-period"],
    ["AMR burden", "Unique AMR annotations per genome"],
    ["Interpretation", "Sampled genome-level AMR burden"],
], columns=["Parameter", "Value"])

with pd.ExcelWriter(
    XLS,
    engine="openpyxl"
) as writer:

    genome_export.to_excel(
        writer,
        sheet_name="ALL_9157_GENOMES",
        index=False
    )

    summary.to_excel(
        writer,
        sheet_name="SPECIES_INCOME_PERIOD",
        index=False
    )

    omnibus.to_excel(
        writer,
        sheet_name="KRUSKAL_WALLIS_FDR",
        index=False
    )

    pairwise.to_excel(
        writer,
        sheet_name="PAIRWISE_HOLM",
        index=False
    )

    qc.to_excel(
        writer,
        sheet_name="QC_METHODS",
        index=False
    )

    style_workbook(writer)

print()
print("=" * 100)
print("FIGURE 6B FINISHED")
print("=" * 100)
print("Master genomes:", len(df))
print("Eligible income-stratified genomes:", len(work))
print("Year column:", year_col)
print("Year type:", year_type)
print("Statistical tests:", int(valid.sum()))
print("Significant omnibus tests:",
      int((omnibus["Q_BH"] < 0.05).sum()))
print("JPEG:", JPG)
print("Excel:", XLS)
print("=" * 100)
