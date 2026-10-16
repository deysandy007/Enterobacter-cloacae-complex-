
#!/usr/bin/env python3

# ============================================================
# FIGURE 7 — ECC TEMPORAL PLASMID-REPLICON ANALYSIS
#
# INPUT: 9157-genome ECC master + replicon annotations
#
# OUTPUTS:
#  A. 7 species figures: 3 income-specific stacked plots
#  B. 7 species replicon-prevalence heatmaps
#  C. Complete genome-level and statistical Excel
#
# Statistics:
#  Kruskal-Wallis across temporal periods
#  BH-FDR across species-income omnibus tests
#  Holm-adjusted planned pairwise comparisons
#
# IMPORTANT:
# Distinct replicon types per genome != plasmid count.
# ============================================================

from pathlib import Path
import re
import gc
from collections import Counter

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import Normalize

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

INCOME_DB = (
    ROOT / "Fig_01A_Global_Geographical_Distribution"
    / "Fig_01A_Country_Income_Database.xlsx"
)

OUT = ROOT / "Fig_07_ECC_Temporal_Replicon_Burden"
OUT.mkdir(parents=True, exist_ok=True)

TEMPORAL_DIR = OUT / "01_TEMPORAL_BURDEN"
HEATMAP_DIR = OUT / "02_REPLICON_PREVALENCE"

TEMPORAL_DIR.mkdir(exist_ok=True)
HEATMAP_DIR.mkdir(exist_ok=True)

EXCEL = OUT / "Fig07_REPLICON_COMPLETE_STATISTICS.xlsx"

# ------------------------------------------------------------
# OPTIONAL EXTERNAL REPLICON TABLE
#
# If your master already has a replicon column, leave None.
#
# If results are in a separate Excel/TSV/CSV file, supply
# its actual full Windows path here.
#
# The external table may have:
# 1. One accession per row + semicolon-separated replicons
# 2. Multiple rows per accession, one replicon per row
# ------------------------------------------------------------

REPLICON_FILE = None

# If the external workbook has multiple sheets, put the
# exact sheet name here; otherwise the first sheet is used.
REPLICON_SHEET = 0

# Optional explicit column names.
REPLICON_COLUMN = None
ACCESSION_COLUMN = None

YEAR_COLUMN = None

DPI = 600
MIN_N = 5
MAX_BRACKETS = 3
TOP_REPLICONS = 18

# ============================================================
# CONSTANTS
# ============================================================

SPECIES = [
    "Enterobacter hormaechei",
    "Enterobacter cloacae",
    "Enterobacter roggenkampii",
    "Enterobacter asburiae",
    "Enterobacter kobei",
    "Enterobacter ludwigii",
    "Enterobacter bugandensis"
]

SPECIES_COLORS = {
    "Enterobacter hormaechei": "#0072B2",
    "Enterobacter cloacae": "#D55E00",
    "Enterobacter roggenkampii": "#009E73",
    "Enterobacter asburiae": "#CC79A7",
    "Enterobacter kobei": "#7851A9",
    "Enterobacter ludwigii": "#A6761D",
    "Enterobacter bugandensis": "#56B4E9"
}

INCOMES = ["H", "M", "L+LM"]

INCOME_NAMES = {
    "H": "High-income countries",
    "M": "Upper-middle-income countries",
    "L+LM": "Lower-middle- and low-income countries"
}

INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82"
}

PERIODS = [
    (2000, 2004),
    (2005, 2009),
    (2010, 2014),
    (2015, 2019),
    (2020, 2024),
    (2025, 2026)
]

PERIOD_LABELS = [
    f"{a}\u2013{b}" for a, b in PERIODS
]

PLANNED_PAIRS = [
    (0, 1), (1, 2), (2, 3),
    (3, 4), (4, 5), (0, 5)
]

YEAR_CANDIDATES = [
    "Collection_Year",
    "Isolation_Year",
    "Isolate_Collection_Year",
    "Year_of_Collection",
    "Collection_Date",
    "Isolation_Date",
    "Year",
    "Assembly_Release_Year",
    "Release_Year",
    "Submission_Year",
    "Release_Date"
]

REPLICON_CANDIDATES = [
    "Replicon_Types",
    "Plasmid_Replicons",
    "Plasmid_Replicon_Types",
    "Replicons",
    "Replicon_Type",
    "Replicon",
    "PlasmidFinder_Replicons",
    "PlasmidFinder_Results",
    "Replicon_Names"
]

ACCESSION_CANDIDATES = [
    "Assembly_Accession",
    "assembly_accession",
    "Accession",
    "Genome_Accession",
    "Assembly",
    "Sample",
    "Sample_ID"
]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# ============================================================
# FUNCTIONS
# ============================================================

def match_column(frame, candidates, manual=None):

    if manual is not None:
        if manual not in frame.columns:
            raise ValueError(
                f"Requested column not found: {manual}\n"
                f"Available: {list(frame.columns)}"
            )
        return manual

    lookup = {
        str(c).strip().lower(): c
        for c in frame.columns
    }

    for name in candidates:
        if name.lower() in lookup:
            return lookup[name.lower()]

    return None


def parse_year(value):

    if pd.isna(value):
        return np.nan

    if isinstance(value, pd.Timestamp):
        return int(value.year)

    s = str(value).strip()

    match = re.search(
        r"(?<!\d)(19\d{2}|20\d{2})(?!\d)", s
    )

    return int(match.group(1)) if match else np.nan


def find_year_column(frame):

    if YEAR_COLUMN is not None:
        if YEAR_COLUMN not in frame.columns:
            raise ValueError(
                "Specified YEAR_COLUMN not found."
            )
        return YEAR_COLUMN

    for candidate in YEAR_CANDIDATES:
        col = match_column(frame, [candidate])

        if col is not None:
            years = frame[col].map(parse_year)

            if years.between(2000, 2026).any():
                return col

    raise ValueError(
        "No usable year column found.\n"
        + "\n".join(map(str, frame.columns))
    )


def period_index(year):

    if pd.isna(year):
        return np.nan

    for i, (a, b) in enumerate(PERIODS):
        if a <= int(year) <= b:
            return i

    return np.nan


def normalize_income(value):

    s = re.sub(
        r"[\s_\-]+", "",
        str(value).strip().upper()
    )

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
        "LOWERMIDDLEANDLOW": "L+LM"
    }

    return mapping.get(s, "Unclassified")


def parse_replicons(value):
    """
    Splits named replicon annotations only.
    Counts unique labels, not sequence copies.
    """

    if pd.isna(value):
        return set()

    s = str(value).strip()

    if s.lower() in {
        "", "nan", "none", "na", "n/a",
        "-", "no hit", "no hits", "not detected",
        "negative", "0"
    }:
        return set()

    # Semicolon, pipe and comma are common
    # delimiters in summary columns.
    tokens = re.split(r"[;|,\n]+", s)

    return {
        t.strip()
        for t in tokens
        if t.strip()
    }


def mean_sem(values):

    a = np.asarray(values, dtype=float)

    if len(a) == 0:
        return np.nan, np.nan

    mean = float(np.mean(a))

    if len(a) == 1:
        return mean, np.nan

    return (
        mean,
        float(a.std(ddof=1) / np.sqrt(len(a)))
    )


def star(p):

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


def ptext(p):

    if not np.isfinite(p):
        return "NA"

    if p < 0.0001:
        return f"{p:.2e}"

    return f"{p:.4f}"


def save_jpg(fig, filename):

    fig.savefig(
        filename,
        dpi=DPI,
        format="jpeg",
        facecolor="white",
        pil_kwargs={
            "quality": 93,
            "subsampling": 1
        }
    )

    plt.close(fig)
    gc.collect()

    print("[SAVED]", filename)


def read_external_table(path):

    suffix = path.suffix.lower()

    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(
            path,
            sheet_name=REPLICON_SHEET,
            engine="openpyxl"
            if suffix == ".xlsx" else None
        )

    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(
            path, sep="\t", low_memory=False
        )

    if suffix == ".csv":
        return pd.read_csv(
            path, low_memory=False
        )

    raise ValueError(
        "Replicon file must be XLSX, XLS, TSV or CSV."
    )


def style_excel(writer):

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

        for j, cell in enumerate(ws[1], 1):
            ws.column_dimensions[
                get_column_letter(j)
            ].width = min(
                45,
                max(16, len(str(cell.value)) + 3)
            )


# ============================================================
# LOAD MASTER
# ============================================================

print("=" * 95)
print("FIGURE 7 — TEMPORAL REPLICON BURDEN")
print("=" * 95)

if not MASTER.is_file():
    raise FileNotFoundError(MASTER)

if not INCOME_DB.is_file():
    raise FileNotFoundError(INCOME_DB)

df = pd.read_excel(
    MASTER,
    sheet_name="FINAL_MASTER_2000_2026",
    engine="openpyxl",
    dtype={"Assembly_Accession": str}
)

if len(df) != 9157:
    raise ValueError(
        f"Expected 9157 genomes; found {len(df)}"
    )

required = [
    "Assembly_Accession",
    "Normalized_Species",
    "Country_Harmonized"
]

for col in required:
    if col not in df.columns:
        raise ValueError(
            f"Missing master column: {col}"
        )

if (
    df["Assembly_Accession"].isna().any()
    or df["Assembly_Accession"].duplicated().any()
):
    raise ValueError(
        "Missing or duplicate master accession."
    )

# ============================================================
# LOAD REPLICON INFORMATION
# ============================================================

rep_col = match_column(
    df,
    REPLICON_CANDIDATES,
    REPLICON_COLUMN
)

if REPLICON_FILE is not None:

    external_path = Path(REPLICON_FILE)

    if not external_path.is_file():
        raise FileNotFoundError(
            f"External replicon file not found: "
            f"{external_path}"
        )

    rep = read_external_table(external_path)

    acc_col = match_column(
        rep,
        ACCESSION_CANDIDATES,
        ACCESSION_COLUMN
    )

    external_rep_col = match_column(
        rep,
        REPLICON_CANDIDATES,
        REPLICON_COLUMN
    )

    if acc_col is None or external_rep_col is None:
        raise ValueError(
            "Cannot identify accession and replicon "
            "columns in external file.\n"
            f"Columns: {list(rep.columns)}"
        )

    rep = rep[[acc_col, external_rep_col]].copy()

    rep["Join_Accession"] = (
        rep[acc_col].fillna("").astype(str).str.strip()
    )

    # Multiple hits/rows for the same genome
    # are combined into one unique set.
    rep["Rep_Set"] = (
        rep[external_rep_col].map(parse_replicons)
    )

    rep_by_accession = {}

    for accession, group in rep.groupby(
        "Join_Accession"
    ):
        if not accession:
            continue

        merged = set()

        for value in group["Rep_Set"]:
            merged.update(value)

        rep_by_accession[accession] = merged

    df["Rep_Set"] = (
        df["Assembly_Accession"]
        .astype(str)
        .str.strip()
        .map(rep_by_accession)
    )

    df["Replicon_Data_Available"] = (
        df["Rep_Set"].notna()
    )

    df["Rep_Set"] = df["Rep_Set"].map(
        lambda x: x if isinstance(x, set) else set()
    )

    REPLICON_SOURCE = str(external_path)

else:

    if rep_col is None:
        print("\nNO REPLICON COLUMN DETECTED.")
        print("\nMaster columns:")
        for c in df.columns:
            print(" ", c)

        print(
            "\nSet REPLICON_FILE near the top of "
            "the script to your actual replicon "
            "result table, then rerun."
        )

        raise ValueError(
            "No verified replicon source was supplied."
        )

    df["Rep_Set"] = (
        df[rep_col].map(parse_replicons)
    )

    # Where only a master column exists, missing/blank
    # replicon annotations are not automatically proof
    # that the annotation tool ran successfully.
    # Retain a separate data-availability flag.
    df["Replicon_Data_Available"] = (
        df[rep_col].notna()
    )

    REPLICON_SOURCE = (
        f"Master column: {rep_col}"
    )

# ============================================================
# FROZEN COUNTRY-INCOME MAPPING
# ============================================================

income = pd.read_excel(
    INCOME_DB,
    sheet_name="COUNTRY_INCOME_COUNTS",
    engine="openpyxl"
)

for col in [
    "Country_Harmonized",
    "Income_Group"
]:
    if col not in income.columns:
        raise ValueError(
            f"Missing income column: {col}"
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
        "Duplicate country in frozen income database"
    )

income_lookup = (
    income.set_index("Country_Key")["Income_Group"]
)

df["Country_Key"] = (
    df["Country_Harmonized"]
    .fillna("").astype(str).str.strip()
)

df["Income"] = (
    df["Country_Key"]
    .map(income_lookup)
    .fillna("Unclassified")
)

df["Species"] = (
    df["Normalized_Species"]
    .fillna("").astype(str).str.strip()
)

YEAR_USED = find_year_column(df)

df["Year"] = (
    df[YEAR_USED].map(parse_year)
)

df["Period_Index"] = (
    df["Year"].map(period_index)
)

df["Replicon_Types_Per_Genome"] = (
    df["Rep_Set"].map(len)
)

df["Replicon_List"] = (
    df["Rep_Set"].map(
        lambda x: "; ".join(sorted(x))
    )
)

df["Eligible"] = (
    df["Species"].isin(SPECIES)
    & df["Income"].isin(INCOMES)
    & df["Period_Index"].notna()
    & df["Replicon_Data_Available"]
)

work = df.loc[df["Eligible"]].copy()

if work.empty:
    raise ValueError(
        "No eligible genomes with replicon data."
    )

work["Period_Index"] = (
    work["Period_Index"].astype(int)
)

year_lower = YEAR_USED.lower()

if (
    "collection" in year_lower
    or "isolation" in year_lower
):
    YEAR_TYPE = "Collection/isolation year"
elif (
    "release" in year_lower
    or "submi" in year_lower
):
    YEAR_TYPE = "NCBI submission/release year"
else:
    YEAR_TYPE = "Year meaning unverified"

print("Master genomes:", len(df))
print("Eligible genomes:", len(work))
print("Replicon source:", REPLICON_SOURCE)
print("Year field:", YEAR_USED)
print("Year type:", YEAR_TYPE)

# ============================================================
# SUMMARY AND STATISTICS
# ============================================================

summary_rows = []
omnibus_rows = []
pairwise_rows = []

for species in SPECIES:

    for inc in INCOMES:

        sub = work[
            work["Species"].eq(species)
            & work["Income"].eq(inc)
        ]

        groups = {}

        for i, label in enumerate(PERIOD_LABELS):

            a = sub.loc[
                sub["Period_Index"].eq(i),
                "Replicon_Types_Per_Genome"
            ].to_numpy(float)

            groups[i] = a

            avg, sem = mean_sem(a)

            summary_rows.append({
                "Species": species,
                "Income": inc,
                "Period": label,
                "Period_Index": i,
                "N": len(a),
                "Mean_Replicons": avg,
                "SEM": sem,
                "Median": (
                    float(np.median(a))
                    if len(a) else np.nan
                ),
                "Q1": (
                    float(np.quantile(a, .25))
                    if len(a) else np.nan
                ),
                "Q3": (
                    float(np.quantile(a, .75))
                    if len(a) else np.nan
                )
            })

        adequate = [
            a for a in groups.values()
            if len(a) >= MIN_N
        ]

        if len(adequate) >= 2:
            try:
                h, p = kruskal(*adequate)
            except ValueError:
                h, p = 0.0, 1.0
        else:
            h, p = np.nan, np.nan

        omnibus_rows.append({
            "Species": species,
            "Income": inc,
            "Total_N": len(sub),
            "Kruskal_Wallis_H": h,
            "P_raw": p,
            "Testable_Periods": len(adequate)
        })

        family = []

        for i, j in PLANNED_PAIRS:

            a, b = groups[i], groups[j]

            if len(a) < MIN_N or len(b) < MIN_N:
                continue

            p_pair = float(
                mannwhitneyu(
                    a, b,
                    alternative="two-sided"
                ).pvalue
            )

            family.append({
                "Species": species,
                "Income": inc,
                "Index_1": i,
                "Index_2": j,
                "Period_1": PERIOD_LABELS[i],
                "Period_2": PERIOD_LABELS[j],
                "N_1": len(a),
                "N_2": len(b),
                "Mean_Difference": float(
                    np.mean(b) - np.mean(a)
                ),
                "P_raw": p_pair
            })

        if family:

            adjusted = multipletests(
                [r["P_raw"] for r in family],
                method="holm"
            )[1]

            for row, value in zip(family, adjusted):
                row["P_Holm"] = float(value)
                pairwise_rows.append(row)

summary = pd.DataFrame(summary_rows)
omnibus = pd.DataFrame(omnibus_rows)

omnibus["Q_BH"] = np.nan

valid = omnibus["P_raw"].notna()

if valid.any():
    omnibus.loc[valid, "Q_BH"] = (
        multipletests(
            omnibus.loc[valid, "P_raw"],
            method="fdr_bh"
        )[1]
    )

pair_columns = [
    "Species", "Income", "Index_1", "Index_2",
    "Period_1", "Period_2", "N_1", "N_2",
    "Mean_Difference", "P_raw", "P_Holm"
]

pairwise = pd.DataFrame(
    pairwise_rows,
    columns=pair_columns
)

pairwise = pairwise.merge(
    omnibus[["Species", "Income", "Q_BH"]],
    on=["Species", "Income"],
    how="left",
    validate="many_to_one"
)

pairwise["Significant"] = (
    (pairwise["P_Holm"] < .05)
    & (pairwise["Q_BH"] < .05)
)

pairwise["Stars"] = [
    star(p) if ok else "ns"
    for p, ok in zip(
        pairwise["P_Holm"],
        pairwise["Significant"]
    )
]

# ============================================================
# DRAW THREE STACKED REPLICON BURDEN GRAPHS
# ============================================================

def plot_burden(ax, stat_ax, species, inc):

    stat_ax.axis("off")

    color = INCOME_COLORS[inc]

    sub = summary[
        summary["Species"].eq(species)
        & summary["Income"].eq(inc)
    ].sort_values("Period_Index")

    x = sub["Period_Index"].to_numpy(int)
    y = sub["Mean_Replicons"].to_numpy(float)
    sem = sub["SEM"].to_numpy(float)
    ns = sub["N"].to_numpy(int)

    good = (ns >= MIN_N) & np.isfinite(y)
    sparse = (ns > 0) & ~good & np.isfinite(y)

    ax.errorbar(
        x,
        np.where(good, y, np.nan),
        yerr=np.where(
            good,
            np.nan_to_num(sem, nan=0),
            np.nan
        ),
        color=color,
        linewidth=2.4,
        marker="o",
        markersize=7,
        markeredgecolor="white",
        markeredgewidth=.9,
        capsize=4,
        elinewidth=1.6,
        zorder=5
    )

    if sparse.any():
        ax.scatter(
            x[sparse],
            y[sparse],
            marker="x",
            s=52,
            linewidth=1.6,
            color=color,
            zorder=6
        )

    tests = pairwise[
        pairwise["Species"].eq(species)
        & pairwise["Income"].eq(inc)
        & pairwise["Significant"]
    ].sort_values("P_Holm").head(MAX_BRACKETS)

    tops = y + np.nan_to_num(sem, nan=0)
    tops = tops[np.isfinite(tops)]

    observed = (
        max(1.0, float(tops.max()))
        if len(tops) else 1.0
    )

    step = max(.38, observed * .14)
    base = observed + step * .8

    ax.set_ylim(
        0,
        max(
            observed * 1.45,
            base + (len(tests) + 1)*step
        )
    )

    for level, row in enumerate(
        tests.itertuples(index=False)
    ):

        yy = base + level * step

        ax.plot(
            [
                row.Index_1, row.Index_1,
                row.Index_2, row.Index_2
            ],
            [
                yy-step*.16, yy,
                yy, yy-step*.16
            ],
            color="#263B50",
            linewidth=1.0,
            zorder=10
        )

        ax.text(
            (row.Index_1+row.Index_2)/2,
            yy+step*.05,
            f"{row.Stars}  p={ptext(row.P_Holm)}",
            ha="center",
            va="bottom",
            fontsize=8.7,
            fontweight="bold",
            color="#263B50"
        )

    ax.set_xticks(range(6))
    ax.set_xticklabels(
        [
            f"{period}\nn={n:,}"
            for period, n in zip(PERIOD_LABELS, ns)
        ],
        fontsize=8.9,
        linespacing=1.6
    )

    ax.set_xlim(-.4, 5.4)

    ax.set_xlabel(
        "Five-year interval",
        fontsize=10,
        labelpad=9
    )

    ax.set_ylabel(
        "Replicon types per genome\n(mean ± SEM)",
        fontsize=10
    )

    ax.set_title(
        INCOME_NAMES[inc],
        loc="left",
        fontsize=13,
        fontweight="bold",
        color=color,
        pad=12
    )

    ax.yaxis.set_major_locator(
        MaxNLocator(nbins=5)
    )

    ax.grid(
        axis="y",
        alpha=.25,
        linestyle=":"
    )

    ax.spines["left"].set_color(color)
    ax.spines["left"].set_linewidth(2)

    row = omnibus[
        omnibus["Species"].eq(species)
        & omnibus["Income"].eq(inc)
    ].iloc[0]

    stat_ax.text(
        .03, .90,
        "Kruskal–Wallis",
        fontsize=11,
        fontweight="bold",
        color="#173352",
        va="top",
        transform=stat_ax.transAxes
    )

    if np.isfinite(row["Kruskal_Wallis_H"]):
        info = (
            f"H = {row['Kruskal_Wallis_H']:.3f}\n"
            f"p = {ptext(row['P_raw'])}\n"
            f"FDR q = {ptext(row['Q_BH'])}"
        )
    else:
        info = "Not testable\nInsufficient groups"

    stat_ax.text(
        .03, .72,
        info,
        fontsize=10,
        color="#344A60",
        va="top",
        linespacing=1.8,
        transform=stat_ax.transAxes
    )

    stat_ax.text(
        .03, .20,
        f"Total n = {int(row['Total_N']):,}",
        fontsize=9,
        color="#50667A",
        transform=stat_ax.transAxes
    )


# ============================================================
# SEVEN SPECIES TEMPORAL FIGURES
# ============================================================

print("\nRendering temporal burden figures...")

for index, species in enumerate(SPECIES, 1):

    fig = plt.figure(figsize=(12.7, 17.7))

    grid = GridSpec(
        3, 2,
        figure=fig,
        width_ratios=[4.8, 1.35],
        left=.12,
        right=.97,
        top=.92,
        bottom=.12,
        hspace=.67,
        wspace=.18
    )

    for j, inc in enumerate(INCOMES):

        ax = fig.add_subplot(grid[j, 0])
        stats_ax = fig.add_subplot(grid[j, 1])

        plot_burden(ax, stats_ax, species, inc)

    fig.suptitle(
        species,
        fontsize=21,
        fontweight="bold",
        fontstyle="italic",
        color=SPECIES_COLORS[species],
        y=.982
    )

    fig.text(
        .5, .075,
        "Brackets: Holm-adjusted pairwise comparisons "
        "after Kruskal–Wallis FDR screening.",
        ha="center",
        fontsize=9,
        color="#3B5066"
    )

    fig.text(
        .5, .052,
        "* p<0.05   ** p<0.01   *** p<0.001   "
        "**** p<0.0001",
        ha="center",
        fontsize=9,
        color="#3B5066"
    )

    fig.text(
        .5, .030,
        f"Year field: {YEAR_USED} ({YEAR_TYPE}). "
        "Values are detected replicon types, not plasmid counts.",
        ha="center",
        fontsize=8.2,
        color="#65788A"
    )

    filename = (
        TEMPORAL_DIR
        / f"Fig07_{index:02d}_{species.split()[-1]}_Replicon_Temporal.jpg"
    )

    save_jpg(fig, filename)

# ============================================================
# PREVALENCE OF INDIVIDUAL REPLICON TYPES
# ============================================================

print("\nCalculating replicon-type prevalence...")

type_rows = []

for species in SPECIES:

    sp = work[work["Species"].eq(species)]

    for inc in INCOMES:

        si = sp[sp["Income"].eq(inc)]

        for period_i, period_label in enumerate(
            PERIOD_LABELS
        ):

            part = si[
                si["Period_Index"].eq(period_i)
            ]

            n = len(part)

            carriers = Counter()

            for rep_set in part["Rep_Set"]:
                carriers.update(rep_set)

            for rep_name, count in carriers.items():

                type_rows.append({
                    "Species": species,
                    "Income": inc,
                    "Period": period_label,
                    "Period_Index": period_i,
                    "Replicon": rep_name,
                    "Carriers": count,
                    "Total_Genomes": n,
                    "Prevalence_Percent": (
                        100.0 * count/n
                        if n else np.nan
                    )
                })

type_prevalence = pd.DataFrame(
    type_rows,
    columns=[
        "Species", "Income", "Period",
        "Period_Index", "Replicon",
        "Carriers", "Total_Genomes",
        "Prevalence_Percent"
    ]
)

# ============================================================
# SPECIES-LEVEL TEMPORAL REPLICON HEATMAPS
# Pooled over income groups for an uncluttered overview.
# ============================================================

print("\nRendering replicon prevalence heatmaps...")

for index, species in enumerate(SPECIES, 1):

    part = work[
        work["Species"].eq(species)
    ]

    total_carriers = Counter()

    for reps in part["Rep_Set"]:
        total_carriers.update(reps)

    selected = [
        rep for rep, count in
        total_carriers.most_common(TOP_REPLICONS)
    ]

    if not selected:
        print("[SKIP] No replicons:", species)
        continue

    matrix = np.full(
        (len(selected), len(PERIOD_LABELS)),
        np.nan
    )

    n_per_period = []

    for j, period in enumerate(PERIOD_LABELS):

        sub = part[
            part["Period_Index"].eq(j)
        ]

        n = len(sub)
        n_per_period.append(n)

        if n == 0:
            continue

        period_counts = Counter()

        for reps in sub["Rep_Set"]:
            period_counts.update(reps)

        for i, name in enumerate(selected):
            matrix[i, j] = (
                100.0 * period_counts[name] / n
            )

    fig, ax = plt.subplots(
        figsize=(11.8, max(7.5, len(selected)*.40+2.2))
    )

    image = ax.imshow(
        matrix,
        aspect="auto",
        cmap="YlGnBu",
        vmin=0,
        vmax=100,
        interpolation="nearest"
    )

    ax.set_yticks(np.arange(len(selected)))
    ax.set_yticklabels(
        selected,
        fontsize=9.2
    )

    ax.set_xticks(range(6))
    ax.set_xticklabels(
        [
            f"{p}\nn={n:,}"
            for p, n in zip(
                PERIOD_LABELS,
                n_per_period
            )
        ],
        fontsize=9.1
    )

    for i in range(len(selected)):
        for j in range(6):

            value = matrix[i, j]

            if not np.isfinite(value):
                label = "NA"
                text_color = "#344A60"
            else:
                label = f"{value:.1f}"
                text_color = (
                    "white" if value > 60
                    else "#1E344A"
                )

            ax.text(
                j, i,
                label,
                ha="center",
                va="center",
                fontsize=8.0,
                color=text_color
            )

    ax.set_xticks(
        np.arange(-.5, 6, 1),
        minor=True
    )

    ax.set_yticks(
        np.arange(-.5, len(selected), 1),
        minor=True
    )

    ax.grid(
        which="minor",
        color="white",
        linewidth=1
    )

    ax.tick_params(
        which="minor",
        left=False,
        bottom=False
    )

    ax.set_title(
        species,
        fontstyle="italic",
        fontweight="bold",
        fontsize=16,
        color=SPECIES_COLORS[species],
        pad=18
    )

    cbar = fig.colorbar(
        image,
        ax=ax,
        fraction=.037,
        pad=.025
    )

    cbar.set_label(
        "Genomes carrying replicon (%)",
        fontsize=10
    )

    fig.text(
        .5, .025,
        "Cells denote detected replicon carriage among "
        "genomes sampled during each period. "
        "NA = no genomes in period.",
        ha="center",
        fontsize=8.3,
        color="#50677C"
    )

    fig.subplots_adjust(
        left=.23,
        right=.91,
        top=.92,
        bottom=.14
    )

    filename = (
        HEATMAP_DIR
        / f"Fig07_{index:02d}_{species.split()[-1]}_Replicon_Prevalence.jpg"
    )

    save_jpg(fig, filename)

# ============================================================
# EXPORT EXCEL
# ============================================================

print("\nWriting comprehensive Excel...")

genome_export = df[[
    "Assembly_Accession",
    "Species",
    "Country_Harmonized",
    "Income",
    "Year",
    "Period_Index",
    "Replicon_List",
    "Replicon_Types_Per_Genome",
    "Replicon_Data_Available",
    "Eligible"
]].copy()

genome_export["Period"] = (
    genome_export["Period_Index"].map(
        lambda x: (
            PERIOD_LABELS[int(x)]
            if pd.notna(x) else ""
        )
    )
)

audit = pd.DataFrame([
    ["Master genomes", len(df)],
    ["Eligible genomes with replicon data", len(work)],
    ["Replicon data source", REPLICON_SOURCE],
    ["Year column", YEAR_USED],
    ["Year meaning", YEAR_TYPE],
    ["Metric", "Number of distinct replicon names per genome"],
    ["Not equivalent to", "Number of complete plasmids"],
    ["Omnibus test", "Kruskal-Wallis"],
    ["Omnibus correction", "BH-FDR across 21 species-income tests"],
    ["Pairwise test", "Mann-Whitney U"],
    ["Pairwise correction", "Holm within species-income"],
    ["Min n", MIN_N],
    ["Replicon prevalence denominator",
     "Genomes with available replicon results in group"],
    ["Replicon prevalence",
     "Genome carriers / eligible genomes x 100"],
    ["Caveat",
     "Replicon detection is assembly-dependent; absence requires a successfully run annotation"],
], columns=[
    "Parameter",
    "Description"
])

with pd.ExcelWriter(
    EXCEL,
    engine="openpyxl"
) as writer:

    genome_export.to_excel(
        writer,
        sheet_name="ALL_9157_GENOMES",
        index=False
    )

    summary.to_excel(
        writer,
        sheet_name="REPLICON_BURDEN_SUMMARY",
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

    type_prevalence.to_excel(
        writer,
        sheet_name="REPLICON_TYPE_PREVALENCE",
        index=False
    )

    audit.to_excel(
        writer,
        sheet_name="QC_METHODS",
        index=False
    )

    style_excel(writer)

print()
print("=" * 95)
print("FIGURE 7 COMPLETED")
print("=" * 95)
print("Master genomes:", len(df))
print("Eligible genomes:", len(work))
print("Replicon source:", REPLICON_SOURCE)
print("Temporal figures:", TEMPORAL_DIR)
print("Replicon heatmaps:", HEATMAP_DIR)
print("Excel:", EXCEL)
print("=" * 95)
