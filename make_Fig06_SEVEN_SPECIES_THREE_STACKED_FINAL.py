
#!/usr/bin/env python3

# ============================================================
# FIGURE 6 — FINAL SEVEN SPECIES × THREE INCOME GROUPS
#
# Each species has THREE LARGE VERTICALLY STACKED GRAPHS:
#   1. High-income countries
#   2. Upper-middle-income countries
#   3. Lower-middle- and low-income countries
#
# Each graph includes:
#   - Full five-year intervals
#   - Individual n values under every period
#   - Mean AMR gene burden ± SEM
#   - Kruskal-Wallis H, raw p and BH-adjusted q
#   - Holm-adjusted pairwise significance brackets
#   - The exact pairwise adjusted p beside each bracket
#
# Input:
#   Final 9,157 ECC genome master
#   Frozen country-income database
#
# Outputs:
#   Seven separate 600-DPI JPEGs
#   One detailed Excel workbook
#
# No dependencies on previous Figure 6 scripts.
# ============================================================

from pathlib import Path
import re
import math
import gc

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter

# ============================================================
# INPUT AND OUTPUT PATHS
# ============================================================

ROOT = Path(r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES")

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
    ROOT
    / "Fig_06_ECC_Temporal_AMR_Burden"
    / "SEVEN_SPECIES_THREE_STACKED_FINAL"
)

OUT.mkdir(parents=True, exist_ok=True)

EXCEL = OUT / "Fig06_Seven_Species_Income_Temporal_Statistics.xlsx"

DPI = 600

# Only groups with at least this many genomes
# enter the statistical testing.
MIN_N = 5

# Maximum number of significant brackets displayed
# in a graph. All comparisons remain in Excel.
MAX_BRACKETS = 4

# Specify the actual collection year field if necessary.
# Leave None for automatic identification.
YEAR_COLUMN = None

YEAR_CANDIDATES = [
    "Collection_Year",
    "Isolation_Year",
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
# SPECIES AND INCOME PALETTE
# ============================================================

SPECIES = [
    "Enterobacter hormaechei",
    "Enterobacter cloacae",
    "Enterobacter roggenkampii",
    "Enterobacter asburiae",
    "Enterobacter kobei",
    "Enterobacter ludwigii",
    "Enterobacter bugandensis",
]

SPECIES_COLORS = {
    "Enterobacter hormaechei": "#0072B2",
    "Enterobacter cloacae": "#D55E00",
    "Enterobacter roggenkampii": "#009E73",
    "Enterobacter asburiae": "#CC79A7",
    "Enterobacter kobei": "#7851A9",
    "Enterobacter ludwigii": "#A6761D",
    "Enterobacter bugandensis": "#56B4E9",
}

INCOMES = ["H", "M", "L+LM"]

INCOME_NAMES = {
    "H": "High-income countries",
    "M": "Upper-middle-income countries",
    "L+LM": "Lower-middle- and low-income countries",
}

INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82",
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
    f"{a}\N{EN DASH}{b}"
    for a, b in PERIODS
]

# Five adjacent comparisons plus first versus last.
PLANNED_COMPARISONS = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 4),
    (4, 5),
    (0, 5),
]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

# ============================================================
# HELPERS
# ============================================================

def year_from_cell(value):
    if pd.isna(value):
        return np.nan

    if isinstance(value, pd.Timestamp):
        return int(value.year)

    text = str(value).strip()

    if re.fullmatch(r"\d{4}(?:\.0+)?", text):
        return int(float(text))

    match = re.search(
        r"(?<!\d)(19\d{2}|20\d{2})(?!\d)",
        text
    )

    return int(match.group(1)) if match else np.nan


def find_year_column(data):
    if YEAR_COLUMN is not None:
        if YEAR_COLUMN not in data.columns:
            raise ValueError(
                "Specified year column not found: "
                + YEAR_COLUMN
            )
        return YEAR_COLUMN

    lookup = {c.lower(): c for c in data.columns}

    for candidate in YEAR_CANDIDATES:
        key = candidate.lower()

        if key not in lookup:
            continue

        actual = lookup[key]
        values = data[actual].map(year_from_cell)

        if values.between(2000, 2026).any():
            return actual

    raise ValueError(
        "No usable year column found.\n"
        "Available columns:\n"
        + "\n".join(data.columns)
    )


def period_index(year):
    if pd.isna(year):
        return np.nan

    for i, (start, end) in enumerate(PERIODS):
        if start <= int(year) <= end:
            return i

    return np.nan


def gene_set(value):
    if pd.isna(value):
        return set()

    text = str(value).strip()

    if text.lower() in {
        "", "nan", "none", "na", "n/a",
        "-", "no hits", "not detected"
    }:
        return set()

    return {
        g.strip()
        for g in text.split(";")
        if g.strip()
    }


def normalize_income(value):
    text = str(value).strip().upper()
    text = re.sub(r"[\s_\-]+", "", text)

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

    return mapping.get(text, "Unclassified")


def mean_sem(values):
    a = np.asarray(values, dtype=float)

    if len(a) == 0:
        return np.nan, np.nan

    mean = float(np.mean(a))

    if len(a) < 2:
        return mean, np.nan

    sem = float(
        np.std(a, ddof=1) / np.sqrt(len(a))
    )

    return mean, sem


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


def format_p(p):
    if not np.isfinite(p):
        return "NA"

    if p < 0.0001:
        return f"{p:.2e}"

    return f"{p:.4f}"


def excel_style(writer):
    for ws in writer.book.worksheets:

        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions

        for cell in ws[1]:
            cell.fill = PatternFill(
                "solid", fgColor="173352"
            )
            cell.font = Font(
                color="FFFFFF", bold=True
            )

        for i, cell in enumerate(ws[1], start=1):
            ws.column_dimensions[
                get_column_letter(i)
            ].width = min(
                48,
                max(16, len(str(cell.value)) + 3)
            )


def save_jpeg(fig, filename):
    fig.savefig(
        filename,
        dpi=DPI,
        format="jpeg",
        facecolor="white",
        pil_kwargs={
            "quality": 94,
            "subsampling": 1
        }
    )

    plt.close(fig)
    gc.collect()

    if not filename.is_file():
        raise RuntimeError(
            "Figure not saved: " + str(filename)
        )

    print(
        "[SAVED]",
        filename.name,
        f"({filename.stat().st_size / 1024**2:.1f} MB)"
    )


# ============================================================
# LOAD GENOME MASTER
# ============================================================

print("=" * 100)
print("FIGURE 6 — THREE FULL TEMPORAL GRAPHS PER SPECIES")
print("=" * 100)

for path in [MASTER, INCOME_FILE]:
    if not path.is_file():
        raise FileNotFoundError(path)

data = pd.read_excel(
    MASTER,
    sheet_name="FINAL_MASTER_2000_2026",
    engine="openpyxl",
    dtype={"Assembly_Accession": str}
)

if len(data) != 9157:
    raise ValueError(
        f"Expected 9157 genomes; found {len(data)}"
    )

required = [
    "Assembly_Accession",
    "Normalized_Species",
    "Country_Harmonized",
    "AMR_Genes",
]

for column in required:
    if column not in data.columns:
        raise ValueError(
            "Missing master column: " + column
        )

if (
    data["Assembly_Accession"].isna().any()
    or data["Assembly_Accession"].duplicated().any()
):
    raise ValueError(
        "Missing or duplicate genome accession"
    )

income_db = pd.read_excel(
    INCOME_FILE,
    sheet_name="COUNTRY_INCOME_COUNTS",
    engine="openpyxl"
)

for col in ["Country_Harmonized", "Income_Group"]:
    if col not in income_db.columns:
        raise ValueError(
            "Missing income database column: " + col
        )

income_db["Country_Key"] = (
    income_db["Country_Harmonized"]
    .fillna("")
    .astype(str)
    .str.strip()
)

income_db["Income_Group"] = (
    income_db["Income_Group"].map(
        normalize_income
    )
)

if income_db["Country_Key"].duplicated().any():
    raise ValueError(
        "Duplicate country labels in income database"
    )

income_lookup = income_db.set_index(
    "Country_Key"
)["Income_Group"]

data["Country_Key"] = (
    data["Country_Harmonized"]
    .fillna("")
    .astype(str)
    .str.strip()
)

data["Income"] = (
    data["Country_Key"]
    .map(income_lookup)
    .fillna("Unclassified")
)

data["Species"] = (
    data["Normalized_Species"]
    .fillna("")
    .astype(str)
    .str.strip()
)

# ============================================================
# YEAR AND AMR BURDEN
# ============================================================

YEAR_USED = find_year_column(data)

data["Year"] = data[YEAR_USED].map(
    year_from_cell
)

data["Period_Index"] = data["Year"].map(
    period_index
)

data["AMR_Genes_Per_Genome"] = (
    data["AMR_Genes"].map(
        lambda x: len(gene_set(x))
    )
)

if "AMRFinder_Status" in data.columns:

    status = (
        data["AMRFinder_Status"]
        .fillna("")
        .astype(str)
        .str.lower()
    )

    data["Annotation_Failed"] = status.str.contains(
        r"fail|error|not.run|not processed",
        regex=True
    )

else:
    data["Annotation_Failed"] = False

data["Eligible"] = (
    data["Species"].isin(SPECIES)
    & data["Income"].isin(INCOMES)
    & data["Period_Index"].notna()
    & ~data["Annotation_Failed"]
)

work = data[data["Eligible"]].copy()

if work.empty:
    raise ValueError(
        "No genomes eligible for temporal analysis"
    )

work["Period_Index"] = (
    work["Period_Index"].astype(int)
)

year_lower = YEAR_USED.lower()

if "collection" in year_lower or "isolation" in year_lower:
    year_type = "Collection/isolation year"
elif "release" in year_lower or "submi" in year_lower:
    year_type = "NCBI submission/release year"
else:
    year_type = "Year field not independently verified"

print("Master genomes:", len(data))
print("Eligible genomes:", len(work))
print("Year field:", YEAR_USED)
print("Year type:", year_type)

# ============================================================
# SUMMARY, KRUSKAL-WALLIS, PAIRWISE TESTS
# ============================================================

summary_rows = []
omnibus_rows = []
pairwise_rows = []

for species in SPECIES:

    for income in INCOMES:

        subset = work[
            work["Species"].eq(species)
            & work["Income"].eq(income)
        ]

        period_data = {}

        for period_i, period_name in enumerate(PERIOD_LABELS):

            arr = subset.loc[
                subset["Period_Index"].eq(period_i),
                "AMR_Genes_Per_Genome"
            ].to_numpy(dtype=float)

            period_data[period_i] = arr

            mean, sem = mean_sem(arr)

            summary_rows.append({
                "Species": species,
                "Income": income,
                "Period": period_name,
                "Period_Index": period_i,
                "N": len(arr),
                "Mean": mean,
                "SEM": sem,
                "Median": (
                    float(np.median(arr))
                    if len(arr) else np.nan
                ),
                "Q1": (
                    float(np.percentile(arr, 25))
                    if len(arr) else np.nan
                ),
                "Q3": (
                    float(np.percentile(arr, 75))
                    if len(arr) else np.nan
                ),
                "Eligible_For_Testing": (
                    len(arr) >= MIN_N
                )
            })

        available = [
            arr for arr in period_data.values()
            if len(arr) >= MIN_N
        ]

        if len(available) >= 2:
            try:
                h, p = kruskal(*available)
            except ValueError:
                h, p = 0.0, 1.0
        else:
            h, p = np.nan, np.nan

        omnibus_rows.append({
            "Species": species,
            "Income": income,
            "Total_Genomes": len(subset),
            "Testable_Periods": len(available),
            "Kruskal_Wallis_H": h,
            "Kruskal_Wallis_P_raw": p,
        })

        family = []

        for i, j in PLANNED_COMPARISONS:

            a = period_data[i]
            b = period_data[j]

            if len(a) < MIN_N or len(b) < MIN_N:
                continue

            p_raw = float(
                mannwhitneyu(
                    a, b,
                    alternative="two-sided"
                ).pvalue
            )

            family.append({
                "Species": species,
                "Income": income,
                "Period_1": PERIOD_LABELS[i],
                "Period_2": PERIOD_LABELS[j],
                "Index_1": i,
                "Index_2": j,
                "N_1": len(a),
                "N_2": len(b),
                "Mean_Difference": (
                    float(np.mean(b) - np.mean(a))
                ),
                "Pairwise_P_raw": p_raw
            })

        if family:

            corrected = multipletests(
                [
                    r["Pairwise_P_raw"]
                    for r in family
                ],
                method="holm"
            )[1]

            for row, adjusted in zip(
                family, corrected
            ):
                row["Pairwise_P_Holm"] = float(adjusted)
                pairwise_rows.append(row)

summary = pd.DataFrame(summary_rows)
omnibus = pd.DataFrame(omnibus_rows)

omnibus["Kruskal_Wallis_Q_BH"] = np.nan

valid = omnibus[
    "Kruskal_Wallis_P_raw"
].notna()

if valid.any():

    omnibus.loc[
        valid,
        "Kruskal_Wallis_Q_BH"
    ] = multipletests(
        omnibus.loc[
            valid,
            "Kruskal_Wallis_P_raw"
        ],
        method="fdr_bh"
    )[1]

pair_columns = [
    "Species",
    "Income",
    "Period_1",
    "Period_2",
    "Index_1",
    "Index_2",
    "N_1",
    "N_2",
    "Mean_Difference",
    "Pairwise_P_raw",
    "Pairwise_P_Holm",
]

pairwise = pd.DataFrame(
    pairwise_rows,
    columns=pair_columns
)

pairwise = pairwise.merge(
    omnibus[[
        "Species",
        "Income",
        "Kruskal_Wallis_Q_BH"
    ]],
    on=["Species", "Income"],
    how="left",
    validate="many_to_one"
)

pairwise["Significant"] = (
    (pairwise["Pairwise_P_Holm"] < .05)
    & (pairwise["Kruskal_Wallis_Q_BH"] < .05)
)

pairwise["Stars"] = [
    stars(p) if significant else "ns"
    for p, significant in zip(
        pairwise["Pairwise_P_Holm"],
        pairwise["Significant"]
    )
]

# ============================================================
# PLOT THREE STACKED GRAPHS FOR ONE SPECIES
# ============================================================

def draw_income_panel(ax, species, income):

    color = INCOME_COLORS[income]

    ss = summary[
        summary["Species"].eq(species)
        & summary["Income"].eq(income)
    ].sort_values("Period_Index")

    x = ss["Period_Index"].to_numpy(dtype=float)
    y = ss["Mean"].to_numpy(dtype=float)
    sem = ss["SEM"].to_numpy(dtype=float)
    n = ss["N"].to_numpy(dtype=int)

    sufficient = (
        (n >= MIN_N)
        & np.isfinite(y)
    )

    sparse = (
        (n > 0)
        & ~sufficient
        & np.isfinite(y)
    )

    # Mean ± SEM.
    ax.errorbar(
        x,
        np.where(sufficient, y, np.nan),
        yerr=np.where(
            sufficient,
            np.nan_to_num(sem, nan=0),
            np.nan
        ),
        color=color,
        linewidth=2.8,
        marker="o",
        markersize=8,
        markeredgecolor="white",
        markeredgewidth=1.0,
        capsize=5,
        capthick=1.5,
        elinewidth=1.7,
        zorder=5
    )

    if sparse.any():
        ax.scatter(
            x[sparse],
            y[sparse],
            marker="x",
            s=64,
            color=color,
            linewidth=1.8,
            zorder=7
        )

    # Graph statistic box.
    stat = omnibus[
        omnibus["Species"].eq(species)
        & omnibus["Income"].eq(income)
    ].iloc[0]

    h = stat["Kruskal_Wallis_H"]
    p = stat["Kruskal_Wallis_P_raw"]
    q = stat["Kruskal_Wallis_Q_BH"]

    if np.isfinite(h):
        stat_text = (
            f"Kruskal–Wallis H = {h:.3f}\n"
            f"p = {format_p(p)}\n"
            f"FDR q = {format_p(q)}"
        )
    else:
        stat_text = (
            "Kruskal–Wallis: not testable\n"
            "Fewer than two periods with n ≥ 5"
        )

    # Significant comparisons for this specific graph.
    tests = pairwise[
        pairwise["Species"].eq(species)
        & pairwise["Income"].eq(income)
        & pairwise["Significant"]
    ].sort_values(
        ["Pairwise_P_Holm", "Index_1", "Index_2"]
    ).head(MAX_BRACKETS)

    top_values = y + np.nan_to_num(
        sem, nan=0
    )

    finite = top_values[np.isfinite(top_values)]

    observed_max = (
        float(finite.max()) if len(finite) else 2.0
    )

    observed_max = max(2.0, observed_max)

    gap = max(0.8, observed_max * 0.14)

    bracket_start = observed_max + gap * 0.65

    # Reserve a large, clean region for brackets.
    upper_limit = (
        observed_max
        + gap * (len(tests) + 2.1)
    )

    upper_limit = max(
        upper_limit,
        observed_max * 1.48
    )

    ax.set_ylim(0, upper_limit)

    # Draw actual pairwise brackets.
    for level, row in enumerate(
        tests.itertuples(index=False)
    ):

        height = (
            bracket_start
            + level * gap
        )

        left = row.Index_1
        right = row.Index_2

        ax.plot(
            [
                left, left,
                right, right
            ],
            [
                height-gap*.16,
                height,
                height,
                height-gap*.16
            ],
            linewidth=1.35,
            color="#263747",
            zorder=11
        )

        ax.text(
            (left + right)/2,
            height + gap*.07,
            (
                f"{row.Stars}  "
                f"Holm p={format_p(row.Pairwise_P_Holm)}"
            ),
            ha="center",
            va="bottom",
            fontsize=10.3,
            fontweight="bold",
            color="#263747",
            zorder=12
        )

    # Exactly six period positions.
    ax.set_xlim(-.46, 5.46)

    # Each period label contains its own genome count.
    # This keeps n directly associated with its year,
    # unlike a distant row of tiny numbers.
    labels_with_n = [
        f"{label}\nn = {count:,}"
        for label, count in zip(
            PERIOD_LABELS, n
        )
    ]

    ax.set_xticks(range(6))
    ax.set_xticklabels(
        labels_with_n,
        fontsize=10.5,
        linespacing=1.65
    )

    ax.tick_params(
        axis="x",
        pad=12,
        length=3
    )

    ax.tick_params(
        axis="y",
        labelsize=10.5
    )

    ax.set_xlabel(
        "Year interval and number of genomes",
        fontsize=11,
        fontweight="medium",
        labelpad=13
    )

    ax.set_ylabel(
        "AMR genes per genome\n(mean ± SEM)",
        fontsize=11,
        fontweight="medium",
        labelpad=10
    )

    ax.yaxis.set_major_locator(
        MaxNLocator(nbins=6)
    )

    ax.grid(
        axis="y",
        color="#B6C2CD",
        linestyle=(0, (2, 3)),
        linewidth=.85,
        alpha=.42
    )

    ax.set_axisbelow(True)

    ax.spines["left"].set_color(color)
    ax.spines["left"].set_linewidth(2.3)

    ax.spines["bottom"].set_color("#728599")

    # Each graph has its own independent title.
    ax.set_title(
        INCOME_NAMES[income],
        loc="left",
        fontsize=15,
        fontweight="bold",
        color=color,
        pad=20
    )

    # Statistics placed within the right-upper
    # graph corner with a solid white background.
    ax.text(
        .986,
        .963,
        stat_text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10.1,
        color="#24384D",
        linespacing=1.55,
        bbox={
            "boxstyle": "round,pad=0.55",
            "facecolor": "white",
            "edgecolor": "#CBD6E0",
            "linewidth": .9,
            "alpha": .97
        },
        zorder=20
    )


# ============================================================
# GENERATE SEVEN LARGE THREE-PANEL FIGURES
# ============================================================

print()
print("Generating large species-specific figures...")

for species_number, species in enumerate(SPECIES, start=1):

    fig, axes = plt.subplots(
        3, 1,
        figsize=(17.0, 19.0),
        facecolor="white"
    )

    for ax, income in zip(axes, INCOMES):
        draw_income_panel(
            ax,
            species,
            income
        )

    species_n = int(
        work["Species"].eq(species).sum()
    )

    display_name = species.replace(
        "Enterobacter", "E."
    )

    fig.suptitle(
        display_name,
        fontsize=23,
        fontweight="bold",
        fontstyle="italic",
        color=SPECIES_COLORS[species],
        y=.987
    )

    fig.text(
        .5, .962,
        (
            "Income-stratified temporal AMR gene burden"
            f"  |  n = {species_n:,} genomes"
        ),
        ha="center",
        fontsize=13,
        color="#3B5167"
    )

    fig.text(
        .5,
        .030,
        (
            "Each bracket compares the two connected periods "
            "within the same income category. "
            "* p < 0.05; ** p < 0.01; "
            "*** p < 0.001; **** p < 0.0001 "
            "(Holm-adjusted)."
        ),
        ha="center",
        fontsize=9.8,
        color="#42576A"
    )

    fig.text(
        .5,
        .014,
        (
            f"Year field: {YEAR_USED} ({year_type}). "
            "Burden = distinct annotated AMR genes/genome; "
            "intrinsic-associated annotations included."
        ),
        ha="center",
        fontsize=9.1,
        color="#66798B"
    )

    # Large gap ensures statistics, x-axis values
    # and period-specific n remain visible.
    fig.subplots_adjust(
        left=.090,
        right=.967,
        top=.927,
        bottom=.078,
        hspace=.64
    )

    species_short = species.split()[-1]

    filename = (
        OUT
        / (
            f"Fig06_{species_number:02d}_"
            f"{species_short}_H_M_LLM_STACKED_600DPI.jpg"
        )
    )

    save_jpeg(fig, filename)

# ============================================================
# EXPORT FULL EXCEL
# ============================================================

genome_columns = [
    "Assembly_Accession",
    "Species",
    "Country_Harmonized",
    "Income",
    "Year",
    "Period_Index",
    "AMR_Genes_Per_Genome",
    "AMR_Genes",
    "Annotation_Failed",
    "Eligible"
]

genome_export = data[
    genome_columns
].copy()

genome_export["Period"] = (
    genome_export["Period_Index"].map(
        lambda x: (
            PERIOD_LABELS[int(x)]
            if pd.notna(x)
            else ""
        )
    )
)

qc = pd.DataFrame([
    ["Total master genomes", len(data)],
    ["Eligible temporal genomes", len(work)],
    ["Seven species", len(SPECIES)],
    ["Income categories", len(INCOMES)],
    ["Time periods", len(PERIOD_LABELS)],
    ["Year column", YEAR_USED],
    ["Year field type", year_type],
    ["Minimum n for statistical testing", MIN_N],
    [
        "Omnibus test",
        "Kruskal-Wallis within each species × income"
    ],
    [
        "Omnibus correction",
        "Benjamini-Hochberg FDR across 21 tests"
    ],
    [
        "Pairwise test",
        "Two-sided Mann-Whitney U"
    ],
    [
        "Pairwise correction",
        "Holm within each species × income"
    ],
    [
        "Pairwise comparisons",
        "Five adjacent periods and first versus last"
    ],
    [
        "Figure significance",
        "Only comparisons passing omnibus q and Holm p < 0.05"
    ],
    [
        "AMR burden",
        "Number of distinct AMR gene names per genome"
    ],
    [
        "Interpretation",
        "Annotated genomic AMR burden, not phenotypic resistance"
    ],
    [
        "Data limitation",
        "Intrinsic resistance-associated genes are included"
    ],
    [
        "Sampling limitation",
        "Public genome sampling is not population-representative"
    ],
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
        sheet_name="SPECIES_INCOME_PERIOD",
        index=False
    )

    omnibus.to_excel(
        writer,
        sheet_name="KRUSKAL_WALLIS_P_Q",
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

    excel_style(writer)

print()
print("=" * 100)
print("FIGURE 6 COMPLETED")
print("=" * 100)
print("Master genomes:", len(data))
print("Eligible genomes:", len(work))
print("Year column:", YEAR_USED)
print("Year type:", year_type)
print("Separate JPEG figures:", 7)
print("Omnibus tests:", int(valid.sum()))
print(
    "Significant omnibus comparisons:",
    int(
        (omnibus["Kruskal_Wallis_Q_BH"] < .05).sum()
    )
)
print(
    "Significant pairwise comparisons:",
    int(pairwise["Significant"].sum())
)
print("Output folder:", OUT)
print("Excel:", EXCEL)
print("=" * 100)
