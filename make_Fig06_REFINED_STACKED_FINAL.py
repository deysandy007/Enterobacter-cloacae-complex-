
#!/usr/bin/env python3

# ================================================================
# FIGURE 6 — REFINED INCOME-STRATIFIED TEMPORAL AMR BURDEN
#
# Seven separate species figures
# Three vertically stacked income panels in each figure
#
# IMPROVEMENTS:
# - Narrower canvas: 12.5 x 17.5 inches
# - Full species name ONLY at figure top
# - Kruskal-Wallis statistics in separate right column
# - No statistics positioned over plotted lines
# - Full five-year labels and sample sizes on every graph
# - Significance brackets with dedicated vertical space
# - Star significance legend moved to bottom footer
# - 600-DPI JPEG and complete statistical Excel
#
# Data: final 9,157-genome ECC master
# ================================================================

from pathlib import Path
import re
import gc

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.gridspec import GridSpec

from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter

# ================================================================
# PATHS AND SETTINGS
# ================================================================

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
    / "REFINED_THREE_STACKED"
)

OUT.mkdir(parents=True, exist_ok=True)

XLSX = OUT / "Fig06_Refined_Temporal_Statistics.xlsx"

DPI = 600
MIN_N = 5
MAX_BRACKETS = 3

# Force collection year if needed; otherwise identify automatically.
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
    f"{a}\u2013{b}" for a, b in PERIODS
]

# Predefined temporal comparisons:
# Five consecutive comparisons, plus first vs final.
COMPARISONS = [
    (0, 1), (1, 2), (2, 3),
    (3, 4), (4, 5), (0, 5)
]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
})

# ================================================================
# HELPER FUNCTIONS
# ================================================================

def parse_year(value):
    if pd.isna(value):
        return np.nan
    if isinstance(value, pd.Timestamp):
        return int(value.year)

    text = str(value).strip()

    if re.fullmatch(r"\d{4}(?:\.0+)?", text):
        return int(float(text))

    m = re.search(
        r"(?<!\d)(19\d{2}|20\d{2})(?!\d)",
        text
    )
    return int(m.group(1)) if m else np.nan


def find_year_column(frame):
    if YEAR_COLUMN is not None:
        if YEAR_COLUMN not in frame.columns:
            raise ValueError(
                f"Year column absent: {YEAR_COLUMN}"
            )
        return YEAR_COLUMN

    lookup = {str(c).lower(): c for c in frame.columns}

    for name in YEAR_CANDIDATES:
        key = name.lower()
        if key not in lookup:
            continue
        col = lookup[key]
        years = frame[col].map(parse_year)
        if years.between(2000, 2026).any():
            return col

    raise ValueError(
        "No valid year column found. Available columns:\n"
        + "\n".join(map(str, frame.columns))
    )


def period_index(year):
    if pd.isna(year):
        return np.nan
    for i, (first, last) in enumerate(PERIODS):
        if first <= int(year) <= last:
            return i
    return np.nan


def gene_set(value):
    if pd.isna(value):
        return set()
    text = str(value).strip()

    if text.lower() in {
        "", "nan", "na", "n/a", "none",
        "-", "no hits", "not detected"
    }:
        return set()

    return {
        s.strip()
        for s in text.split(";")
        if s.strip()
    }


def normalize_income(value):
    text = re.sub(
        r"[\s_\-]+", "",
        str(value).strip().upper()
    )

    lookup = {
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
    return lookup.get(text, "Unclassified")


def mean_sem(values):
    a = np.asarray(values, dtype=float)
    n = len(a)

    if n == 0:
        return np.nan, np.nan
    if n == 1:
        return float(a.mean()), np.nan

    return (
        float(a.mean()),
        float(a.std(ddof=1) / np.sqrt(n))
    )


def p_label(p):
    if not np.isfinite(p):
        return "not estimable"
    if p < 0.0001:
        return f"{p:.2e}"
    return f"{p:.4f}"


def significance_stars(p):
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


def save_jpg(fig, filepath):
    fig.savefig(
        filepath,
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
    print("[SAVED]", filepath)


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

        for j, cell in enumerate(ws[1], 1):
            ws.column_dimensions[
                get_column_letter(j)
            ].width = min(
                48,
                max(15, len(str(cell.value)) + 3)
            )


# ================================================================
# READ GENOME AND INCOME DATA
# ================================================================

print("=" * 90)
print("FIGURE 6 — REFINED SEVEN SPECIES TEMPORAL ANALYSIS")
print("=" * 90)

for f in (MASTER, INCOME_FILE):
    if not f.is_file():
        raise FileNotFoundError(f)

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

needed = [
    "Assembly_Accession",
    "Normalized_Species",
    "Country_Harmonized",
    "AMR_Genes",
]

for c in needed:
    if c not in df.columns:
        raise ValueError(f"Missing master column: {c}")

if (
    df["Assembly_Accession"].isna().any()
    or df["Assembly_Accession"].duplicated().any()
):
    raise ValueError(
        "Genome accessions missing or duplicated."
    )

income = pd.read_excel(
    INCOME_FILE,
    sheet_name="COUNTRY_INCOME_COUNTS",
    engine="openpyxl"
)

for c in ("Country_Harmonized", "Income_Group"):
    if c not in income.columns:
        raise ValueError(f"Missing income column: {c}")

income["Country_Key"] = (
    income["Country_Harmonized"]
    .fillna("").astype(str).str.strip()
)

income["Income_Group"] = (
    income["Income_Group"].map(normalize_income)
)

if income["Country_Key"].duplicated().any():
    raise ValueError(
        "Duplicate country in frozen income database."
    )

df["Country_Key"] = (
    df["Country_Harmonized"]
    .fillna("").astype(str).str.strip()
)

income_lookup = (
    income.set_index("Country_Key")["Income_Group"]
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

df["Year"] = df[YEAR_USED].map(parse_year)
df["Period_Index"] = df["Year"].map(period_index)

df["AMR_Genes_Per_Genome"] = (
    df["AMR_Genes"].map(
        lambda value: len(gene_set(value))
    )
)

if "AMRFinder_Status" in df.columns:
    status = (
        df["AMRFinder_Status"]
        .fillna("").astype(str).str.lower()
    )
    df["Annotation_Failed"] = status.str.contains(
        r"fail|error|not.run|not processed",
        regex=True
    )
else:
    df["Annotation_Failed"] = False

df["Eligible"] = (
    df["Species"].isin(SPECIES)
    & df["Income"].isin(INCOMES)
    & df["Period_Index"].notna()
    & ~df["Annotation_Failed"]
)

work = df.loc[df["Eligible"]].copy()

if work.empty:
    raise ValueError("No eligible genomes.")

work["Period_Index"] = (
    work["Period_Index"].astype(int)
)

year_lower = YEAR_USED.lower()

if "collection" in year_lower or "isolation" in year_lower:
    YEAR_TYPE = "Collection/isolation year"
elif "release" in year_lower or "submi" in year_lower:
    YEAR_TYPE = "Submission/release year"
else:
    YEAR_TYPE = "Year type not independently verified"

print("Master genomes:", len(df))
print("Eligible genomes:", len(work))
print("Year column:", YEAR_USED)
print("Year interpretation:", YEAR_TYPE)

# ================================================================
# SUMMARY AND STATISTICAL TESTS
# ================================================================

summary_rows = []
omnibus_rows = []
pairwise_rows = []

for species in SPECIES:
    for inc in INCOMES:

        sub = work[
            work["Species"].eq(species)
            & work["Income"].eq(inc)
        ]

        arrays = {}

        for k, label in enumerate(PERIOD_LABELS):
            a = sub.loc[
                sub["Period_Index"].eq(k),
                "AMR_Genes_Per_Genome"
            ].to_numpy(float)

            arrays[k] = a
            avg, sem = mean_sem(a)

            summary_rows.append({
                "Species": species,
                "Income": inc,
                "Period": label,
                "Period_Index": k,
                "N": len(a),
                "Mean": avg,
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
                ),
            })

        testable = [
            a for a in arrays.values()
            if len(a) >= MIN_N
        ]

        if len(testable) >= 2:
            try:
                h, p = kruskal(*testable)
            except ValueError:
                h, p = 0.0, 1.0
        else:
            h, p = np.nan, np.nan

        omnibus_rows.append({
            "Species": species,
            "Income": inc,
            "Total_N": len(sub),
            "Testable_Periods": len(testable),
            "Kruskal_Wallis_H": h,
            "P_raw": p,
        })

        family = []

        for i, j in COMPARISONS:
            a, b = arrays[i], arrays[j]

            if len(a) < MIN_N or len(b) < MIN_N:
                continue

            raw_p = float(
                mannwhitneyu(
                    a, b,
                    alternative="two-sided"
                ).pvalue
            )

            family.append({
                "Species": species,
                "Income": inc,
                "Period_1": PERIOD_LABELS[i],
                "Period_2": PERIOD_LABELS[j],
                "Index_1": i,
                "Index_2": j,
                "N_1": len(a),
                "N_2": len(b),
                "Mean_Difference": float(
                    b.mean() - a.mean()
                ),
                "P_raw": raw_p
            })

        if family:
            corrected = multipletests(
                [r["P_raw"] for r in family],
                method="holm"
            )[1]

            for record, p_holm in zip(
                family, corrected
            ):
                record["P_Holm"] = float(p_holm)
                pairwise_rows.append(record)

summary = pd.DataFrame(summary_rows)
omnibus = pd.DataFrame(omnibus_rows)

omnibus["Q_BH"] = np.nan
valid = omnibus["P_raw"].notna()

if valid.any():
    omnibus.loc[valid, "Q_BH"] = multipletests(
        omnibus.loc[valid, "P_raw"],
        method="fdr_bh"
    )[1]

pair_cols = [
    "Species", "Income", "Period_1",
    "Period_2", "Index_1", "Index_2",
    "N_1", "N_2", "Mean_Difference",
    "P_raw", "P_Holm",
]

pairwise = pd.DataFrame(
    pairwise_rows,
    columns=pair_cols
)

pairwise = pairwise.merge(
    omnibus[["Species", "Income", "Q_BH"]],
    on=["Species", "Income"],
    how="left",
    validate="many_to_one"
)

pairwise["Significant"] = (
    (pairwise["P_Holm"] < 0.05)
    & (pairwise["Q_BH"] < 0.05)
)

pairwise["Stars"] = [
    significance_stars(p) if ok else "ns"
    for p, ok in zip(
        pairwise["P_Holm"],
        pairwise["Significant"]
    )
]

# ================================================================
# PLOT: LARGER DATA REGION, SEPARATE STATISTICS COLUMN
# ================================================================

def draw_panel(ax, stats_ax, species, inc):

    color = INCOME_COLORS[inc]

    stats_ax.axis("off")

    sub = summary[
        summary["Species"].eq(species)
        & summary["Income"].eq(inc)
    ].sort_values("Period_Index")

    x = sub["Period_Index"].to_numpy(float)
    y = sub["Mean"].to_numpy(float)
    sem = sub["SEM"].to_numpy(float)
    n = sub["N"].to_numpy(int)

    enough = (n >= MIN_N) & np.isfinite(y)
    sparse = (n > 0) & (~enough) & np.isfinite(y)

    # Data display, using original measurements.
    ax.errorbar(
        x,
        np.where(enough, y, np.nan),
        yerr=np.where(
            enough,
            np.nan_to_num(sem, nan=0),
            np.nan
        ),
        color=color,
        marker="o",
        markersize=7,
        linewidth=2.5,
        capsize=4.5,
        elinewidth=1.6,
        markeredgecolor="white",
        markeredgewidth=0.8,
        zorder=5
    )

    if sparse.any():
        ax.scatter(
            x[sparse],
            y[sparse],
            marker="x",
            s=55,
            linewidth=1.7,
            color=color,
            zorder=6
        )

    tests = pairwise[
        pairwise["Species"].eq(species)
        & pairwise["Income"].eq(inc)
        & pairwise["Significant"]
    ].sort_values(
        ["P_Holm", "Index_1", "Index_2"]
    ).head(MAX_BRACKETS)

    max_data = y + np.nan_to_num(
        sem, nan=0
    )

    finite = max_data[np.isfinite(max_data)]

    observed_top = max(
        2.0,
        float(np.max(finite)) if len(finite) else 2.0
    )

    step = max(.8, .13 * observed_top)

    bracket_first = observed_top + step * .75

    ax.set_ylim(
        0,
        max(
            observed_top * 1.33,
            bracket_first + (len(tests)+.9)*step
        )
    )

    # Brackets sit ABOVE data, not through the points.
    for level, r in enumerate(
        tests.itertuples(index=False)
    ):
        yy = bracket_first + level*step

        ax.plot(
            [r.Index_1, r.Index_1, r.Index_2, r.Index_2],
            [yy-step*.14, yy, yy, yy-step*.14],
            color="#27394B",
            linewidth=1.15,
            clip_on=False,
            zorder=9
        )

        ax.text(
            (r.Index_1+r.Index_2)/2,
            yy+step*.06,
            f"{r.Stars}  p={p_label(r.P_Holm)}",
            ha="center",
            va="bottom",
            fontsize=9.3,
            fontweight="bold",
            color="#27394B",
            zorder=10
        )

    # Full years plus individual n on EVERY graph.
    ticklabels = [
        f"{label}\nn = {count:,}"
        for label, count in zip(PERIOD_LABELS, n)
    ]

    ax.set_xticks(range(6))
    ax.set_xticklabels(
        ticklabels,
        fontsize=9.1,
        linespacing=1.5
    )

    ax.set_xlim(-.4, 5.4)

    ax.set_xlabel(
        "Five-year interval",
        fontsize=10,
        labelpad=10
    )

    ax.set_ylabel(
        "AMR genes per genome\n(mean ± SEM)",
        fontsize=10.2,
        labelpad=12
    )

    ax.set_title(
        INCOME_NAMES[inc],
        loc="left",
        fontsize=13,
        fontweight="bold",
        color=color,
        pad=14
    )

    ax.yaxis.set_major_locator(
        MaxNLocator(nbins=5)
    )

    ax.grid(
        axis="y",
        linestyle=":",
        linewidth=.7,
        alpha=.28
    )

    ax.spines["left"].set_color(color)
    ax.spines["left"].set_linewidth(2.0)
    ax.spines["bottom"].set_color("#A0ACB9")
    ax.tick_params(axis="y", labelsize=9.5)

    # Statistics go OUTSIDE plotting axes.
    row = omnibus[
        omnibus["Species"].eq(species)
        & omnibus["Income"].eq(inc)
    ].iloc[0]

    h, p, q = (
        row["Kruskal_Wallis_H"],
        row["P_raw"],
        row["Q_BH"]
    )

    stats_ax.text(
        .07, .90,
        "Kruskal–Wallis",
        ha="left",
        va="top",
        fontsize=11,
        fontweight="bold",
        color="#173352",
        transform=stats_ax.transAxes
    )

    if np.isfinite(h):
        stat_text = (
            f"H = {h:.3f}\n"
            f"p = {p_label(p)}\n"
            f"FDR q = {p_label(q)}"
        )
    else:
        stat_text = (
            "Not testable\n"
            "Fewer than two\n"
            "eligible periods"
        )

    stats_ax.text(
        .07, .73,
        stat_text,
        ha="left",
        va="top",
        fontsize=10.2,
        linespacing=1.9,
        color="#344A60",
        transform=stats_ax.transAxes
    )

    stats_ax.text(
        .07, .22,
        f"Total n = {int(row['Total_N']):,}",
        ha="left",
        va="top",
        fontsize=9.7,
        color="#50667A",
        transform=stats_ax.transAxes
    )


# ================================================================
# SEVEN SPECIES FIGURES
# ================================================================

for idx, species in enumerate(SPECIES, start=1):

    fig = plt.figure(
        figsize=(12.5, 17.5),
        facecolor="white"
    )

    # Three rows, each containing chart + separate stats area.
    grid = GridSpec(
        3, 2,
        figure=fig,
        width_ratios=[4.8, 1.35],
        left=.11,
        right=.97,
        top=.91,
        bottom=.125,
        hspace=.65,
        wspace=.17
    )

    for row_idx, inc in enumerate(INCOMES):
        ax = fig.add_subplot(grid[row_idx, 0])
        stat_ax = fig.add_subplot(grid[row_idx, 1])

        draw_panel(
            ax,
            stat_ax,
            species,
            inc
        )

    # ONLY full species name at the top.
    fig.suptitle(
        species,
        fontsize=21,
        fontweight="bold",
        fontstyle="italic",
        color=SPECIES_COLORS[species],
        y=.980
    )

    # Star explanation in a DEDICATED bottom footer.
    fig.text(
        .5, .075,
        "Significance: * p < 0.05   "
        "** p < 0.01   "
        "*** p < 0.001   "
        "**** p < 0.0001",
        ha="center",
        fontsize=10.4,
        color="#26394C"
    )

    fig.text(
        .5, .056,
        "Brackets identify comparisons between the connected "
        "time intervals; reported bracket p-values are Holm-adjusted.",
        ha="center",
        fontsize=9.0,
        color="#53697C"
    )

    fig.text(
        .5, .037,
        "Omnibus q: Benjamini–Hochberg FDR across 21 "
        "species–income tests. Sparse periods (n < 5) "
        "are excluded from significance testing.",
        ha="center",
        fontsize=8.6,
        color="#607385"
    )

    fig.text(
        .5, .019,
        f"Time field: {YEAR_USED} ({YEAR_TYPE}). "
        "Burden includes all distinct annotated AMR genes.",
        ha="center",
        fontsize=8.1,
        color="#65788A"
    )

    output = OUT / (
        f"Fig06_{idx:02d}_"
        f"{species.split()[-1]}_REFINED_600DPI.jpg"
    )

    save_jpg(fig, output)

# ================================================================
# EXPORT RESULTS
# ================================================================

genomes = df[[
    "Assembly_Accession",
    "Species",
    "Country_Harmonized",
    "Income",
    "Year",
    "Period_Index",
    "AMR_Genes_Per_Genome",
    "AMR_Genes",
    "Annotation_Failed",
    "Eligible",
]].copy()

genomes["Period"] = genomes["Period_Index"].map(
    lambda x: (
        PERIOD_LABELS[int(x)]
        if pd.notna(x) else ""
    )
)

methods = pd.DataFrame([
    ["Final master genomes", len(df)],
    ["Eligible temporal genomes", len(work)],
    ["Species", 7],
    ["Income groups", "H, M, L+LM"],
    ["Year column", YEAR_USED],
    ["Year field type", YEAR_TYPE],
    ["AMR burden", "Distinct annotated gene names per genome"],
    ["Error bars", "Mean +/- standard error of mean"],
    ["Minimum group n for tests", MIN_N],
    ["Omnibus", "Kruskal-Wallis within species-income"],
    ["Omnibus adjustment", "BH FDR across 21 tests"],
    ["Pairwise", "Mann-Whitney U, predefined temporal pairs"],
    ["Pairwise adjustment", "Holm within each species-income"],
    ["Bracket criterion", "Omnibus q < 0.05 and pairwise Holm p < 0.05"],
    ["Maximum brackets shown", MAX_BRACKETS],
    ["Interpretation", "Sampled-genome AMR burden; intrinsic genes included"],
    ["Caveat", "Country and lineage clustering are not modeled"],
], columns=["Parameter", "Value"])

with pd.ExcelWriter(
    XLSX,
    engine="openpyxl"
) as writer:

    genomes.to_excel(
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

    methods.to_excel(
        writer,
        sheet_name="METHODS_QC",
        index=False
    )

    excel_style(writer)

print()
print("=" * 90)
print("FIGURE 6 — REFINED LAYOUT COMPLETED")
print("=" * 90)
print("Master genomes:", len(df))
print("Eligible genomes:", len(work))
print("Year field:", YEAR_USED)
print("Individual species JPEGs:", 7)
print("Output folder:", OUT)
print("Statistics Excel:", XLSX)
print("=" * 90)
