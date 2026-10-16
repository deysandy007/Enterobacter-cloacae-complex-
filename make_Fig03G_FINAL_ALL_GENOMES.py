
from pathlib import Path
from itertools import combinations
from datetime import datetime
import gc

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ============================================================
# FIGURE 3G - FINAL ALL-GENOME ONE HEALTH AMR BURDEN
#
# NO COUNTRY ELIGIBILITY THRESHOLD
#
# INCLUSION:
# Every genome from seven resolved ECC species
# with H / M / L+LM income and
# Human / Animal / Environment source.
#
# TESTS:
# Human vs Animal
# Human vs Environment
# Animal vs Environment
#
# Test if BOTH groups have at least one genome.
# Empty comparisons are excluded.
#
# Two-sided genome-level Mann-Whitney U.
# BH FDR across ALL performed comparisons.
#
# EFFECT:
# Cliff's delta, comparison minus reference
# Difference in median AMR burden
# Bootstrap 95% CI of mean difference
#
# OUTPUT:
# Seven 600-DPI JPEG figures
# One comprehensive Excel workbook
#
# Note:
# Genome-level P-values are unadjusted for
# country, BioProject and lineage clustering.
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

EXCEL = OUT / "Fig03G_FINAL_ALL_GENOMES_STATISTICS.xlsx"

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
PAIRS = list(combinations(SOURCES, 2))

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

INCOME_NAMES = {
    "H": "High income",
    "M": "Upper-middle income",
    "L+LM": "Low + lower-middle income"
}

SOURCE_POSITIONS = {
    "Human": 1,
    "Animal": 2,
    "Environment": 3
}

BOOTSTRAPS = 3000
MAX_SCATTER_POINTS = 350
DPI = 600

RNG = np.random.default_rng(20261010)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})

# ============================================================
# FUNCTIONS
# ============================================================

def fdr_bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    q = np.full(len(p), np.nan)
    valid = np.flatnonzero(np.isfinite(p))

    if len(valid) == 0:
        return q

    order = valid[np.argsort(p[valid])]
    n = len(order)
    adj = p[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    q[order] = np.clip(adj, 0, 1)
    return q


def sig_label(q):
    if not np.isfinite(q):
        return "EMPTY"
    if q < 0.0001:
        return "****"
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def cliff_delta(a, b):
    # Positive: comparison group tends to be higher.
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    if not len(a) or not len(b):
        return np.nan

    # Efficient calculation without generating
    # a large pairwise difference matrix.
    a_sorted = np.sort(a)
    b_sorted = np.sort(b)

    fewer = np.searchsorted(
        a_sorted, b_sorted, side="left"
    ).sum()

    greater = (
        len(a_sorted)
        - np.searchsorted(
            a_sorted, b_sorted, side="right"
        )
    ).sum()

    return float(
        (fewer - greater) / (len(a) * len(b))
    )


def bootstrap_mean_difference(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    # Memory-safe iterative bootstrap.
    effects = np.empty(BOOTSTRAPS, dtype=float)

    for i in range(BOOTSTRAPS):
        sampled_a = RNG.choice(
            a, size=len(a), replace=True
        )
        sampled_b = RNG.choice(
            b, size=len(b), replace=True
        )
        effects[i] = (
            sampled_b.mean() - sampled_a.mean()
        )

    return (
        float(np.quantile(effects, .025)),
        float(np.quantile(effects, .975))
    )


def short_species(name):
    return name.replace("Enterobacter ", "E. ")


def species_tag(name):
    return name.replace("Enterobacter ", "E_")


# ============================================================
# LOAD MASTER
# ============================================================

print("=" * 110)
print("FIGURE 3G - FINAL GENOME-LEVEL AMR BURDEN ANALYSIS")
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

if len(df) != 9157:
    raise ValueError(
        f"Expected 9157 genomes; found {len(df)}"
    )

required = [
    "Assembly_Accession",
    "Normalized_Species",
    "Final_Source_Category",
    "Country_Harmonized",
    "Unique_AMR_Gene_Count",
    "AMRFinder_Status"
]

missing = [c for c in required if c not in df.columns]

if missing:
    raise ValueError("Missing master fields: " + str(missing))

if (
    df["Assembly_Accession"].isna().any()
    or df["Assembly_Accession"].duplicated().any()
):
    raise ValueError("Invalid accession identifiers")

status = (
    df["AMRFinder_Status"]
    .fillna("")
    .astype(str)
    .str.strip()
    .str.upper()
)

if not status.eq("COMPLETE").all():
    raise ValueError("Incomplete AMRFinder results")

df["Species"] = (
    df["Normalized_Species"]
    .fillna("").astype(str).str.strip()
)

df["Source"] = (
    df["Final_Source_Category"]
    .fillna("").astype(str).str.strip()
)

df["Country"] = (
    df["Country_Harmonized"]
    .fillna("").astype(str).str.strip()
)

df["AMR"] = pd.to_numeric(
    df["Unique_AMR_Gene_Count"],
    errors="coerce"
)

if (
    df["AMR"].isna().any()
    or (df["AMR"] < 0).any()
    or (df["AMR"] % 1 != 0).any()
):
    raise ValueError("Invalid AMR gene counts")

# ============================================================
# INCOME CLASSIFICATION
# ============================================================

income = pd.read_excel(
    INCOME_FILE,
    sheet_name="COUNTRY_INCOME_COUNTS",
    engine="openpyxl"
)

for c in [
    "Country_Harmonized",
    "Income_Group",
    "Genome_Count"
]:
    if c not in income.columns:
        raise ValueError("Missing income column: " + c)

income["Country"] = (
    income["Country_Harmonized"]
    .fillna("").astype(str).str.strip()
)

income["Income_Group"] = (
    income["Income_Group"]
    .fillna("Unclassified")
    .astype(str).str.strip()
)

if income["Country"].duplicated().any():
    raise ValueError("Duplicate countries in income map")

master_counts = (
    df["Country"].value_counts().rename("Master")
)

income_counts = (
    income.set_index("Country")[
        "Genome_Count"
    ].rename("Reference")
)

check = pd.concat(
    [master_counts, income_counts],
    axis=1
)

if (
    check.isna().any().any()
    or not (
        check["Master"].astype(int)
        == check["Reference"].astype(int)
    ).all()
):
    raise ValueError(
        "Frozen country-income map does not "
        "match the current master"
    )

df = df.merge(
    income[["Country", "Income_Group"]],
    on="Country",
    how="left",
    validate="many_to_one"
)

if df["Income_Group"].isna().any():
    raise ValueError("Income merge failed")

# ============================================================
# FILTER
# ============================================================

data = df[
    df["Species"].isin(SPECIES)
    & df["Source"].isin(SOURCES)
    & df["Income_Group"].isin(INCOMES)
].copy()

if data.empty:
    raise ValueError("No genomes meet analysis filters")

print("Master genomes:", len(df))
print("Genomes used:", len(data))

# ============================================================
# SUMMARY
# ============================================================

full_groups = pd.MultiIndex.from_product(
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
        Mean_AMR=("AMR", "mean"),
        Median_AMR=("AMR", "median"),
        SD_AMR=("AMR", "std"),
        Q1_AMR=("AMR", lambda x: x.quantile(.25)),
        Q3_AMR=("AMR", lambda x: x.quantile(.75))
    )
    .reindex(full_groups)
    .reset_index()
)

summary["Genome_N"] = (
    summary["Genome_N"].fillna(0).astype(int)
)

summary["Country_N"] = (
    summary["Country_N"].fillna(0).astype(int)
)

if summary["Genome_N"].sum() != len(data):
    raise RuntimeError("Summary counts do not match")

# ============================================================
# 63 POTENTIAL PAIRWISE COMPARISONS
# ============================================================

print("\nRunning Mann-Whitney U comparisons...")

results = []

for species in SPECIES:

    for inc in INCOMES:

        sub = data[
            data["Species"].eq(species)
            & data["Income_Group"].eq(inc)
        ]

        arrays = {
            source: sub.loc[
                sub["Source"].eq(source),
                "AMR"
            ].to_numpy(dtype=float)
            for source in SOURCES
        }

        countries = {
            source: sub.loc[
                sub["Source"].eq(source),
                "Country"
            ].nunique()
            for source in SOURCES
        }

        for source_a, source_b in PAIRS:

            a = arrays[source_a]
            b = arrays[source_b]

            n_a, n_b = len(a), len(b)

            row = {
                "Species": species,
                "Income_Group": inc,
                "Source_A": source_a,
                "Source_B": source_b,
                "N_A": n_a,
                "N_B": n_b,
                "Countries_A": countries[source_a],
                "Countries_B": countries[source_b],
                "Mean_A": (
                    float(a.mean()) if n_a else np.nan
                ),
                "Mean_B": (
                    float(b.mean()) if n_b else np.nan
                ),
                "Median_A": (
                    float(np.median(a))
                    if n_a else np.nan
                ),
                "Median_B": (
                    float(np.median(b))
                    if n_b else np.nan
                ),
                "Mean_Difference_B_minus_A": np.nan,
                "Median_Difference_B_minus_A": np.nan,
                "Cliffs_Delta": np.nan,
                "CI95_MeanDiff_Lower": np.nan,
                "CI95_MeanDiff_Upper": np.nan,
                "Mann_Whitney_U": np.nan,
                "Raw_P": np.nan,
                "FDR": np.nan,
                "Significance": "EMPTY",
                "Status": (
                    "TESTED"
                    if n_a > 0 and n_b > 0
                    else "EMPTY_GROUP"
                )
            }

            if n_a > 0 and n_b > 0:

                if (
                    np.all(a == a[0])
                    and np.all(b == b[0])
                    and a[0] == b[0]
                ):
                    u = n_a * n_b / 2
                    p = 1.0

                else:
                    test = mannwhitneyu(
                        a,
                        b,
                        alternative="two-sided",
                        method="auto"
                    )
                    u = float(test.statistic)
                    p = float(test.pvalue)

                ci_low, ci_high = (
                    bootstrap_mean_difference(a, b)
                )

                row.update({
                    "Mean_Difference_B_minus_A": (
                        float(b.mean() - a.mean())
                    ),
                    "Median_Difference_B_minus_A": (
                        float(np.median(b) - np.median(a))
                    ),
                    "Cliffs_Delta": cliff_delta(a, b),
                    "CI95_MeanDiff_Lower": ci_low,
                    "CI95_MeanDiff_Upper": ci_high,
                    "Mann_Whitney_U": u,
                    "Raw_P": p
                })

            results.append(row)

tests = pd.DataFrame(results)

if len(tests) != 63:
    raise RuntimeError(
        f"Expected 63 comparisons, found {len(tests)}"
    )

tests["FDR"] = fdr_bh(
    tests["Raw_P"].to_numpy()
)

tests["Significance"] = (
    tests["FDR"].apply(sig_label)
)

tested = tests[
    tests["Status"] == "TESTED"
].copy()

print("Total possible comparisons:", len(tests))
print("Tests performed:", len(tested))
print(
    "FDR-significant:",
    int((tested["FDR"] < .05).sum())
)
print(
    "Excluded empty comparisons:",
    int(tests["Status"].eq("EMPTY_GROUP").sum())
)

# ============================================================
# SAVE EXCEL BEFORE FIGURES
# ============================================================

print("\nSaving Excel...")

methods = pd.DataFrame([
    ("Primary outcome", "Unique_AMR_Gene_Count per genome"),
    ("Species", "All seven resolved ECC species"),
    ("Sources", "Human; Animal; Environment"),
    ("Income", "H; M; L+LM"),
    ("Exclusion", "Only empty source comparisons"),
    ("Primary test", "Two-sided Mann-Whitney U on genomes"),
    ("Independence assumption", "Genomes treated as independent observations"),
    ("Multiple testing", "BH across all nonempty comparisons"),
    ("Effect", "Cliff's delta, positive means group B tends higher"),
    ("Additional effect", "Mean and median burden difference: B minus A"),
    ("CI", "Genome-level bootstrap 95% CI for mean difference"),
    ("Important limitation", "Unadjusted for country and BioProject clustering"),
    ("Outcome interpretation", "Annotated AMR gene repertoire, not phenotype"),
    ("Generated", datetime.now().isoformat())
], columns=["Parameter", "Description"])

with pd.ExcelWriter(
    EXCEL,
    engine="openpyxl"
) as writer:

    tests.to_excel(
        writer,
        sheet_name="ALL_63_COMPARISONS",
        index=False
    )

    tested.to_excel(
        writer,
        sheet_name="TESTED_ONLY",
        index=False
    )

    summary.to_excel(
        writer,
        sheet_name="SOURCE_SUMMARY",
        index=False
    )

    data[[
        "Assembly_Accession",
        "Species",
        "Income_Group",
        "Source",
        "Country",
        "AMR"
    ]].to_excel(
        writer,
        sheet_name="GENOME_AMR_DATA",
        index=False
    )

    methods.to_excel(
        writer,
        sheet_name="METHODS",
        index=False
    )

print("Excel saved:", EXCEL)

# ============================================================
# CREATE SEVEN SIGNIFICANCE-ANNOTATED FIGURES
# ============================================================

print("\nGenerating seven species figures...")

global_max = float(data["AMR"].max())
ymax = global_max + 17

brackets = [
    ("Human", "Animal", global_max + 6),
    ("Animal", "Environment", global_max + 10),
    ("Human", "Environment", global_max + 14)
]

for species in SPECIES:

    sub_species = data[
        data["Species"] == species
    ]

    fig, axes = plt.subplots(
        1, 3,
        figsize=(18, 8),
        sharey=True,
        facecolor="white"
    )

    fig.subplots_adjust(
        left=.077,
        right=.985,
        top=.80,
        bottom=.19,
        wspace=.20
    )

    for panel, inc in enumerate(INCOMES):

        ax = axes[panel]

        sub = sub_species[
            sub_species["Income_Group"] == inc
        ]

        for source in SOURCES:

            group = sub[
                sub["Source"] == source
            ]

            values = group[
                "AMR"
            ].to_numpy(dtype=float)

            x = SOURCE_POSITIONS[source]

            if len(values) == 0:
                ax.text(
                    x,
                    global_max * .5,
                    "No data",
                    rotation=90,
                    ha="center",
                    va="center",
                    color="#999999",
                    fontsize=9
                )
                continue

            bp = ax.boxplot(
                [values],
                positions=[x],
                widths=.54,
                patch_artist=True,
                showfliers=False,
                medianprops={
                    "color": "#111111",
                    "linewidth": 1.5
                },
                boxprops={
                    "edgecolor": "#444444",
                    "linewidth": 1.0
                },
                whiskerprops={
                    "color": "#444444"
                },
                capprops={
                    "color": "#444444"
                }
            )

            bp["boxes"][0].set_facecolor(
                SOURCE_COLORS[source]
            )
            bp["boxes"][0].set_alpha(.55)

            if len(values) > MAX_SCATTER_POINTS:
                selected = RNG.choice(
                    values,
                    MAX_SCATTER_POINTS,
                    replace=False
                )
            else:
                selected = values

            jitter = RNG.uniform(
                -.19, .19,
                size=len(selected)
            )

            ax.scatter(
                np.full(len(selected), x) + jitter,
                selected,
                s=6,
                color=SOURCE_COLORS[source],
                alpha=.25,
                edgecolor="none",
                zorder=3
            )

            ax.text(
                x,
                global_max + 1.0,
                f"n={len(values):,}\n"
                f"C={group['Country'].nunique()}",
                ha="center",
                va="bottom",
                fontsize=9
            )

        panel_tests = tests[
            tests["Species"].eq(species)
            & tests["Income_Group"].eq(inc)
        ]

        for source_a, source_b, height in brackets:

            match = panel_tests[
                panel_tests["Source_A"].eq(source_a)
                & panel_tests["Source_B"].eq(source_b)
            ]

            if len(match) != 1:
                raise RuntimeError(
                    "Missing pairwise comparison"
                )

            result = match.iloc[0]
            label = result["Significance"]

            # Do not draw brackets for empty groups.
            if result["Status"] == "EMPTY_GROUP":
                continue

            x1 = SOURCE_POSITIONS[source_a]
            x2 = SOURCE_POSITIONS[source_b]

            ax.plot(
                [x1, x1, x2, x2],
                [
                    height - .55,
                    height,
                    height,
                    height - .55
                ],
                color="#333333",
                linewidth=.95
            )

            ax.text(
                (x1 + x2) / 2,
                height + .35,
                label,
                ha="center",
                va="bottom",
                fontsize=11,
                fontweight=(
                    "bold" if "*" in label
                    else "normal"
                ),
                color=(
                    "#B22222" if "*" in label
                    else "#444444"
                )
            )

        ax.set_xlim(.45, 3.55)
        ax.set_ylim(0, ymax)

        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(SOURCES, fontsize=11)

        ax.set_title(
            INCOME_NAMES[inc],
            fontsize=14,
            fontweight="bold",
            color=INCOME_COLORS[inc],
            pad=12
        )

        ax.grid(
            axis="y",
            color="#E6E6E6",
            linewidth=.65
        )

        ax.set_axisbelow(True)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if panel == 0:
            ax.set_ylabel(
                "Unique AMR genes per genome",
                fontsize=12
            )

    fig.suptitle(
        short_species(species),
        fontsize=22,
        fontweight="bold",
        fontstyle="italic",
        y=.970
    )

    fig.text(
        .5,
        .917,
        "Human–Animal–Environment AMR burden "
        "within country-income categories",
        ha="center",
        fontsize=12
    )

    legend = [
        Patch(
            facecolor=SOURCE_COLORS[s],
            edgecolor="#444444",
            label=s,
            alpha=.65
        )
        for s in SOURCES
    ]

    fig.legend(
        handles=legend,
        loc="lower center",
        bbox_to_anchor=(.5, .115),
        ncol=3,
        frameon=False
    )

    fig.text(
        .5,
        .065,
        "Two-sided genome-level Mann–Whitney U; "
        "BH-FDR across all nonempty comparisons. "
        "ns = not significant; empty groups not tested.",
        ha="center",
        fontsize=9.5
    )

    fig.text(
        .5,
        .042,
        "Unadjusted for country/BioProject clustering; "
        "interpret P-values as exploratory.",
        ha="center",
        fontsize=9.2,
        color="#666666"
    )

    output = OUT / (
        f"Fig03G_{species_tag(species)}_"
        f"ALL_GENOMES_SIGNIFICANCE.jpg"
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

    print(
        "[SAVED]",
        output.name,
        flush=True
    )

# ============================================================
# FINAL RESULTS
# ============================================================

print("\n" + "=" * 110)
print("FIGURE 3G - FINAL RESULTS")
print("=" * 110)

print("Genomes used:", len(data))
print("Potential comparisons:", len(tests))
print("Tests performed:", len(tested))
print(
    "FDR significant:",
    int((tested["FDR"] < .05).sum())
)

print("\nSIGNIFICANCE COUNTS")
print(tests["Significance"].value_counts().to_string())

print("\nALL PERFORMED COMPARISONS")
print(
    tested[[
        "Species",
        "Income_Group",
        "Source_A",
        "Source_B",
        "N_A",
        "N_B",
        "Mean_Difference_B_minus_A",
        "Median_Difference_B_minus_A",
        "Cliffs_Delta",
        "Raw_P",
        "FDR",
        "Significance"
    ]].round(5).to_string(index=False)
)

print("\nOUTPUT EXCEL:", EXCEL)
print("OUTPUT FOLDER:", OUT)
print("JPEG figures: 7")
print("Resolution: 600 DPI")
print("=" * 110)
