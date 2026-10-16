
from pathlib import Path
from itertools import combinations, product
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
# FIGURE 3G - ONE HEALTH AMR BURDEN WITH SIGNIFICANCE
#
# Species: all 7 resolved ECC species
# Income: H, M, L+LM
# Sources: Human, Animal, Environment
#
# Genome-level descriptive boxplots.
#
# COUNTRY-LEVEL INFERENCE:
# 1. Paired country comparison where >=5 countries
#    have both sources with >=5 genomes/source.
# 2. Otherwise, unpaired Mann-Whitney U using only
#    countries exclusive to each source, provided
#    each group has >=5 eligible countries.
# 3. Otherwise NT (not tested).
#
# This prevents countries contributing to both
# samples of an unpaired test.
#
# BH FDR across all performed comparisons.
# Plot all results: stars, ns, NT.
#
# Seven JPEGs at 600 DPI + full Excel workbook.
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

EXCEL = OUT / "Fig03G_AMR_Burden_WITH_STATISTICS.xlsx"

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

MIN_GENOMES_COUNTRY = 5
MIN_COUNTRIES_GROUP = 5

N_BOOT = 5000
N_SIGNFLIP = 9999
MAX_DISPLAY_POINTS = 300

RNG = np.random.default_rng(20261010)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    q = np.full(len(p), np.nan)

    valid = np.flatnonzero(np.isfinite(p))
    if not len(valid):
        return q

    order = valid[np.argsort(p[valid])]
    n = len(order)

    values = p[order] * n / np.arange(1, n + 1)
    values = np.minimum.accumulate(values[::-1])[::-1]

    q[order] = np.clip(values, 0, 1)
    return q


def significance(q):
    if not np.isfinite(q):
        return "NT"
    if q < 0.0001:
        return "****"
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def paired_test(a, b):
    d = np.asarray(b - a, dtype=float)
    effect = float(np.mean(d))

    if np.all(d == 0):
        return effect, 1.0

    if len(d) <= 16:
        signs = np.asarray(
            list(product([-1, 1], repeat=len(d))),
            dtype=float
        )
        null = (signs * d).mean(axis=1)
        p = np.mean(
            np.abs(null) >= abs(effect) - 1e-12
        )
    else:
        extreme = 0
        for _ in range(N_SIGNFLIP):
            signs = RNG.choice([-1, 1], size=len(d))
            null = float(np.mean(signs * d))
            extreme += (
                abs(null) >= abs(effect) - 1e-12
            )
        p = (extreme + 1) / (N_SIGNFLIP + 1)

    return effect, float(p)


def bootstrap_ci(a, b, paired):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    if paired:
        d = b - a
        idx = RNG.integers(
            0, len(d), size=(N_BOOT, len(d))
        )
        estimates = d[idx].mean(axis=1)
    else:
        ia = RNG.integers(
            0, len(a), size=(N_BOOT, len(a))
        )
        ib = RNG.integers(
            0, len(b), size=(N_BOOT, len(b))
        )
        estimates = (
            b[ib].mean(axis=1)
            - a[ia].mean(axis=1)
        )

    return (
        float(np.quantile(estimates, .025)),
        float(np.quantile(estimates, .975))
    )


# ============================================================
# LOAD MASTER AND INCOME DATA
# ============================================================

print("=" * 100)
print("FIGURE 3G - AMR BURDEN WITH STATISTICAL SIGNIFICANCE")
print("=" * 100)

df = pd.read_excel(
    MASTER,
    sheet_name="FINAL_MASTER_2000_2026",
    engine="openpyxl",
    dtype={"Assembly_Accession": str}
)

if len(df) != 9157:
    raise ValueError("Expected 9157 master genomes")

required = [
    "Assembly_Accession",
    "Normalized_Species",
    "Final_Source_Category",
    "Country_Harmonized",
    "Unique_AMR_Gene_Count",
    "AMRFinder_Status"
]

for field in required:
    if field not in df:
        raise ValueError("Missing field: " + field)

if (
    df["Assembly_Accession"].isna().any()
    or df["Assembly_Accession"].duplicated().any()
):
    raise ValueError("Invalid assembly accessions")

if not (
    df["AMRFinder_Status"]
    .astype(str)
    .str.strip()
    .str.upper()
    == "COMPLETE"
).all():
    raise ValueError("Incomplete AMRFinder data")

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

if df["AMR"].isna().any() or (df["AMR"] < 0).any():
    raise ValueError("Invalid AMR burden values")

income = pd.read_excel(
    INCOME_FILE,
    sheet_name="COUNTRY_INCOME_COUNTS",
    engine="openpyxl"
)

for field in [
    "Country_Harmonized",
    "Income_Group",
    "Genome_Count"
]:
    if field not in income:
        raise ValueError("Missing income field: " + field)

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
    raise ValueError("Duplicate country income assignments")

reference = income.set_index("Country")["Genome_Count"]
observed = df["Country"].value_counts()

check = pd.concat(
    [observed.rename("Observed"),
     reference.rename("Reference")],
    axis=1
)

if (
    check.isna().any().any()
    or not (
        check["Observed"].astype(int)
        == check["Reference"].astype(int)
    ).all()
):
    raise ValueError("Income reference count mismatch")

df = df.merge(
    income[["Country", "Income_Group"]],
    on="Country",
    how="left",
    validate="many_to_one"
)

data = df[
    df["Species"].isin(SPECIES)
    & df["Source"].isin(SOURCES)
    & df["Income_Group"].isin(INCOMES)
].copy()

print("Genomes analyzed:", len(data))

# ============================================================
# COUNTRY-LEVEL SUMMARY
# ============================================================

country = (
    data.groupby([
        "Species", "Income_Group",
        "Source", "Country"
    ])
    .agg(
        Genome_N=("AMR", "size"),
        Mean_AMR=("AMR", "mean"),
        Median_AMR=("AMR", "median")
    )
    .reset_index()
)

country["Eligible"] = (
    country["Genome_N"] >= MIN_GENOMES_COUNTRY
)

eligible = country[country["Eligible"]].copy()

group_summary = (
    data.groupby(["Species", "Income_Group", "Source"])
    .agg(
        Genomes=("AMR", "size"),
        Countries=("Country", "nunique"),
        Mean_AMR=("AMR", "mean"),
        Median_AMR=("AMR", "median"),
        SD_AMR=("AMR", "std"),
        Q1=("AMR", lambda x: x.quantile(.25)),
        Q3=("AMR", lambda x: x.quantile(.75))
    )
    .reset_index()
)

# ============================================================
# ALL 63 SOURCE PAIRS
# ============================================================

print("\nRunning country-aware pairwise comparisons...")

stats = []

for species in SPECIES:
    for inc in INCOMES:
        for source_a, source_b in PAIRS:

            a = eligible[
                eligible["Species"].eq(species)
                & eligible["Income_Group"].eq(inc)
                & eligible["Source"].eq(source_a)
            ].set_index("Country")["Mean_AMR"]

            b = eligible[
                eligible["Species"].eq(species)
                & eligible["Income_Group"].eq(inc)
                & eligible["Source"].eq(source_b)
            ].set_index("Country")["Mean_AMR"]

            common = a.index.intersection(b.index)
            only_a = a.index.difference(b.index)
            only_b = b.index.difference(a.index)

            method = "NOT_TESTED"
            used_a = np.array([], dtype=float)
            used_b = np.array([], dtype=float)
            effect = np.nan
            p = np.nan
            ci_low = np.nan
            ci_high = np.nan

            if len(common) >= MIN_COUNTRIES_GROUP:
                method = "PAIRED_SIGN_FLIP"
                used_a = a.loc[common].to_numpy()
                used_b = b.loc[common].to_numpy()
                effect, p = paired_test(used_a, used_b)
                ci_low, ci_high = bootstrap_ci(
                    used_a, used_b, paired=True
                )

            elif (
                len(only_a) >= MIN_COUNTRIES_GROUP
                and len(only_b) >= MIN_COUNTRIES_GROUP
            ):
                method = "UNPAIRED_MANN_WHITNEY"
                used_a = a.loc[only_a].to_numpy()
                used_b = b.loc[only_b].to_numpy()

                effect = float(
                    used_b.mean() - used_a.mean()
                )

                p = float(
                    mannwhitneyu(
                        used_a,
                        used_b,
                        alternative="two-sided",
                        method="auto"
                    ).pvalue
                )

                ci_low, ci_high = bootstrap_ci(
                    used_a, used_b, paired=False
                )

            stats.append({
                "Species": species,
                "Income_Group": inc,
                "Source_A": source_a,
                "Source_B": source_b,
                "Comparison": (
                    source_b + " minus " + source_a
                ),
                "Method": method,
                "Common_Eligible_Countries": len(common),
                "Exclusive_Countries_A": len(only_a),
                "Exclusive_Countries_B": len(only_b),
                "Countries_Used_A": len(used_a),
                "Countries_Used_B": len(used_b),
                "Mean_A_Used": (
                    float(used_a.mean())
                    if len(used_a) else np.nan
                ),
                "Mean_B_Used": (
                    float(used_b.mean())
                    if len(used_b) else np.nan
                ),
                "Difference_B_minus_A": effect,
                "CI95_Lower": ci_low,
                "CI95_Upper": ci_high,
                "Raw_P": p
            })

tests = pd.DataFrame(stats)

if len(tests) != 63:
    raise RuntimeError("Expected 63 comparisons")

tests["FDR"] = bh_adjust(tests["Raw_P"].to_numpy())
tests["Significance"] = tests["FDR"].apply(significance)

print("Possible comparisons:", len(tests))
print("Performed:", int(tests["Raw_P"].notna().sum()))
print("FDR significant:", int((tests["FDR"] < .05).sum()))

# ============================================================
# SAVE STATISTICAL MASTER BEFORE PLOTTING
# ============================================================

all_groups = pd.MultiIndex.from_product(
    [SPECIES, INCOMES, SOURCES],
    names=["Species", "Income_Group", "Source"]
)

group_summary = (
    group_summary.set_index(
        ["Species", "Income_Group", "Source"]
    )
    .reindex(all_groups)
    .reset_index()
)

group_summary["Genomes"] = (
    group_summary["Genomes"].fillna(0).astype(int)
)
group_summary["Countries"] = (
    group_summary["Countries"].fillna(0).astype(int)
)

methods = pd.DataFrame([
    ("Outcome", "Unique annotated AMR genes per genome"),
    ("Descriptive plots", "Genome-level boxplots"),
    ("Inference unit", "Country-level mean AMR burden"),
    ("Eligible country", "At least 5 genomes per source/species/country"),
    ("Paired rule", "At least 5 overlapping eligible countries"),
    ("Fallback rule", "At least 5 source-exclusive eligible countries per group"),
    ("Paired test", "Two-sided sign-flip permutation"),
    ("Unpaired test", "Two-sided Mann-Whitney U"),
    ("Effect", "Difference in country-equal means, B minus A"),
    ("CI", "Country-level bootstrap 95% interval"),
    ("FDR", "Global BH across all performed comparisons"),
    ("NT", "Insufficient eligible country representation"),
    ("Interpretation", "Associations, not causal effects"),
    ("Important limitation", "Different country compositions can confound unpaired comparisons"),
    ("Another limitation", "BioProject, year and lineage not adjusted"),
    ("Generated", datetime.now().isoformat())
], columns=["Parameter", "Description"])

EXCEL.parent.mkdir(parents=True, exist_ok=True)

with pd.ExcelWriter(EXCEL, engine="openpyxl") as writer:
    tests.to_excel(
        writer,
        sheet_name="ALL_63_PAIRWISE_FDR",
        index=False
    )
    tests[tests["Raw_P"].notna()].to_excel(
        writer,
        sheet_name="TESTED_ONLY",
        index=False
    )
    group_summary.to_excel(
        writer,
        sheet_name="GENOME_GROUP_SUMMARY",
        index=False
    )
    country.to_excel(
        writer,
        sheet_name="COUNTRY_LEVEL_DATA",
        index=False
    )
    data[[
        "Assembly_Accession",
        "Species",
        "Source",
        "Income_Group",
        "Country",
        "AMR"
    ]].to_excel(
        writer,
        sheet_name="GENOME_LEVEL_DATA",
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

print("\nGenerating seven 600-DPI significance figures...")

source_positions = {
    "Human": 1,
    "Animal": 2,
    "Environment": 3
}

income_names = {
    "H": "High income",
    "M": "Upper-middle income",
    "L+LM": "Low + lower-middle income"
}

overall_max = float(data["AMR"].max())

# Shared scale across all species and panels.
plot_top = overall_max + 15

for species in SPECIES:

    sub_species = data[data["Species"] == species]

    fig, axes = plt.subplots(
        1, 3,
        figsize=(18, 8),
        sharey=True,
        facecolor="white"
    )

    fig.subplots_adjust(
        left=.075,
        right=.985,
        top=.80,
        bottom=.20,
        wspace=.18
    )

    for panel, inc in enumerate(INCOMES):

        ax = axes[panel]
        sub = sub_species[
            sub_species["Income_Group"] == inc
        ]

        for source in SOURCES:

            group = sub[sub["Source"] == source]
            values = group["AMR"].to_numpy(dtype=float)
            x = source_positions[source]

            if not len(values):
                ax.text(
                    x, overall_max * .5,
                    "No data",
                    ha="center",
                    rotation=90,
                    color="#999999"
                )
                continue

            bp = ax.boxplot(
                [values],
                positions=[x],
                widths=.5,
                patch_artist=True,
                showfliers=False,
                medianprops={
                    "color": "#111111",
                    "linewidth": 1.6
                },
                boxprops={
                    "edgecolor": "#444444"
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

            selected = values

            if len(selected) > 300:
                selected = RNG.choice(
                    selected,
                    size=300,
                    replace=False
                )

            jitter = RNG.uniform(
                -.18, .18, size=len(selected)
            )

            ax.scatter(
                np.full(len(selected), x) + jitter,
                selected,
                s=6,
                color=SOURCE_COLORS[source],
                alpha=.26,
                edgecolor="none",
                zorder=3
            )

            ax.text(
                x,
                overall_max + 2.2,
                f"n={len(values):,}\n"
                f"C={group['Country'].nunique()}",
                ha="center",
                va="bottom",
                fontsize=9
            )

        # ----------------------------------------------------
        # THREE SIGNIFICANCE BRACKETS
        # ----------------------------------------------------

        panel_tests = tests[
            tests["Species"].eq(species)
            & tests["Income_Group"].eq(inc)
        ]

        bracket_specs = [
            ("Human", "Animal", overall_max + 6),
            ("Animal", "Environment", overall_max + 9),
            ("Human", "Environment", overall_max + 12)
        ]

        for source_a, source_b, height in bracket_specs:

            result = panel_tests[
                panel_tests["Source_A"].eq(source_a)
                & panel_tests["Source_B"].eq(source_b)
            ]

            if len(result) != 1:
                raise RuntimeError("Missing statistical result")

            result = result.iloc[0]
            label = result["Significance"]

            x1 = source_positions[source_a]
            x2 = source_positions[source_b]

            ax.plot(
                [x1, x1, x2, x2],
                [
                    height - .5,
                    height,
                    height,
                    height - .5
                ],
                color="#444444",
                linewidth=.9
            )

            ax.text(
                (x1 + x2) / 2,
                height + .25,
                label,
                ha="center",
                va="bottom",
                fontsize=11,
                fontweight=(
                    "bold" if "*" in label
                    else "normal"
                ),
                color=(
                    "#A71919" if "*" in label
                    else "#444444"
                )
            )

        ax.set_xlim(.45, 3.55)
        ax.set_ylim(0, plot_top)

        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(SOURCES, fontsize=11)

        ax.set_title(
            income_names[inc],
            fontsize=14,
            fontweight="bold",
            color=INCOME_COLORS[inc],
            pad=10
        )

        ax.grid(
            axis="y",
            color="#E5E5E5",
            linewidth=.6
        )

        ax.set_axisbelow(True)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if panel == 0:
            ax.set_ylabel(
                "Unique annotated AMR genes per genome",
                fontsize=12
            )

    fig.suptitle(
        species.replace("Enterobacter ", "E. "),
        fontsize=22,
        fontweight="bold",
        fontstyle="italic",
        y=.97
    )

    fig.text(
        .5,
        .918,
        "One Health AMR burden comparisons "
        "with country-level statistical significance",
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
        .062,
        "Stars and ns reflect BH-adjusted country-level tests; "
        "NT = insufficient eligible countries. "
        "Full method, effect sizes, P and FDR: Excel.",
        ha="center",
        fontsize=9.5
    )

    fig.text(
        .5,
        .041,
        "Boxplots display all genomes; country means form "
        "the inferential observations. "
        "Paired and unpaired country designs are identified in Excel.",
        ha="center",
        fontsize=9,
        color="#666666"
    )

    tag = species.replace(
        "Enterobacter ", "E_"
    ).replace(" ", "_")

    output = OUT / (
        f"Fig03G_{tag}_WITH_SIGNIFICANCE.jpg"
    )

    fig.savefig(
        output,
        dpi=600,
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

    print("[SAVED]", output.name, flush=True)

# ============================================================
# RESULTS
# ============================================================

print("\n" + "=" * 100)
print("FIGURE 3G - COMPLETE STATISTICAL RESULTS")
print("=" * 100)

print("Possible comparisons:", len(tests))
print("Performed:", tests["Raw_P"].notna().sum())
print("FDR significant:", (tests["FDR"] < .05).sum())

print("\nTEST TYPES")
print(tests["Method"].value_counts().to_string())

print("\nALL PERFORMED COMPARISONS")
performed = tests[tests["Raw_P"].notna()]

if performed.empty:
    print("None passed country eligibility")
else:
    print(
        performed[[
            "Species", "Income_Group",
            "Source_A", "Source_B",
            "Method",
            "Countries_Used_A",
            "Countries_Used_B",
            "Difference_B_minus_A",
            "CI95_Lower", "CI95_Upper",
            "Raw_P", "FDR",
            "Significance"
        ]].round(5).to_string(index=False)
    )

print("\nExcel:", EXCEL)
print("Seven 600-DPI species figures saved in:", OUT)
print("Original master unchanged")
print("=" * 100)
