
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.special import gammaln
from scipy.stats import kruskal
from matplotlib import font_manager
from matplotlib.lines import Line2D

# ============================================================
# ECC FIGURE 2C
# COUNTRY-LEVEL SPECIES DIVERSITY AND RAREFACTION
#
# Input: Audited country-level species counts
# Dataset: 9157 genomes, 97 countries
#
# Primary rarefaction depth: 20 resolved genomes
# Sensitivity depth: 10 resolved genomes
#
# Exactly:
# - ONE final JPEG, 600 DPI
# - ONE final Excel workbook
#
# Existing audit workbook is preserved.
# Final genomic master is unchanged.
# ============================================================

ROOT = Path(
    r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES"
)

OUT = ROOT / "Fig_02C_Country_Species_Diversity"

AUDIT = OUT / "Fig_02C_Country_Diversity_Audit.xlsx"

JPG = OUT / "Fig_02C_Country_Species_Diversity.jpg"
XLSX = OUT / "Fig_02C_Country_Diversity_Data.xlsx"

GROUPS = ["H", "M", "L+LM"]

NAMES = {
    "H": "High income",
    "M": "Upper-middle",
    "L+LM": "Low + lower-middle"
}

COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82"
}

PRIMARY_DEPTH = 20
SENSITIVITY_DEPTH = 10

UNRESOLVED = (
    "Enterobacter cloacae complex - unresolved"
)

RNG = np.random.default_rng(2026)

# ============================================================
# FONT
# ============================================================

fonts = {
    f.name for f in font_manager.fontManager.ttflist
}

FONT = next(
    (
        f for f in [
            "Arial", "Helvetica",
            "Liberation Sans", "DejaVu Sans"
        ] if f in fonts
    ),
    "DejaVu Sans"
)

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "axes.linewidth": 0.8,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})

# ============================================================
# LOAD AUDIT WORKBOOK
# ============================================================

print("=" * 95)
print("FIGURE 2C: COUNTRY-LEVEL DIVERSITY")
print("=" * 95)

if not AUDIT.exists():
    raise FileNotFoundError(AUDIT)

metrics = pd.read_excel(
    AUDIT,
    sheet_name="COUNTRY_DIVERSITY"
)

matrix = pd.read_excel(
    AUDIT,
    sheet_name="COUNTRY_SPECIES_MATRIX"
)

required_metrics = [
    "Country",
    "Income_Group",
    "Total_Genomes",
    "Resolved_Genomes",
    "Resolved_Species_Richness",
    "Shannon_H"
]

for col in required_metrics:
    if col not in metrics.columns:
        raise ValueError("Missing audit column: " + col)

if metrics["Country"].duplicated().any():
    raise ValueError("Duplicate country in audit")

if matrix["Country"].duplicated().any():
    raise ValueError("Duplicate country in species matrix")

if len(metrics) != 97:
    raise ValueError(
        f"Expected 97 countries, found {len(metrics)}"
    )

if int(metrics["Total_Genomes"].sum()) != 9157:
    raise ValueError("Final genome count mismatch")

species = [
    col for col in matrix.columns
    if col not in [
        "Country",
        "Income_Group",
        "Total_Genomes",
        "Resolved_Genomes",
        UNRESOLVED
    ]
]

if len(species) != 7:
    raise ValueError(
        f"Expected 7 resolved species, found {len(species)}"
    )

matrix = matrix.set_index("Country")

# Reconcile matrix with audit.
for row in metrics.itertuples():
    if row.Country not in matrix.index:
        raise ValueError(
            "Country absent from species matrix: "
            + row.Country
        )

    n = int(matrix.loc[row.Country, species].sum())

    if n != int(row.Resolved_Genomes):
        raise ValueError(
            "Resolved count mismatch: " + row.Country
        )

# ============================================================
# EXACT EXPECTED RAREFIED SPECIES RICHNESS
#
# E[S_m] = sum_i {
#     1 - C(N-n_i,m) / C(N,m)
# }
#
# This is sampling without replacement.
# ============================================================

def log_choose(n, k):
    if k < 0 or k > n:
        return -np.inf

    return (
        gammaln(n + 1)
        - gammaln(k + 1)
        - gammaln(n - k + 1)
    )


def rarefied_richness(counts, depth):

    counts = np.asarray(counts, dtype=int)

    N = int(counts.sum())

    if N < depth or N == 0:
        return np.nan

    denominator = log_choose(N, depth)

    richness = 0.0

    for n_i in counts:

        if n_i <= 0:
            continue

        if N - n_i < depth:
            absent_probability = 0.0
        else:
            absent_probability = np.exp(
                log_choose(N - n_i, depth)
                - denominator
            )

        richness += 1.0 - absent_probability

    return float(richness)


# ============================================================
# CALCULATE RAREFACTION
# ============================================================

print("\nCalculating exact expected richness...")

records = []

for country in metrics["Country"]:

    row = metrics.loc[
        metrics["Country"] == country
    ].iloc[0]

    counts = matrix.loc[
        country, species
    ].to_numpy(dtype=int)

    result = {
        "Country": country,
        "Income_Group": row["Income_Group"],
        "Total_Genomes": int(row["Total_Genomes"]),
        "Resolved_Genomes": int(row["Resolved_Genomes"]),
        "Observed_Richness": int(
            row["Resolved_Species_Richness"]
        ),
        "Shannon_H": row["Shannon_H"],
        "Rarefied_Richness_10": rarefied_richness(
            counts, SENSITIVITY_DEPTH
        ),
        "Rarefied_Richness_20": rarefied_richness(
            counts, PRIMARY_DEPTH
        )
    }

    records.append(result)

rarefaction = pd.DataFrame(records)

primary = rarefaction[
    rarefaction["Income_Group"].isin(GROUPS)
    & rarefaction["Rarefied_Richness_20"].notna()
].copy()

sensitivity = rarefaction[
    rarefaction["Income_Group"].isin(GROUPS)
    & rarefaction["Rarefied_Richness_10"].notna()
].copy()

counts20 = (
    primary["Income_Group"]
    .value_counts()
    .reindex(GROUPS, fill_value=0)
)

counts10 = (
    sensitivity["Income_Group"]
    .value_counts()
    .reindex(GROUPS, fill_value=0)
)

if counts20.to_dict() != {
    "H": 24, "M": 9, "L+LM": 9
}:
    raise ValueError(
        "Depth-20 country counts differ from audit: "
        + str(counts20.to_dict())
    )

print("\nCountries eligible at depth 20:")
print(counts20.to_string())

print("\nCountries eligible at depth 10:")
print(counts10.to_string())

# ============================================================
# SAMPLING-DEPTH ELIGIBILITY
# ============================================================

depths = [1, 5, 10, 20, 30, 50, 100, 200]

eligibility_rows = []

for depth in depths:
    for group in GROUPS:
        sub = rarefaction[
            rarefaction["Income_Group"] == group
        ]

        eligible = int(
            (sub["Resolved_Genomes"] >= depth).sum()
        )

        eligibility_rows.append({
            "Depth": depth,
            "Income_Group": group,
            "Eligible_Countries": eligible
        })

eligibility = pd.DataFrame(eligibility_rows)

# ============================================================
# DESCRIPTIVE SUMMARY
# ============================================================

summary_rows = []

for group in GROUPS:

    sub = primary[
        primary["Income_Group"] == group
    ]

    summary_rows.append({
        "Income_Group": group,
        "Eligible_Countries": len(sub),
        "Median_Rarefied_Richness_20": float(
            sub["Rarefied_Richness_20"].median()
        ),
        "Mean_Rarefied_Richness_20": float(
            sub["Rarefied_Richness_20"].mean()
        ),
        "Median_Shannon_H": float(
            sub["Shannon_H"].median()
        ),
        "Median_Resolved_Genomes": float(
            sub["Resolved_Genomes"].median()
        )
    })

summary = pd.DataFrame(summary_rows)

# Exploratory country-level omnibus comparisons.
# No individual-genome independence assumption.
# Interpret cautiously given country heterogeneity.
stat_rows = []

for metric in [
    "Rarefied_Richness_20",
    "Shannon_H"
]:
    arrays = [
        primary.loc[
            primary["Income_Group"] == g, metric
        ].dropna().to_numpy()
        for g in GROUPS
    ]

    if all(len(x) > 0 for x in arrays):
        if len(np.unique(np.concatenate(arrays))) > 1:
            statistic, p = kruskal(*arrays)
        else:
            statistic, p = 0.0, 1.0

        stat_rows.append({
            "Metric": metric,
            "Test": "Kruskal-Wallis, country-level",
            "Statistic": statistic,
            "P_Value": p,
            "Countries": sum(len(x) for x in arrays)
        })

stats = pd.DataFrame(stat_rows)

# ============================================================
# FIGURE DESIGN
# ============================================================

print("\nRendering Figure 2C...")

fig, axes = plt.subplots(
    1, 3,
    figsize=(18.2, 7.7),
    gridspec_kw={
        "width_ratios": [1.25, 1.0, 1.0],
        "wspace": 0.38
    }
)

fig.suptitle(
    "Country-level diversity of the "
    "Enterobacter cloacae complex",
    fontsize=18,
    weight="bold",
    y=0.975
)

fig.text(
    0.5, 0.930,
    "Genome sampling coverage, standardized "
    "species richness and Shannon diversity",
    ha="center",
    fontsize=11.5,
    color="#555555"
)

# ============================================================
# PANEL 1: ELIGIBLE COUNTRIES BY DEPTH
# ============================================================

ax = axes[0]

for group in GROUPS:

    sub = eligibility[
        eligibility["Income_Group"] == group
    ]

    ax.plot(
        sub["Depth"],
        sub["Eligible_Countries"],
        marker="o",
        markersize=5,
        linewidth=2.0,
        color=COLORS[group],
        label=NAMES[group]
    )

ax.axvline(
    PRIMARY_DEPTH,
    color="#555555",
    linestyle="--",
    linewidth=1
)

ax.set_xscale("log")

ax.set_xticks(depths)
ax.set_xticklabels(
    [str(x) for x in depths],
    rotation=35
)

ax.set_xlabel("Minimum resolved genomes per country")
ax.set_ylabel("Eligible countries")

ax.set_title(
    "Sampling-depth coverage",
    loc="left",
    weight="bold",
    pad=12
)

ax.legend(
    frameon=False,
    fontsize=9,
    loc="upper right"
)

# ============================================================
# PANELS 2 AND 3: COUNTRY-LEVEL DISTRIBUTIONS
# ============================================================

def draw_country_metric(ax, metric, title, ylabel):

    positions = np.arange(len(GROUPS))

    for idx, group in enumerate(GROUPS):

        values = primary.loc[
            primary["Income_Group"] == group,
            metric
        ].dropna().to_numpy(dtype=float)

        # Light box with sample points to show
        # each country rather than hide variation.
        ax.boxplot(
            [values],
            positions=[idx],
            widths=0.48,
            showfliers=False,
            patch_artist=True,
            medianprops={
                "color": "#222222",
                "linewidth": 1.7
            },
            boxprops={
                "facecolor": COLORS[group],
                "alpha": 0.17,
                "edgecolor": COLORS[group],
                "linewidth": 1.2
            },
            whiskerprops={
                "color": COLORS[group],
                "linewidth": 1.1
            },
            capprops={
                "color": COLORS[group],
                "linewidth": 1.1
            }
        )

        jitter = RNG.uniform(
            -0.14, 0.14,
            len(values)
        )

        ax.scatter(
            np.full(len(values), idx) + jitter,
            values,
            color=COLORS[group],
            s=34,
            alpha=0.84,
            edgecolor="white",
            linewidth=0.45,
            zorder=3
        )

        ax.text(
            idx, 1.035,
            f"n={len(values)}",
            ha="center",
            transform=ax.get_xaxis_transform(),
            fontsize=9.5,
            color="#555555"
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(
        ["High", "Upper-middle", "Low +\nlower-middle"],
        fontsize=10
    )

    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", weight="bold", pad=12)
    ax.set_xlim(-0.6, 2.6)


draw_country_metric(
    axes[1],
    "Rarefied_Richness_20",
    "Species richness after rarefaction",
    "Expected resolved species at n = 20"
)

draw_country_metric(
    axes[2],
    "Shannon_H",
    "Observed Shannon diversity",
    "Shannon diversity (H)"
)

axes[1].set_ylim(0, 7.25)
axes[2].set_ylim(bottom=0)

for ax in axes:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#666666")
    ax.spines["bottom"].set_color("#666666")

    ax.grid(
        axis="y",
        color="#E6E6E6",
        linewidth=0.55,
        alpha=0.8
    )

    ax.set_axisbelow(True)

fig.subplots_adjust(
    left=0.075,
    right=0.975,
    top=0.83,
    bottom=0.185
)

fig.text(
    0.5, 0.060,
    "Rarefaction uses 20 resolved genomes per country "
    "(24 high-income, 9 upper-middle, 9 low/lower-middle). "
    "Unresolved ECC annotations are excluded from "
    "species-diversity calculations.",
    ha="center",
    fontsize=10,
    color="#555555"
)

fig.savefig(
    JPG,
    dpi=600,
    format="jpeg",
    facecolor="white",
    pil_kwargs={
        "quality": 95,
        "subsampling": 0,
        "optimize": True
    }
)

plt.close(fig)

# ============================================================
# EXCEL SOURCE DATA
# ============================================================

methods = pd.DataFrame([
    (
        "Master",
        "Final 9157 ECC genomes, 2000-2026"
    ),
    (
        "Species categories",
        "Seven resolved species; unresolved ECC excluded"
    ),
    (
        "Primary rarefaction depth",
        PRIMARY_DEPTH
    ),
    (
        "Sensitivity rarefaction depth",
        SENSITIVITY_DEPTH
    ),
    (
        "Rarefaction method",
        "Exact expected richness under sampling "
        "without replacement"
    ),
    (
        "Country comparison unit",
        "Country/territory"
    ),
    (
        "Country eligibility",
        "At least 20 resolved genomes"
    ),
    (
        "Shannon",
        "Observed Shannon H; not rarefied or "
        "bias-corrected"
    ),
    (
        "Statistical comparison",
        "Exploratory country-level Kruskal-Wallis"
    ),
    (
        "Interpretation",
        "Sequencing effort, sampling selection and "
        "country composition remain potential biases"
    ),
    (
        "Created",
        datetime.now().isoformat()
    )
], columns=["Parameter", "Value"])

with pd.ExcelWriter(
    XLSX,
    engine="openpyxl"
) as writer:

    rarefaction.to_excel(
        writer,
        sheet_name="ALL_COUNTRIES",
        index=False
    )

    primary.to_excel(
        writer,
        sheet_name="PRIMARY_DEPTH20",
        index=False
    )

    sensitivity.to_excel(
        writer,
        sheet_name="SENSITIVITY_DEPTH10",
        index=False
    )

    eligibility.to_excel(
        writer,
        sheet_name="DEPTH_ELIGIBILITY",
        index=False
    )

    summary.to_excel(
        writer,
        sheet_name="INCOME_SUMMARY",
        index=False
    )

    stats.to_excel(
        writer,
        sheet_name="COUNTRY_LEVEL_TESTS",
        index=False
    )

    methods.to_excel(
        writer,
        sheet_name="METHODS",
        index=False
    )

# ============================================================
# FINAL REPORT
# ============================================================

if not JPG.exists() or not XLSX.exists():
    raise RuntimeError("Output files missing")

print("\n" + "=" * 95)
print("FIGURE 2C COMPLETED")
print("=" * 95)

print("\nCOUNTRY-LEVEL DIVERSITY SUMMARY")
print(summary.to_string(index=False))

print("\nEXPLORATORY COUNTRY-LEVEL TESTS")
print(stats.to_string(index=False))

print("\nOUTPUTS")
print(JPG)
print(XLSX)

print("\nFigure: 600 DPI JPEG")
print("Excel: 7 worksheets")
print("Original audit and final master unchanged")
print("=" * 95)
