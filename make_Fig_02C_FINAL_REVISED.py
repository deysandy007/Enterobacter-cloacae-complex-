from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.special import gammaln
from scipy.stats import kruskal
from matplotlib import font_manager

# ============================================================
# ECC FIGURE 2C - REVISED FINAL
# COUNTRY-LEVEL SPECIES DIVERSITY AND RAREFACTION
#
# FIXES IN THIS VERSION:
# - keep same analysis as previous final version
# - move eligible-country counts (n) to x-axis bottom
# - remove overlapping n labels from top of panels
# - slightly cleaner publication-style layout
#
# OUTPUTS:
# 1 JPEG (600 DPI)
# 1 Excel workbook
# ============================================================

ROOT = Path(r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES")
OUT = ROOT / "Fig_02C_Country_Species_Diversity"

AUDIT = OUT / "Fig_02C_Country_Diversity_Audit.xlsx"

JPG = OUT / "Fig_02C_Country_Species_Diversity.jpg"
XLSX = OUT / "Fig_02C_Country_Diversity_Data.xlsx"

GROUPS = ["H", "M", "L+LM"]

GROUP_NAMES = {
    "H": "High income",
    "M": "Upper-middle",
    "L+LM": "Low + lower-middle"
}

GROUP_SHORT = {
    "H": "High income",
    "M": "Upper-middle",
    "L+LM": "Low +\nlower-middle"
}

GROUP_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82"
}

PRIMARY_DEPTH = 20
SENSITIVITY_DEPTH = 10

UNRESOLVED = "Enterobacter cloacae complex - unresolved"
RNG = np.random.default_rng(2026)

# ============================================================
# FONT SETUP
# ============================================================

available_fonts = {f.name for f in font_manager.fontManager.ttflist}
FONT = next(
    (f for f in ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"] if f in available_fonts),
    "DejaVu Sans"
)

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})

# ============================================================
# LOAD AUDIT WORKBOOK
# ============================================================

print("=" * 95)
print("FIGURE 2C: COUNTRY-LEVEL DIVERSITY (REVISED)")
print("=" * 95)

if not AUDIT.exists():
    raise FileNotFoundError(AUDIT)

metrics = pd.read_excel(AUDIT, sheet_name="COUNTRY_DIVERSITY")
matrix = pd.read_excel(AUDIT, sheet_name="COUNTRY_SPECIES_MATRIX")

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
        raise ValueError(f"Missing audit column: {col}")

if metrics["Country"].duplicated().any():
    raise ValueError("Duplicate country in COUNTRY_DIVERSITY")

if matrix["Country"].duplicated().any():
    raise ValueError("Duplicate country in COUNTRY_SPECIES_MATRIX")

if len(metrics) != 97:
    raise ValueError(f"Expected 97 countries, found {len(metrics)}")

if int(metrics["Total_Genomes"].sum()) != 9157:
    raise ValueError("Total genome count mismatch; expected 9157")

species_cols = [
    c for c in matrix.columns
    if c not in [
        "Country",
        "Income_Group",
        "Total_Genomes",
        "Resolved_Genomes",
        UNRESOLVED
    ]
]

if len(species_cols) != 7:
    raise ValueError(f"Expected 7 resolved species, found {len(species_cols)}")

matrix = matrix.set_index("Country")

# Validate resolved counts
for row in metrics.itertuples():
    if row.Country not in matrix.index:
        raise ValueError(f"{row.Country} absent from species matrix")
    resolved_sum = int(matrix.loc[row.Country, species_cols].sum())
    if resolved_sum != int(row.Resolved_Genomes):
        raise ValueError(
            f"Resolved genome mismatch for {row.Country}: "
            f"{resolved_sum} vs {row.Resolved_Genomes}"
        )

# ============================================================
# EXACT INDIVIDUAL-BASED RAREFACTION
# Expected richness at depth m without replacement
# ============================================================

def log_choose(n, k):
    if k < 0 or k > n:
        return -np.inf
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)

def rarefied_richness(counts, depth):
    counts = np.asarray(counts, dtype=int)
    N = int(counts.sum())

    if N < depth or N == 0:
        return np.nan

    denom = log_choose(N, depth)
    expected_richness = 0.0

    for n_i in counts:
        if n_i <= 0:
            continue

        if (N - n_i) < depth:
            absent_prob = 0.0
        else:
            absent_prob = np.exp(log_choose(N - n_i, depth) - denom)

        expected_richness += 1.0 - absent_prob

    return float(expected_richness)

# ============================================================
# CALCULATE RAREFIED RICHNESS
# ============================================================

print("\nCalculating exact expected richness...")

records = []

for country in metrics["Country"]:
    row = metrics.loc[metrics["Country"] == country].iloc[0]
    counts = matrix.loc[country, species_cols].to_numpy(dtype=int)

    records.append({
        "Country": country,
        "Income_Group": row["Income_Group"],
        "Total_Genomes": int(row["Total_Genomes"]),
        "Resolved_Genomes": int(row["Resolved_Genomes"]),
        "Observed_Richness": int(row["Resolved_Species_Richness"]),
        "Shannon_H": row["Shannon_H"],
        "Rarefied_Richness_10": rarefied_richness(counts, SENSITIVITY_DEPTH),
        "Rarefied_Richness_20": rarefied_richness(counts, PRIMARY_DEPTH)
    })

rarefaction = pd.DataFrame(records)

primary = rarefaction[
    rarefaction["Income_Group"].isin(GROUPS)
    & rarefaction["Rarefied_Richness_20"].notna()
].copy()

sensitivity = rarefaction[
    rarefaction["Income_Group"].isin(GROUPS)
    & rarefaction["Rarefied_Richness_10"].notna()
].copy()

counts20 = primary["Income_Group"].value_counts().reindex(GROUPS, fill_value=0)
counts10 = sensitivity["Income_Group"].value_counts().reindex(GROUPS, fill_value=0)

print("\nCountries eligible at depth 20:")
print(counts20.to_string())

print("\nCountries eligible at depth 10:")
print(counts10.to_string())

expected_depth20 = {"H": 24, "M": 9, "L+LM": 9}
if counts20.to_dict() != expected_depth20:
    raise ValueError(
        "Depth-20 country counts do not match expected audit values: "
        + str(counts20.to_dict())
    )

# ============================================================
# ELIGIBILITY CURVE TABLE
# ============================================================

depths = [1, 5, 10, 20, 30, 50, 100, 200]
eligibility_rows = []

for depth in depths:
    for group in GROUPS:
        sub = rarefaction[rarefaction["Income_Group"] == group]
        eligible = int((sub["Resolved_Genomes"] >= depth).sum())
        eligibility_rows.append({
            "Depth": depth,
            "Income_Group": group,
            "Eligible_Countries": eligible
        })

eligibility = pd.DataFrame(eligibility_rows)

# ============================================================
# SUMMARY TABLES
# ============================================================

summary_rows = []
for group in GROUPS:
    sub = primary[primary["Income_Group"] == group]
    summary_rows.append({
        "Income_Group": group,
        "Eligible_Countries": len(sub),
        "Median_Rarefied_Richness_20": float(sub["Rarefied_Richness_20"].median()),
        "Mean_Rarefied_Richness_20": float(sub["Rarefied_Richness_20"].mean()),
        "Median_Shannon_H": float(sub["Shannon_H"].median()),
        "Median_Resolved_Genomes": float(sub["Resolved_Genomes"].median())
    })

summary = pd.DataFrame(summary_rows)

stats_rows = []
for metric in ["Rarefied_Richness_20", "Shannon_H"]:
    arrays = [
        primary.loc[primary["Income_Group"] == g, metric].dropna().to_numpy()
        for g in GROUPS
    ]

    if all(len(x) > 0 for x in arrays):
        if len(np.unique(np.concatenate(arrays))) > 1:
            statistic, p_value = kruskal(*arrays)
        else:
            statistic, p_value = 0.0, 1.0

        stats_rows.append({
            "Metric": metric,
            "Test": "Kruskal-Wallis, country-level",
            "Statistic": statistic,
            "P_Value": p_value,
            "Countries": int(sum(len(x) for x in arrays))
        })

stats = pd.DataFrame(stats_rows)

# ============================================================
# FIGURE
# ============================================================

print("\nRendering revised Figure 2C...")

fig, axes = plt.subplots(
    1, 3,
    figsize=(18.2, 7.9),
    gridspec_kw={"width_ratios": [1.25, 1.0, 1.0], "wspace": 0.38}
)

fig.suptitle(
    "Country-level diversity of the Enterobacter cloacae complex",
    fontsize=18,
    weight="bold",
    y=0.975
)

fig.text(
    0.5, 0.932,
    "Sampling-depth coverage, standardized species richness and Shannon diversity",
    ha="center",
    fontsize=11.5,
    color="#555555"
)

# ------------------------------------------------------------
# PANEL 1: ELIGIBILITY CURVES
# ------------------------------------------------------------

ax = axes[0]

for group in GROUPS:
    sub = eligibility[eligibility["Income_Group"] == group]
    ax.plot(
        sub["Depth"],
        sub["Eligible_Countries"],
        marker="o",
        markersize=5,
        linewidth=2.0,
        color=GROUP_COLORS[group],
        label=GROUP_NAMES[group]
    )

ax.axvline(PRIMARY_DEPTH, color="#555555", linestyle="--", linewidth=1.0)

ax.set_xscale("log")
ax.set_xticks(depths)
ax.set_xticklabels([str(x) for x in depths], rotation=35)
ax.set_xlabel("Minimum resolved genomes per country")
ax.set_ylabel("Eligible countries")
ax.set_title("Sampling-depth coverage", loc="left", weight="bold", pad=12)

ax.legend(frameon=False, fontsize=9, loc="upper right")

# ------------------------------------------------------------
# HELPER FOR PANELS 2 AND 3
# ------------------------------------------------------------

def draw_box_scatter(ax, metric, title, ylabel):
    positions = np.arange(len(GROUPS))
    values_by_group = []

    for idx, group in enumerate(GROUPS):
        values = primary.loc[
            primary["Income_Group"] == group, metric
        ].dropna().to_numpy(dtype=float)
        values_by_group.append(values)

        ax.boxplot(
            [values],
            positions=[idx],
            widths=0.48,
            showfliers=False,
            patch_artist=True,
            medianprops={"color": "#222222", "linewidth": 1.7},
            boxprops={
                "facecolor": GROUP_COLORS[group],
                "alpha": 0.18,
                "edgecolor": GROUP_COLORS[group],
                "linewidth": 1.2
            },
            whiskerprops={"color": GROUP_COLORS[group], "linewidth": 1.1},
            capprops={"color": GROUP_COLORS[group], "linewidth": 1.1}
        )

        jitter = RNG.uniform(-0.14, 0.14, len(values))
        ax.scatter(
            np.full(len(values), idx) + jitter,
            values,
            color=GROUP_COLORS[group],
            s=34,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.45,
            zorder=3
        )

    n_labels = [
        f"{GROUP_SHORT[g]}\n(n={len(values_by_group[i])})"
        for i, g in enumerate(GROUPS)
    ]

    ax.set_xticks(positions)
    ax.set_xticklabels(n_labels, fontsize=10)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", weight="bold", pad=12)
    ax.set_xlim(-0.6, 2.6)

# ------------------------------------------------------------
# PANEL 2: RAREFIED RICHNESS
# ------------------------------------------------------------

draw_box_scatter(
    axes[1],
    "Rarefied_Richness_20",
    "Rarefied species richness",
    "Expected resolved species at n = 20"
)

axes[1].set_ylim(0, 7.25)

# ------------------------------------------------------------
# PANEL 3: SHANNON DIVERSITY
# ------------------------------------------------------------

draw_box_scatter(
    axes[2],
    "Shannon_H",
    "Observed Shannon diversity",
    "Shannon diversity (H)"
)

axes[2].set_ylim(bottom=0)

# ------------------------------------------------------------
# SHARED STYLING
# ------------------------------------------------------------

for ax in axes:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#666666")
    ax.spines["bottom"].set_color("#666666")
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.55, alpha=0.8)
    ax.set_axisbelow(True)

fig.subplots_adjust(
    left=0.075,
    right=0.975,
    top=0.83,
    bottom=0.205
)

fig.text(
    0.5, 0.070,
    "Rarefaction uses 20 resolved genomes per country. "
    "The numbers shown beneath each group indicate the number of eligible countries. "
    "Unresolved ECC annotations were excluded from species-diversity calculations.",
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
# WRITE EXCEL SOURCE DATA
# ============================================================

methods = pd.DataFrame([
    ("Master", "Final 9157 ECC genomes, 2000-2026"),
    ("Species categories", "Seven resolved species; unresolved ECC excluded"),
    ("Primary rarefaction depth", PRIMARY_DEPTH),
    ("Sensitivity rarefaction depth", SENSITIVITY_DEPTH),
    ("Rarefaction method", "Exact expected richness under sampling without replacement"),
    ("Country comparison unit", "Country/territory"),
    ("Country eligibility", "At least 20 resolved genomes"),
    ("Shannon", "Observed Shannon H; not rarefied"),
    ("Statistical comparison", "Exploratory country-level Kruskal-Wallis"),
    ("Figure revision note", "Eligible country counts moved to x-axis bottom"),
    ("Created", datetime.now().isoformat())
], columns=["Parameter", "Value"])

with pd.ExcelWriter(XLSX, engine="openpyxl") as writer:
    rarefaction.to_excel(writer, sheet_name="ALL_COUNTRIES", index=False)
    primary.to_excel(writer, sheet_name="PRIMARY_DEPTH20", index=False)
    sensitivity.to_excel(writer, sheet_name="SENSITIVITY_DEPTH10", index=False)
    eligibility.to_excel(writer, sheet_name="DEPTH_ELIGIBILITY", index=False)
    summary.to_excel(writer, sheet_name="INCOME_SUMMARY", index=False)
    stats.to_excel(writer, sheet_name="COUNTRY_LEVEL_TESTS", index=False)
    methods.to_excel(writer, sheet_name="METHODS", index=False)

# ============================================================
# FINAL VALIDATION + REPORT
# ============================================================

if not JPG.exists():
    raise RuntimeError("Final JPEG not created")

if not XLSX.exists():
    raise RuntimeError("Final Excel not created")

print("\n" + "=" * 95)
print("REVISED FIGURE 2C COMPLETED")
print("=" * 95)

print("\nCOUNTRY-LEVEL DIVERSITY SUMMARY")
print(summary.to_string(index=False))

print("\nEXPLORATORY COUNTRY-LEVEL TESTS")
print(stats.to_string(index=False))

print("\nFINAL OUTPUTS")
print(JPG)
print(XLSX)

print("\nMain revision: group sample sizes moved to bottom of panels 2 and 3")
print("Figure: 600 DPI JPEG")
print("Excel: 7 worksheets")
print("=" * 95)