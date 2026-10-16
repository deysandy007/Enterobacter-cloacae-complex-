
from pathlib import Path
from datetime import datetime
import tempfile
import shutil

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

# ============================================================
# FIGURE 2D - COUNTRY-EQUAL-WEIGHTED ECC SPECIES COMPOSITION
#
# Primary analysis:
#   42 countries with >=20 RESOLVED ECC genomes
#   H=24, M=9, L+LM=9
#
# Compare:
#   1. Genome-weighted composition
#   2. Equal-country-weighted composition
#
# Sensitivity: >=10 resolved genomes per country
#
# Country-bootstrap 95% percentile intervals
# Country-level exploratory permutation test
#
# Input: final ECC master, 9157 genomes, 2000-2026
# Income assignment: frozen Figure 1A workbook
#
# Output:
#   ONE JPEG at 600 DPI
#   ONE Excel workbook
#
# Does not modify the master, supplementary folder,
# or previous manuscript figures.
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

OUT = ROOT / "Fig_02D_Country_Weighted_Species"
OUT.mkdir(parents=True, exist_ok=True)

JPG_NAME = "Fig_02D_Country_Weighted_Species.jpg"
XLSX_NAME = "Fig_02D_Country_Weighted_Data.xlsx"

EXPECTED_N = 9157
PRIMARY_DEPTH = 20
SENSITIVITY_DEPTH = 10

BOOTSTRAP_REPLICATES = 3000
PERMUTATIONS = 4999
SEED = 20261010

RNG = np.random.default_rng(SEED)

GROUPS = ["H", "M", "L+LM"]

GROUP_LABELS = {
    "H": "High income",
    "M": "Upper-middle",
    "L+LM": "Low + lower-middle"
}

INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82"
}

SPECIES = [
    "Enterobacter hormaechei",
    "Enterobacter cloacae",
    "Enterobacter roggenkampii",
    "Enterobacter asburiae",
    "Enterobacter kobei",
    "Enterobacter ludwigii",
    "Enterobacter bugandensis",
    "Enterobacter cloacae complex - unresolved"
]

UNRESOLVED = SPECIES[-1]

SPECIES_LABELS = {
    "Enterobacter hormaechei": "E. hormaechei",
    "Enterobacter cloacae": "E. cloacae",
    "Enterobacter roggenkampii": "E. roggenkampii",
    "Enterobacter asburiae": "E. asburiae",
    "Enterobacter kobei": "E. kobei",
    "Enterobacter ludwigii": "E. ludwigii",
    "Enterobacter bugandensis": "E. bugandensis",
    UNRESOLVED: "Unresolved ECC"
}

SPECIES_COLORS = {
    "Enterobacter hormaechei": "#0072B2",
    "Enterobacter cloacae": "#D55E00",
    "Enterobacter roggenkampii": "#009E73",
    "Enterobacter asburiae": "#CC79A7",
    "Enterobacter kobei": "#7851A9",
    "Enterobacter ludwigii": "#A6761D",
    "Enterobacter bugandensis": "#56B4E9",
    UNRESOLVED: "#808080"
}

# ============================================================
# FONT / GRAPHIC SETTINGS
# ============================================================

installed_fonts = {
    x.name for x in font_manager.fontManager.ttflist
}

FONT = next(
    (
        f for f in [
            "Arial",
            "Helvetica",
            "Liberation Sans",
            "DejaVu Sans"
        ]
        if f in installed_fonts
    ),
    "DejaVu Sans"
)

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "axes.linewidth": 0.7,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})

# ============================================================
# 1. READ / VALIDATE FINAL MASTER
# ============================================================

print("=" * 100)
print("FIGURE 2D - COUNTRY-WEIGHTED ECC SPECIES COMPOSITION")
print("=" * 100)

if not MASTER.exists():
    raise FileNotFoundError(MASTER)

if not INCOME_FILE.exists():
    raise FileNotFoundError(INCOME_FILE)

df = pd.read_excel(
    MASTER,
    sheet_name="FINAL_MASTER_2000_2026",
    engine="openpyxl",
    dtype={"Assembly_Accession": str}
)

required = [
    "Assembly_Accession",
    "Country_Harmonized",
    "Normalized_Species",
    "Verified_Collection_Year"
]

for col in required:
    if col not in df.columns:
        raise ValueError("Missing master column: " + col)

if len(df) != EXPECTED_N:
    raise ValueError(
        f"Expected {EXPECTED_N} genomes, found {len(df)}"
    )

if df["Assembly_Accession"].isna().any():
    raise ValueError("Missing assembly accession")

if df["Assembly_Accession"].duplicated().any():
    raise ValueError("Duplicate assembly accession")

years = pd.to_numeric(
    df["Verified_Collection_Year"],
    errors="coerce"
)

if (
    years.isna().any()
    or not years.between(2000, 2026).all()
    or not (years % 1 == 0).all()
):
    raise ValueError("Invalid collection years")

df["Year"] = years.astype(int)

df["Country"] = (
    df["Country_Harmonized"]
    .fillna("")
    .astype(str)
    .str.strip()
)

df["Species"] = (
    df["Normalized_Species"]
    .fillna("")
    .astype(str)
    .str.strip()
)

if df["Country"].eq("").any():
    raise ValueError("Blank country in master")

unknown_species = sorted(
    set(df["Species"]) - set(SPECIES)
)

if unknown_species:
    raise ValueError(
        "Unexpected normalized species: "
        + str(unknown_species)
    )

# ============================================================
# 2. REUSE FIGURE 1A INCOME ASSIGNMENTS
# ============================================================

print("\nValidating frozen Figure 1A income assignments...")

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
            "Missing Figure 1A income column: " + col
        )

income["Country"] = (
    income["Country_Harmonized"]
    .fillna("")
    .astype(str)
    .str.strip()
)

if income["Country"].duplicated().any():
    raise ValueError("Duplicate country in Figure 1A")

income["Income_Group"] = (
    income["Income_Group"]
    .fillna("Unclassified")
    .astype(str)
    .str.strip()
)

expected_countries = (
    df["Country"].value_counts()
    .rename("Master_Genomes")
)

saved_countries = (
    income.set_index("Country")["Genome_Count"]
    .rename("Figure1A_Genomes")
)

comparison = pd.concat(
    [expected_countries, saved_countries],
    axis=1
)

if comparison.isna().any().any():
    raise ValueError(
        "Figure 1A and final master have "
        "different country labels"
    )

if not (
    comparison["Master_Genomes"].astype(int)
    == comparison["Figure1A_Genomes"].astype(int)
).all():
    raise ValueError(
        "Figure 1A country counts do not "
        "match the final master"
    )

df = df.merge(
    income[["Country", "Income_Group"]],
    on="Country",
    how="left",
    validate="many_to_one"
)

if df["Income_Group"].isna().any():
    raise ValueError("Income merge failed")

if len(df) != EXPECTED_N:
    raise ValueError("Income merge changed genome count")

# ============================================================
# 3. CREATE COUNTRY BY SPECIES MATRIX
# ============================================================

country_counts = pd.crosstab(
    [
        df["Country"],
        df["Income_Group"]
    ],
    df["Species"]
).reindex(
    columns=SPECIES,
    fill_value=0
)

country_counts.columns.name = None

country_counts["Total_Genomes"] = (
    country_counts[SPECIES].sum(axis=1)
)

resolved_species = [
    x for x in SPECIES if x != UNRESOLVED
]

country_counts["Resolved_Genomes"] = (
    country_counts[resolved_species].sum(axis=1)
)

country_counts["Unresolved_Genomes"] = (
    country_counts[UNRESOLVED]
)

country_counts = (
    country_counts
    .reset_index()
    .sort_values(
        ["Income_Group", "Total_Genomes"],
        ascending=[True, False]
    )
    .reset_index(drop=True)
)

if int(country_counts["Total_Genomes"].sum()) != EXPECTED_N:
    raise ValueError("Country-level count mismatch")

print("Total countries:", len(country_counts))
print("Total genomes:", EXPECTED_N)

# ============================================================
# 4. COUNTRY ELIGIBILITY
#
# Primary: >=20 resolved genomes
# Sensitivity: >=10 resolved genomes
#
# Unclassified income remains in the source audit
# but is excluded from the three-group comparison.
# ============================================================

primary_countries = country_counts[
    country_counts["Income_Group"].isin(GROUPS)
    & (
        country_counts["Resolved_Genomes"]
        >= PRIMARY_DEPTH
    )
].copy()

sensitivity_countries = country_counts[
    country_counts["Income_Group"].isin(GROUPS)
    & (
        country_counts["Resolved_Genomes"]
        >= SENSITIVITY_DEPTH
    )
].copy()

primary_n = (
    primary_countries["Income_Group"]
    .value_counts()
    .reindex(GROUPS, fill_value=0)
)

sensitivity_n = (
    sensitivity_countries["Income_Group"]
    .value_counts()
    .reindex(GROUPS, fill_value=0)
)

if primary_n.to_dict() != {
    "H": 24, "M": 9, "L+LM": 9
}:
    raise ValueError(
        "Primary eligible country counts differ "
        "from Figure 2C: "
        + str(primary_n.to_dict())
    )

if sensitivity_n.to_dict() != {
    "H": 26, "M": 14, "L+LM": 12
}:
    raise ValueError(
        "Sensitivity country counts differ "
        "from Figure 2C: "
        + str(sensitivity_n.to_dict())
    )

print("\nPRIMARY COUNTRY COUNTS (DEPTH 20)")
print(primary_n.to_string())

print("\nSENSITIVITY COUNTRY COUNTS (DEPTH 10)")
print(sensitivity_n.to_string())

# ============================================================
# 5. COMPUTE GENOME AND COUNTRY-EQUAL WEIGHTING
# ============================================================

def calculate_weighting(frame, depth_name):

    rows = []
    country_rows = []

    for group in GROUPS:

        subset = frame[
            frame["Income_Group"] == group
        ].copy()

        n_countries = len(subset)

        if n_countries == 0:
            raise ValueError(
                "No countries for income group: " + group
            )

        n_genomes = int(
            subset["Total_Genomes"].sum()
        )

        genome_share = (
            subset[SPECIES].sum(axis=0)
            / n_genomes * 100
        )

        within_country = (
            subset[SPECIES]
            .div(
                subset["Total_Genomes"],
                axis=0
            )
            * 100
        )

        equal_country_share = (
            within_country.mean(axis=0)
        )

        if not np.isclose(
            genome_share.sum(), 100
        ):
            raise RuntimeError(
                "Genome-weighted percentages do not sum to 100"
            )

        if not np.isclose(
            equal_country_share.sum(), 100
        ):
            raise RuntimeError(
                "Country-weighted percentages do not sum to 100"
            )

        for species in SPECIES:

            rows.append({
                "Analysis": depth_name,
                "Income_Group": group,
                "Eligible_Countries": n_countries,
                "Included_Genomes": n_genomes,
                "Species": species,
                "Genome_Weighted_Percent": float(
                    genome_share[species]
                ),
                "Country_Equal_Percent": float(
                    equal_country_share[species]
                ),
                "Difference_Percentage_Points": float(
                    equal_country_share[species]
                    - genome_share[species]
                )
            })

        for index, country_row in subset.iterrows():

            for species in SPECIES:
                country_rows.append({
                    "Analysis": depth_name,
                    "Country": country_row["Country"],
                    "Income_Group": group,
                    "Total_Genomes": int(
                        country_row["Total_Genomes"]
                    ),
                    "Resolved_Genomes": int(
                        country_row["Resolved_Genomes"]
                    ),
                    "Species": species,
                    "Country_Species_Percent": float(
                        within_country.loc[index, species]
                    ),
                    "Country_Weight_Within_Group": (
                        1.0 / n_countries
                    ),
                    "Genome_Weight_Within_Group": (
                        country_row["Total_Genomes"]
                        / n_genomes
                    )
                })

    return (
        pd.DataFrame(rows),
        pd.DataFrame(country_rows)
    )

primary_result, primary_country_data = (
    calculate_weighting(
        primary_countries,
        "PRIMARY_DEPTH20"
    )
)

sensitivity_result, sensitivity_country_data = (
    calculate_weighting(
        sensitivity_countries,
        "SENSITIVITY_DEPTH10"
    )
)

# ============================================================
# 6. BOOTSTRAP COUNTRY-EQUAL COMPOSITION
#
# Resample countries within each income group.
# This estimates between-country sampling variability,
# conditional on the represented country collections.
# ============================================================

print("\nComputing country bootstrap confidence intervals...")

bootstrap_rows = []

for group in GROUPS:

    subset = primary_countries[
        primary_countries["Income_Group"] == group
    ]

    proportions = (
        subset[SPECIES]
        .div(
            subset["Total_Genomes"],
            axis=0
        )
        .to_numpy(dtype=float)
        * 100
    )

    n = len(proportions)

    sampled_indices = RNG.integers(
        0, n,
        size=(BOOTSTRAP_REPLICATES, n)
    )

    boot_means = proportions[
        sampled_indices
    ].mean(axis=1)

    low = np.percentile(
        boot_means, 2.5, axis=0
    )

    high = np.percentile(
        boot_means, 97.5, axis=0
    )

    for j, species in enumerate(SPECIES):
        bootstrap_rows.append({
            "Income_Group": group,
            "Species": species,
            "Bootstrap_95CI_Lower": float(low[j]),
            "Bootstrap_95CI_Upper": float(high[j]),
            "Bootstrap_Replicates": BOOTSTRAP_REPLICATES,
            "Number_of_Countries": n
        })

bootstrap_df = pd.DataFrame(bootstrap_rows)

# ============================================================
# 7. COUNTRY-LEVEL PERMUTATION TEST
#
# Test whether country-level RESOLVED species proportions
# vary among income groups.
#
# One country = one observation.
# The unresolved label is not treated as a species
# for this test.
#
# Euclidean multivariate pseudo-F, with income labels
# permuted across the 42 eligible countries.
#
# This is an exploratory, unadjusted association test,
# not proof of socioeconomic causation.
# ============================================================

print("\nRunning exploratory country-level permutation test...")

permutation_frame = primary_countries.copy()

X = (
    permutation_frame[resolved_species]
    .div(
        permutation_frame["Resolved_Genomes"],
        axis=0
    )
    .to_numpy(dtype=float)
)

labels = permutation_frame[
    "Income_Group"
].to_numpy()

def pseudo_f(values, group_labels):

    n, variables = values.shape
    unique_groups = np.unique(group_labels)
    k = len(unique_groups)

    grand_mean = values.mean(axis=0)

    ss_between = 0.0
    ss_within = 0.0

    for group in unique_groups:

        sub = values[group_labels == group]
        mean = sub.mean(axis=0)

        ss_between += (
            len(sub)
            * np.sum((mean - grand_mean) ** 2)
        )

        ss_within += np.sum(
            (sub - mean) ** 2
        )

    if ss_within <= 0:
        return np.nan

    return (
        (ss_between / (k - 1))
        / (ss_within / (n - k))
    )

observed_f = pseudo_f(X, labels)

if not np.isfinite(observed_f):
    raise RuntimeError(
        "Country-level pseudo-F cannot be calculated"
    )

permuted_f = np.empty(PERMUTATIONS)

for i in range(PERMUTATIONS):

    shuffled = RNG.permutation(labels)

    permuted_f[i] = pseudo_f(
        X, shuffled
    )

permutation_p = (
    1 + int(np.sum(permuted_f >= observed_f))
) / (PERMUTATIONS + 1)

permutation_stats = pd.DataFrame([
    {
        "Test": "Country-label permutation of multivariate pseudo-F",
        "Metric": "Country-level resolved species proportions",
        "Distance": "Squared Euclidean",
        "Pseudo_F": observed_f,
        "Permutation_P": permutation_p,
        "Number_of_Permutations": PERMUTATIONS,
        "Eligible_Countries": len(primary_countries),
        "Group_H": int(primary_n["H"]),
        "Group_M": int(primary_n["M"]),
        "Group_L_LM": int(primary_n["L+LM"]),
        "Interpretation": (
            "Exploratory unadjusted country-level association; "
            "not a causal test"
        )
    }
])

# ============================================================
# 8. COUNTRY SAMPLING IMBALANCE
# ============================================================

imbalance_rows = []

for group in GROUPS:

    sub = primary_countries[
        primary_countries["Income_Group"] == group
    ].sort_values(
        "Total_Genomes",
        ascending=False
    )

    total = int(sub["Total_Genomes"].sum())

    imbalance_rows.append({
        "Income_Group": group,
        "Eligible_Countries": len(sub),
        "Total_Genomes": total,
        "Largest_Country": sub.iloc[0]["Country"],
        "Largest_Country_Genomes": int(
            sub.iloc[0]["Total_Genomes"]
        ),
        "Largest_Country_Percent": (
            100 * sub.iloc[0]["Total_Genomes"]
            / total
        ),
        "Top_3_Country_Percent": (
            100 * sub.head(3)[
                "Total_Genomes"
            ].sum()
            / total
        )
    })

imbalance = pd.DataFrame(imbalance_rows)

# ============================================================
# 9. PREPARE FIGURE TABLES
# ============================================================

genome_plot = (
    primary_result.pivot(
        index="Income_Group",
        columns="Species",
        values="Genome_Weighted_Percent"
    )
    .reindex(
        index=GROUPS,
        columns=SPECIES
    )
)

country_plot = (
    primary_result.pivot(
        index="Income_Group",
        columns="Species",
        values="Country_Equal_Percent"
    )
    .reindex(
        index=GROUPS,
        columns=SPECIES
    )
)

for table in [genome_plot, country_plot]:
    if not np.allclose(
        table.sum(axis=1).to_numpy(),
        100.0
    ):
        raise RuntimeError(
            "Plot proportions do not sum to 100"
        )

# ============================================================
# 10. NATURE-STYLE TWO-PANEL FIGURE
# ============================================================

print("\nRendering Figure 2D...")

fig, axes = plt.subplots(
    1, 2,
    figsize=(16.5, 8.8),
    gridspec_kw={
        "width_ratios": [1, 1],
        "wspace": 0.32
    }
)

fig.suptitle(
    "Country-weighted species composition "
    "of the Enterobacter cloacae complex",
    fontsize=18,
    fontweight="bold",
    y=0.974
)

fig.text(
    0.5, 0.934,
    "Comparison of genome-weighted and "
    "country-equal-weighted taxonomic composition",
    ha="center",
    fontsize=11.5,
    color="#555555"
)

positions = np.arange(len(GROUPS))

labels = [
    f"{GROUP_LABELS[g]}\n"
    f"(countries={int(primary_n[g])})"
    for g in GROUPS
]

def draw_stacked(
    ax,
    table,
    title
):

    left = np.zeros(len(GROUPS))

    for species in SPECIES:

        values = table[species].to_numpy()

        ax.barh(
            positions,
            values,
            left=left,
            height=0.52,
            color=SPECIES_COLORS[species],
            edgecolor="white",
            linewidth=0.4
        )

        for i, value in enumerate(values):

            if value >= 12:
                ax.text(
                    left[i] + value / 2,
                    positions[i],
                    f"{value:.1f}%",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=10,
                    fontweight="bold"
                )

        left += values

    ax.set_title(
        title,
        loc="left",
        fontsize=14,
        fontweight="bold",
        pad=14
    )

    ax.set_xlim(0, 100)

    ax.set_xticks(
        np.arange(0, 101, 20)
    )

    ax.xaxis.set_major_formatter(
        FuncFormatter(
            lambda v, pos: f"{v:.0f}%"
        )
    )

    ax.set_xlabel(
        "Species composition (%)"
    )

    ax.set_yticks(positions)
    ax.set_yticklabels(
        labels,
        fontsize=10.7
    )

    ax.invert_yaxis()

    ax.set_axisbelow(True)

    ax.grid(
        axis="x",
        color="#E5E5E5",
        linewidth=0.6,
        alpha=0.8
    )

    for spine in [
        "left", "right", "top"
    ]:
        ax.spines[spine].set_visible(False)

    ax.spines["bottom"].set_color(
        "#666666"
    )

    ax.tick_params(
        axis="y",
        length=0,
        pad=8
    )

draw_stacked(
    axes[0],
    genome_plot,
    "Genome-weighted composition"
)

draw_stacked(
    axes[1],
    country_plot,
    "Equal-country-weighted composition"
)

# ============================================================
# 11. SHARED SPECIES LEGEND
# ============================================================

handles = [
    Patch(
        facecolor=SPECIES_COLORS[s],
        edgecolor="none",
        label=SPECIES_LABELS[s]
    )
    for s in SPECIES
]

fig.legend(
    handles=handles,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.12),
    ncol=4,
    frameon=False,
    fontsize=10.8,
    columnspacing=1.6,
    handlelength=1.3
)

fig.subplots_adjust(
    left=0.18,
    right=0.96,
    top=0.84,
    bottom=0.25
)

fig.text(
    0.5,
    0.043,
    "The same 42 countries are included in both panels "
    "(minimum 20 resolved genomes per country). "
    "Each country contributes equally in the right panel. "
    "Unresolved ECC annotations are retained separately.",
    ha="center",
    fontsize=9.3,
    color="#555555"
)

# ============================================================
# 12. EXCEL METADATA
# ============================================================

methods = pd.DataFrame([
    (
        "Final input genome master",
        str(MASTER)
    ),
    (
        "Study years",
        "2000-2026"
    ),
    (
        "Original final genome count",
        EXPECTED_N
    ),
    (
        "Income assignments",
        "Frozen Figure 1A country-income database"
    ),
    (
        "Primary inclusion threshold",
        "At least 20 resolved genomes per country"
    ),
    (
        "Sensitivity threshold",
        "At least 10 resolved genomes per country"
    ),
    (
        "Primary countries",
        len(primary_countries)
    ),
    (
        "Primary H / M / L+LM",
        "24 / 9 / 9"
    ),
    (
        "Genome weighting",
        "Sum species genomes across eligible countries "
        "divided by total genomes in that income group"
    ),
    (
        "Equal-country weighting",
        "Calculate within-country species percentages "
        "then take an unweighted mean across countries"
    ),
    (
        "Unresolved categories",
        "Included in descriptive composition; "
        "excluded from resolved-species permutation test"
    ),
    (
        "Bootstrap",
        f"{BOOTSTRAP_REPLICATES} country resamples per group"
    ),
    (
        "Country-level permutation test",
        f"{PERMUTATIONS} label permutations"
    ),
    (
        "Figure format",
        "600-DPI JPEG"
    ),
    (
        "Font",
        FONT
    ),
    (
        "Created",
        datetime.now().isoformat()
    )
], columns=["Parameter", "Value"])

# ============================================================
# 13. WRITE AND VALIDATE IN A STAGING DIRECTORY
# ============================================================

with tempfile.TemporaryDirectory(
    prefix="FIG2D_",
    dir=str(OUT)
) as tmp:

    stage = Path(tmp)

    staged_jpg = stage / JPG_NAME
    staged_xlsx = stage / XLSX_NAME

    fig.savefig(
        staged_jpg,
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

    with pd.ExcelWriter(
        staged_xlsx,
        engine="openpyxl"
    ) as writer:

        primary_result.to_excel(
            writer,
            sheet_name="PRIMARY_COMPARISON",
            index=False
        )

        sensitivity_result.to_excel(
            writer,
            sheet_name="SENSITIVITY_DEPTH10",
            index=False
        )

        primary_countries.to_excel(
            writer,
            sheet_name="ELIGIBLE_COUNTRIES_20",
            index=False
        )

        sensitivity_countries.to_excel(
            writer,
            sheet_name="ELIGIBLE_COUNTRIES_10",
            index=False
        )

        primary_country_data.to_excel(
            writer,
            sheet_name="COUNTRY_SPECIES_PERCENT",
            index=False
        )

        bootstrap_df.to_excel(
            writer,
            sheet_name="BOOTSTRAP_95CI",
            index=False
        )

        permutation_stats.to_excel(
            writer,
            sheet_name="COUNTRY_PERMUTATION",
            index=False
        )

        imbalance.to_excel(
            writer,
            sheet_name="SAMPLING_IMBALANCE",
            index=False
        )

        genome_plot.rename_axis(
            "Income_Group"
        ).reset_index().to_excel(
            writer,
            sheet_name="GENOME_WEIGHTED_PLOT",
            index=False
        )

        country_plot.rename_axis(
            "Income_Group"
        ).reset_index().to_excel(
            writer,
            sheet_name="COUNTRY_WEIGHTED_PLOT",
            index=False
        )

        country_counts.to_excel(
            writer,
            sheet_name="ALL_COUNTRY_COUNTS",
            index=False
        )

        methods.to_excel(
            writer,
            sheet_name="METHODS",
            index=False
        )

    check = pd.read_excel(
        staged_xlsx,
        sheet_name="PRIMARY_COMPARISON"
    )

    if len(check) != len(GROUPS) * len(SPECIES):
        raise RuntimeError(
            "Excel species-income comparison incomplete"
        )

    if staged_jpg.stat().st_size < 100000:
        raise RuntimeError(
            "JPEG unexpectedly small"
        )

    if staged_xlsx.stat().st_size < 10000:
        raise RuntimeError(
            "Excel workbook unexpectedly small"
        )

    shutil.copy2(
        staged_jpg,
        OUT / JPG_NAME
    )

    shutil.copy2(
        staged_xlsx,
        OUT / XLSX_NAME
    )

# ============================================================
# 14. REPORT
# ============================================================

print("\n" + "=" * 100)
print("FIGURE 2D COMPLETED SUCCESSFULLY")
print("=" * 100)

print("\nPRIMARY DEPTH-20 COMPARISON")
print(
    primary_result[
        [
            "Income_Group",
            "Species",
            "Genome_Weighted_Percent",
            "Country_Equal_Percent",
            "Difference_Percentage_Points"
        ]
    ].round(3).to_string(index=False)
)

print("\nCOUNTRY-LEVEL PERMUTATION TEST")
print(
    permutation_stats.to_string(index=False)
)

print("\nSAMPLING IMBALANCE")
print(
    imbalance.round(3).to_string(index=False)
)

print("\nFINAL OUTPUTS")
print(OUT / JPG_NAME)
print(OUT / XLSX_NAME)

print("\n600-DPI JPEG + Excel (12 worksheets)")
print("Final genomic master unchanged")
print("Existing Figures 2A, 2B and 2C unchanged")
print("=" * 100)
