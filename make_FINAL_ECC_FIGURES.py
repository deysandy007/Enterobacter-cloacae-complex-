
from pathlib import Path
from io import BytesIO
from datetime import datetime
import tempfile
import shutil
import math
import re
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import requests
import pycountry

from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

warnings.filterwarnings("ignore", category=FutureWarning)

# ============================================================
# FINAL ECC FIGURES 1A AND 1B - STANDALONE SCRIPT
#
# INPUT: ALREADY NORMALIZED FINAL MASTER (9157 GENOMES)
# NO TAXONOMIC REPROCESSING
# NO SUBPROCESS CALLS OR RECURSION
#
# OUTPUT PER FIGURE:
# ONE 600-DPI JPEG + ONE EXCEL
# ============================================================

PROJECT = Path(
    r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES"
)

MASTER = (
    PROJECT
    / "00_FINAL_COUNTRY_REPORTED_MASTER"
    / "ECC_FINAL_2000_2026_MASTER.xlsx"
)

SUPPLEMENT = (
    PROJECT / "ECC_SUPPLEMENTARY_AND_SOURCE_DATA"
)

FIG_A = (
    PROJECT / "Fig_01A_Global_Geographical_Distribution"
)

FIG_B = (
    PROJECT / "Fig_01B_Temporal_Distribution"
)

A_JPG = "Fig_01A_Global_Geographical_Distribution.jpg"
A_XLSX = "Fig_01A_Country_Income_Database.xlsx"

B_JPG = "Fig_01B_Temporal_Distribution.jpg"
B_XLSX = "Fig_01B_Temporal_Data.xlsx"

INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82",
    "Unclassified": "#999999",
}

INCOME_LABELS = {
    "H": "High income",
    "M": "Upper-middle income",
    "L+LM": "Low + lower-middle income",
    "Unclassified": "Unclassified",
}

INCOME_ORDER = ["H", "M", "L+LM", "Unclassified"]

SOURCE_COLORS = {
    "Human": "#E69F00",
    "Animal": "#56B4E9",
    "Environment": "#009E73",
    "Food": "#D55E00",
    "Plant-associated": "#CC79A7",
}

SOURCE_ORDER = list(SOURCE_COLORS)

SPECIES_COLORS = [
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#7851A9",
    "#A6761D",
    "#56B4E9",
    "#777777",
]

ALIASES = {
    "United States": "USA",
    "United States of America": "USA",
    "United Kingdom": "GBR",
    "Czech Republic": "CZE",
    "Czechia": "CZE",
    "Taiwan": "TWN",
    "Hong Kong": "HKG",
    "South Korea": "KOR",
    "Russia": "RUS",
    "Turkey": "TUR",
    "Türkiye": "TUR",
    "Iran": "IRN",
    "Viet Nam": "VNM",
    "Bolivia": "BOL",
    "Venezuela": "VEN",
    "Tanzania": "TZA",
    "Puerto Rico": "PRI",
    "Reunion": "REU",
    "New Caledonia": "NCL",
}

# Taiwan is a documented supplemental classification,
# rather than a classification returned by the WB API.
SUPPLEMENTAL_INCOME = {"Taiwan": "H"}

WB_GROUPS = {
    "HIC": "H",
    "UMC": "M",
    "LMC": "L+LM",
    "LIC": "L+LM",
}

WORLD_BANK_URL = (
    "https://api.worldbank.org/v2/country"
    "?format=json&per_page=400"
)

WORLD_MAP_URL = (
    "https://raw.githubusercontent.com/"
    "nvkelso/natural-earth-vector/master/"
    "geojson/ne_110m_admin_0_countries.geojson"
)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 15,
    "axes.labelsize": 13,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

def iso3(name):
    name = str(name).strip()
    if name in ALIASES:
        return ALIASES[name]
    try:
        return pycountry.countries.lookup(name).alpha_3
    except LookupError:
        return None

def save_jpeg(fig, path):
    fig.savefig(
        path,
        dpi=600,
        format="jpeg",
        facecolor="white",
        pil_kwargs={
            "quality": 95,
            "subsampling": 0,
            "optimize": True,
        }
    )
    plt.close(fig)

def annual_table(data, field, categories):
    return (
        data.groupby(["Year", field])
        .size()
        .unstack(fill_value=0)
        .reindex(
            index=range(2000, 2027),
            columns=categories,
            fill_value=0
        )
    )

# ============================================================
# 1. READ AND VALIDATE MASTER
# ============================================================

print("=" * 95)
print("FINAL ECC FIGURES 1A + 1B")
print("=" * 95)

if not MASTER.exists():
    raise FileNotFoundError(MASTER)

df = pd.read_excel(
    MASTER,
    sheet_name="FINAL_MASTER_2000_2026",
    dtype={"Assembly_Accession": str},
    engine="openpyxl",
)

required = [
    "Assembly_Accession",
    "Country_Harmonized",
    "Final_Source_Category",
    "Verified_Collection_Year",
    "Normalized_Species",
    "Original_Species",
]

missing = [x for x in required if x not in df.columns]

if missing:
    raise ValueError("Missing master columns: " + str(missing))

N = len(df)

if N != 9157:
    raise ValueError(
        f"Expected 9157 normalized genomes, found {N}"
    )

if df["Assembly_Accession"].duplicated().any():
    raise ValueError("Duplicate assembly accession detected")

year = pd.to_numeric(
    df["Verified_Collection_Year"],
    errors="coerce"
)

if (
    year.isna().any()
    or not year.between(2000, 2026).all()
    or not (year % 1 == 0).all()
):
    raise ValueError("Invalid collection years")

df["Year"] = year.astype(int)

df["Country_Harmonized"] = (
    df["Country_Harmonized"].astype(str).str.strip()
)

df["Normalized_Species"] = (
    df["Normalized_Species"].astype(str).str.strip()
)

df["Source"] = (
    df["Final_Source_Category"].astype(str).str.strip()
)

if not df["Source"].isin(SOURCE_ORDER).all():
    raise ValueError("Unexpected source classification")

print("Final genomes:", N)
print("Collection years:", df.Year.min(), "-", df.Year.max())
print("Reported species:", df.Normalized_Species.nunique())

# ============================================================
# 2. COUNTRY COUNTS AND WORLD BANK INCOME
# ============================================================

print("\nRetrieving income classifications...")

countries = (
    df.groupby("Country_Harmonized")
    .size()
    .reset_index(name="Genome_Count")
    .sort_values(
        "Genome_Count",
        ascending=False
    )
    .reset_index(drop=True)
)

countries["ISO3"] = (
    countries["Country_Harmonized"].map(iso3)
)

session = requests.Session()

r = session.get(WORLD_BANK_URL, timeout=90)
r.raise_for_status()

data = r.json()

if not isinstance(data, list) or len(data) < 2:
    raise RuntimeError("Invalid World Bank response")

wb = {}

for x in data[1]:
    code = (x.get("incomeLevel") or {}).get("id")
    group = WB_GROUPS.get(code)
    if group:
        wb[str(x["id"]).upper()] = {
            "group": group,
            "code": code,
            "longitude": x.get("longitude"),
            "latitude": x.get("latitude"),
        }

if len(wb) < 180:
    raise RuntimeError("Incomplete World Bank response")

def assign_income(row):
    code = row["ISO3"]
    name = row["Country_Harmonized"]

    if code in wb:
        return pd.Series({
            "Income_Group": wb[code]["group"],
            "Income_Status": "World Bank API",
            "WB_Code": wb[code]["code"],
        })

    if name in SUPPLEMENTAL_INCOME:
        return pd.Series({
            "Income_Group": SUPPLEMENTAL_INCOME[name],
            "Income_Status": "Supplemental assignment",
            "WB_Code": "",
        })

    return pd.Series({
        "Income_Group": "Unclassified",
        "Income_Status": "No verified assignment",
        "WB_Code": "",
    })

countries = pd.concat(
    [countries, countries.apply(assign_income, axis=1)],
    axis=1
)

unclassified_n = int(
    countries.loc[
        countries["Income_Group"].eq("Unclassified"),
        "Genome_Count"
    ].sum()
)

print("Unclassified-income genomes:", unclassified_n)

if unclassified_n > N * 0.10:
    raise RuntimeError(
        "More than 10% of genomes lack income mapping. "
        "Inspect country names before plotting."
    )

income_join = countries[
    [
        "Country_Harmonized",
        "ISO3",
        "Income_Group",
        "Income_Status",
        "WB_Code",
    ]
]

df = df.merge(
    income_join,
    on="Country_Harmonized",
    how="left",
    validate="many_to_one"
)

assert len(df) == N

income_summary = (
    df["Income_Group"]
    .value_counts()
    .reindex(INCOME_ORDER, fill_value=0)
    .rename_axis("Income_Group")
    .reset_index(name="Genome_Count")
)

income_summary["Percent"] = (
    income_summary["Genome_Count"] / N * 100
).round(3)

# ============================================================
# 3. LOAD WORLD MAP
# ============================================================

print("Retrieving Natural Earth world boundaries...")

r = session.get(WORLD_MAP_URL, timeout=120)
r.raise_for_status()

with tempfile.TemporaryDirectory() as tmp:
    map_path = Path(tmp) / "world.geojson"
    map_path.write_bytes(r.content)
    world = gpd.read_file(map_path)

iso_col = next(
    (
        c for c in ["ADM0_A3", "ISO_A3_EH", "ISO_A3"]
        if c in world.columns
    ),
    None
)

if iso_col is None:
    raise ValueError("World boundaries have no ISO3 field")

world["ISO3"] = world[iso_col].astype(str)
world = world[world["ISO3"] != "ATA"].copy()

map_income = dict(
    zip(countries["ISO3"], countries["Income_Group"])
)

world["Income_Group"] = (
    world["ISO3"].map(map_income).fillna("No data")
)

# ============================================================
# 4. COUNTRY COORDINATES
# ============================================================

coords = {}

for code, x in wb.items():
    try:
        lon = float(x["longitude"])
        lat = float(x["latitude"])
        if np.isfinite(lon) and np.isfinite(lat):
            if not (lon == 0 and lat == 0):
                coords[code] = (lon, lat)
    except (ValueError, TypeError):
        pass

for _, row in world.iterrows():
    code = row["ISO3"]

    if code not in coords:
        point = row.geometry.representative_point()
        coords[code] = (point.x, point.y)

coords["TWN"] = (121.0, 23.7)
coords["HKG"] = (114.2, 22.3)

countries["Longitude"] = countries["ISO3"].map(
    lambda x: coords.get(x, (np.nan, np.nan))[0]
)

countries["Latitude"] = countries["ISO3"].map(
    lambda x: coords.get(x, (np.nan, np.nan))[1]
)

countries["Map_Status"] = np.where(
    countries["Longitude"].notna()
    & countries["Latitude"].notna(),
    "Mapped",
    "Unmapped"
)

def bubble_area(n):
    return 25 + 24 * np.sqrt(np.asarray(n, dtype=float))

# ============================================================
# 5. CREATE FIGURE 1A
# ============================================================

def make_figure_a(output):

    print("\nCreating Figure 1A...")

    fig = plt.figure(figsize=(20, 16.5))

    fig.text(
        0.035, 0.973,
        "Global geographical distribution of "
        "Enterobacter cloacae complex genomes",
        fontsize=20,
        weight="bold"
    )

    fig.text(
        0.035, 0.946,
        f"{N:,} genomes | Collection years 2000–2026",
        fontsize=14,
        color="#555555"
    )

    ax = fig.add_axes([0.025, 0.48, 0.95, 0.45])

    world.plot(
        ax=ax,
        color="#EEEEEE",
        edgecolor="#C4C4C4",
        linewidth=0.35,
    )

    for group in INCOME_ORDER:
        sub = world[world["Income_Group"].eq(group)]
        if len(sub):
            sub.plot(
                ax=ax,
                color=INCOME_COLORS[group],
                edgecolor="#777777",
                linewidth=0.35,
                alpha=0.68,
            )

    mapped = countries.dropna(
        subset=["Longitude", "Latitude"]
    )

    for group in INCOME_ORDER:
        sub = mapped[
            mapped["Income_Group"].eq(group)
        ]
        if len(sub):
            ax.scatter(
                sub.Longitude,
                sub.Latitude,
                s=bubble_area(sub.Genome_Count),
                color=INCOME_COLORS[group],
                edgecolor="white",
                linewidth=0.8,
                alpha=0.95,
                zorder=10
            )

    for row in countries.head(10).itertuples():
        if not np.isfinite(row.Longitude):
            continue

        ax.annotate(
            f"{row.Country_Harmonized}\n"
            f"{row.Genome_Count:,}",
            (row.Longitude, row.Latitude),
            xytext=(0, 13),
            textcoords="offset points",
            ha="center",
            fontsize=10,
            weight="bold",
            bbox={
                "facecolor": "white",
                "edgecolor": "#CCCCCC",
                "alpha": 0.94,
                "pad": 2,
            },
            zorder=20,
        )

    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 85)
    ax.set_axis_off()

    handles = [
        Patch(
            facecolor=INCOME_COLORS[g],
            label=INCOME_LABELS[g]
        )
        for g in INCOME_ORDER
    ]

    fig.legend(
        handles=handles,
        ncol=4,
        frameon=False,
        fontsize=12,
        loc="center",
        bbox_to_anchor=(0.5, 0.47),
    )

    fig.text(
        0.5, 0.442,
        "Country-income classification: "
        "World Bank API at analysis date",
        ha="center",
        fontsize=12.5,
        weight="bold"
    )

    # Income-stratified lists under the map.
    for group, left in [
        ("H", 0.03),
        ("M", 0.355),
        ("L+LM", 0.68),
    ]:
        axp = fig.add_axes([
            left, 0.12, 0.29, 0.29
        ])

        axp.set_xlim(0, 1)
        axp.set_ylim(0, 1)
        axp.set_xticks([])
        axp.set_yticks([])

        axp.set_facecolor("#FAFAFA")

        for spine in axp.spines.values():
            spine.set_edgecolor("#DDDDDD")

        axp.axhspan(
            0.92, 1,
            color=INCOME_COLORS[group]
        )

        axp.text(
            0.025, 0.96,
            INCOME_LABELS[group],
            color="white",
            fontsize=12,
            weight="bold",
            va="center"
        )

        sub = countries[
            countries.Income_Group.eq(group)
        ].sort_values(
            "Genome_Count",
            ascending=False
        )

        axp.text(
            0.025, 0.90,
            f"{len(sub)} locations | "
            f"{sub.Genome_Count.sum():,} genomes",
            va="top",
            fontsize=10,
            weight="bold"
        )

        entries = [
            f"{r.Country_Harmonized} "
            f"({r.Genome_Count:,})"
            for r in sub.itertuples()
        ]

        rows = max(1, math.ceil(len(entries) / 2))
        font_size = 10 if rows <= 25 else 8.5

        for i, label in enumerate(entries):
            col = i // rows
            row = i % rows

            axp.text(
                0.02 + col * 0.51,
                0.85 - row * (0.80 / rows),
                label,
                fontsize=font_size,
                va="top",
                clip_on=True,
            )

    fig.text(
        0.035, 0.084,
        "Genome-count bubble scale:",
        fontsize=12,
        weight="bold"
    )

    legax = fig.add_axes([0.29, 0.065, 0.65, 0.07])
    legax.set_xlim(0, 6)
    legax.set_ylim(0, 2)
    legax.axis("off")

    for i, count in enumerate([10, 100, 500, 1000, 2500]):
        x = 0.55 + i * 1.2
        legax.scatter(
            x, 1.15,
            s=bubble_area(count),
            color="#777777",
            edgecolor="white"
        )
        legax.text(
            x, 0.1,
            f"{count:,}",
            ha="center",
            fontsize=11
        )

    fig.text(
        0.035, 0.035,
        "Bubble sizes represent numbers of available "
        "genomes, not disease incidence. "
        "Income assignments and unmatched locations "
        "are documented in the source-data workbook.",
        fontsize=10,
        color="#555555"
    )

    save_jpeg(fig, output / A_JPG)

    with pd.ExcelWriter(
        output / A_XLSX,
        engine="openpyxl"
    ) as writer:

        countries.to_excel(
            writer,
            sheet_name="COUNTRY_INCOME_COUNTS",
            index=False
        )

        income_summary.to_excel(
            writer,
            sheet_name="INCOME_SUMMARY",
            index=False
        )

        df.groupby([
            "Country_Harmonized", "Source"
        ]).size().reset_index(
            name="Genome_Count"
        ).to_excel(
            writer,
            sheet_name="COUNTRY_SOURCE_COUNTS",
            index=False
        )

        df[[
            "Assembly_Accession",
            "Country_Harmonized",
            "Year",
            "Normalized_Species",
            "Source",
            "Income_Group",
            "Income_Status",
        ]].to_excel(
            writer,
            sheet_name="GENOME_SOURCE_DATA",
            index=False
        )

        pd.DataFrame([
            ("Final genomes", N),
            ("Collection years", "2000-2026"),
            ("Income source", WORLD_BANK_URL),
            ("Map source", WORLD_MAP_URL),
            ("Date", datetime.now().isoformat()),
        ], columns=[
            "Parameter", "Value"
        ]).to_excel(
            writer,
            sheet_name="METHODOLOGY",
            index=False
        )

# ============================================================
# 6. CREATE FIGURE 1B
# ============================================================

def make_figure_b(output):

    print("Creating Figure 1B...")

    years = np.arange(2000, 2027)

    annual = (
        df.groupby("Year")
        .size()
        .reindex(years, fill_value=0)
    )

    species_order = (
        df["Normalized_Species"]
        .value_counts()
        .index.tolist()
    )

    by_source = annual_table(
        df, "Source", SOURCE_ORDER
    )

    by_species = annual_table(
        df, "Normalized_Species", species_order
    )

    by_income = annual_table(
        df, "Income_Group", INCOME_ORDER
    )

    assert annual.sum() == N
    assert by_source.values.sum() == N
    assert by_species.values.sum() == N
    assert by_income.values.sum() == N

    species_colors = {
        x: SPECIES_COLORS[i % len(SPECIES_COLORS)]
        for i, x in enumerate(species_order)
    }

    fig, axes = plt.subplots(
        4, 1,
        figsize=(18, 19),
        sharex=True,
        gridspec_kw={
            "height_ratios": [1, 1.3, 1.3, 1.3],
            "hspace": 0.42
        }
    )

    fig.suptitle(
        "Temporal distribution of Enterobacter cloacae "
        "complex genomes",
        fontsize=20,
        weight="bold",
        y=0.983
    )

    fig.text(
        0.5, 0.955,
        f"{N:,} genomes | Collection years 2000–2026",
        ha="center",
        fontsize=13,
        color="#555555"
    )

    axes[0].bar(
        years,
        annual.values,
        color="#444444",
        width=0.85
    )

    axes[0].set_title(
        "Overall annual genome collection",
        loc="left",
        weight="bold",
        pad=12
    )

    def stacked(ax, table, colors, labels):
        bottom = np.zeros(len(years))

        for name in table.columns:
            values = table[name].values
            ax.bar(
                years,
                values,
                bottom=bottom,
                width=0.85,
                color=colors[name],
                label=labels.get(name, name),
                edgecolor="none",
            )
            bottom += values

        ax.legend(
            loc="upper left",
            bbox_to_anchor=(0, 1.01),
            ncol=4,
            frameon=False,
            fontsize=10
        )

    stacked(
        axes[1], by_source,
        SOURCE_COLORS, {}
    )

    axes[1].set_title(
        "Temporal distribution by One Health source",
        loc="left",
        weight="bold",
        pad=12
    )

    stacked(
        axes[2], by_species,
        species_colors,
        {
            x: x.replace("Enterobacter ", "E. ")
            for x in species_order
        }
    )

    axes[2].set_title(
        "Temporal distribution by normalized species",
        loc="left",
        weight="bold",
        pad=12
    )

    stacked(
        axes[3], by_income,
        INCOME_COLORS, INCOME_LABELS
    )

    axes[3].set_title(
        "Temporal distribution by country-income group",
        loc="left",
        weight="bold",
        pad=12
    )

    for ax in axes:
        ax.set_xlim(1999.4, 2026.6)
        ax.set_ylabel("Genomes")
        ax.set_axisbelow(True)
        ax.grid(
            axis="y",
            color="#DDDDDD",
            linewidth=0.6
        )
        ax.yaxis.set_major_locator(
            MaxNLocator(nbins=5, integer=True)
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[-1].set_xlabel("Collection year")
    axes[-1].set_xticks(
        [2000, 2005, 2010, 2015, 2020, 2026]
    )

    fig.subplots_adjust(
        left=0.09,
        right=0.97,
        top=0.925,
        bottom=0.07
    )

    fig.text(
        0.5, 0.025,
        "Counts reflect genomes available in public "
        "databases, not population incidence or prevalence.",
        ha="center",
        fontsize=11,
        color="#555555"
    )

    save_jpeg(fig, output / B_JPG)

    with pd.ExcelWriter(
        output / B_XLSX,
        engine="openpyxl"
    ) as writer:

        annual.rename(
            "Genome_Count"
        ).reset_index().rename(
            columns={"index": "Year"}
        ).to_excel(
            writer,
            sheet_name="ANNUAL_COUNTS",
            index=False
        )

        for name, table, category in [
            ("YEAR_BY_SOURCE", by_source, "Source"),
            ("YEAR_BY_SPECIES", by_species, "Species"),
            ("YEAR_BY_INCOME", by_income, "Income_Group"),
        ]:
            table.rename_axis("Year").reset_index().melt(
                id_vars="Year",
                var_name=category,
                value_name="Genome_Count"
            ).to_excel(
                writer,
                sheet_name=name,
                index=False
            )

        df[[
            "Assembly_Accession",
            "Year",
            "Country_Harmonized",
            "Source",
            "Original_Species",
            "Normalized_Species",
            "Income_Group",
        ]].to_excel(
            writer,
            sheet_name="GENOME_SOURCE_DATA",
            index=False
        )

# ============================================================
# 7. BUILD AND VALIDATE IN TEMPORARY FOLDERS
# ============================================================

print("\nBuilding both figures safely...")

with tempfile.TemporaryDirectory(
    prefix="ECC_FIGURES_",
    dir=str(PROJECT)
) as tmp:

    staging = Path(tmp)

    stage_a = staging / FIG_A.name
    stage_b = staging / FIG_B.name

    stage_a.mkdir()
    stage_b.mkdir()

    make_figure_a(stage_a)
    make_figure_b(stage_b)

    # Verify counts against the normalized master.
    a_data = pd.read_excel(
        stage_a / A_XLSX,
        sheet_name="COUNTRY_INCOME_COUNTS"
    )

    b_data = pd.read_excel(
        stage_b / B_XLSX,
        sheet_name="ANNUAL_COUNTS"
    )

    assert int(a_data.Genome_Count.sum()) == N
    assert int(b_data.Genome_Count.sum()) == N

    for filename in [
        stage_a / A_JPG,
        stage_b / B_JPG
    ]:
        if not filename.exists():
            raise RuntimeError("Missing output: " + str(filename))

        if filename.stat().st_size < 100000:
            raise RuntimeError(
                "Unexpectedly small JPEG: " + str(filename)
            )

    print("Figure validation passed.")

    # Replace previous figure folders only after validation.
    for target, source in [
        (FIG_A, stage_a),
        (FIG_B, stage_b),
    ]:
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(source, target)

# ============================================================
# 8. UPDATE SUPPLEMENTARY FOLDER
# ============================================================

print("\nUpdating supplementary master...")

SUPPLEMENT.mkdir(exist_ok=True)

if (
    SUPPLEMENT.name != "ECC_SUPPLEMENTARY_AND_SOURCE_DATA"
    or SUPPLEMENT.parent.resolve() != PROJECT.resolve()
):
    raise RuntimeError("Unexpected supplementary folder")

if list(SUPPLEMENT.rglob("~$*.xlsx")):
    raise RuntimeError(
        "Close any open supplementary Excel workbook"
    )

# The final master is safe in the primary master folder.
for item in SUPPLEMENT.iterdir():
    if item.is_dir() and not item.is_symlink():
        shutil.rmtree(item)
    else:
        item.unlink()

shutil.copy2(
    MASTER,
    SUPPLEMENT / MASTER.name
)

assert (
    len(list(SUPPLEMENT.iterdir())) == 1
)

# ============================================================
# 9. FINAL REPORT
# ============================================================

print("\n" + "=" * 95)
print("FINAL ECC FIGURES COMPLETED SUCCESSFULLY")
print("=" * 95)

print("\nGENOMES:", N)
print("YEARS: 2000-2026")

print("\nNORMALIZED SPECIES:")
print(
    df.Normalized_Species.value_counts().to_string()
)

print("\nINCOME DISTRIBUTION:")
print(income_summary.to_string(index=False))

print("\nFIGURE 1A:")
for x in FIG_A.iterdir():
    print(x)

print("\nFIGURE 1B:")
for x in FIG_B.iterdir():
    print(x)

print("\nSUPPLEMENTARY:")
for x in SUPPLEMENT.iterdir():
    print(x)

print("\nMASTER WAS NOT REWRITTEN")
print("NO RECURSIVE SCRIPT CALLS")
print("TWO FIGURES, EACH WITH ONE JPEG AND ONE EXCEL")
print("=" * 95)
