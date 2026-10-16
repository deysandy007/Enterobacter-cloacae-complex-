
from pathlib import Path
from datetime import datetime
import tempfile
import shutil
import hashlib
import math
import json
import warnings

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
from matplotlib import font_manager
import geopandas as gpd
import requests
import pycountry

warnings.filterwarnings("ignore", category=FutureWarning)

# ============================================================
# ECC: FINAL FIGURES 1A AND 1B
# Nature Microbiology-oriented visual refinement
#
# Final normalized genomic master: 9157 assemblies
# Years: 2000-2026 inclusive
#
# Figure 1A:
#   World map, income classes, country count bubbles,
#   ALL country names under income class headings,
#   explicit income thresholds and bubble legend
#
# Figure 1B:
#   Annual counts, One Health, normalized species,
#   and income categories
#   Thin bars, consistent axes, clear legends
#
# Outputs per figure: ONE JPEG (600 DPI), ONE Excel
# Supplementary folder: ONLY final normalized master
# ============================================================

PROJECT = Path(
    r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES"
)
MASTER_DIR = PROJECT / "00_FINAL_COUNTRY_REPORTED_MASTER"
MASTER = MASTER_DIR / "ECC_FINAL_2000_2026_MASTER.xlsx"
SUPPLEMENT = PROJECT / "ECC_SUPPLEMENTARY_AND_SOURCE_DATA"

FOLDER_A = PROJECT / "Fig_01A_Global_Geographical_Distribution"
FOLDER_B = PROJECT / "Fig_01B_Temporal_Distribution"

FILE_A = "Fig_01A_Global_Geographical_Distribution.jpg"
DATA_A = "Fig_01A_Country_Income_Database.xlsx"
FILE_B = "Fig_01B_Temporal_Distribution.jpg"
DATA_B = "Fig_01B_Temporal_Data.xlsx"

EXPECTED_N = 9157
BAR_WIDTH = 0.58
START_YEAR = 2000
END_YEAR = 2026

INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82",
    "Unclassified": "#999999"
}
INCOME_ORDER = ["H", "M", "L+LM", "Unclassified"]
INCOME_NAMES = {
    "H": "HIGH INCOME",
    "M": "UPPER-MIDDLE INCOME",
    "L+LM": "LOW + LOWER-MIDDLE INCOME",
    "Unclassified": "UNCLASSIFIED INCOME"
}

SOURCE_COLORS = {
    "Human": "#E69F00",
    "Animal": "#56B4E9",
    "Environment": "#009E73",
    "Food": "#D55E00",
    "Plant-associated": "#CC79A7"
}
SOURCE_ORDER = list(SOURCE_COLORS)

SPECIES_COLORS = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#7851A9", "#A6761D", "#56B4E9", "#777777"
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
    "North Korea": "PRK",
    "Russia": "RUS",
    "Turkey": "TUR",
    "Türkiye": "TUR",
    "Iran": "IRN",
    "Viet Nam": "VNM",
    "Vietnam": "VNM",
    "Bolivia": "BOL",
    "Venezuela": "VEN",
    "Tanzania": "TZA",
    "Puerto Rico": "PRI",
    "Reunion": "REU",
    "New Caledonia": "NCL",
    "Ivory Coast": "CIV",
    "Côte d'Ivoire": "CIV",
    "Democratic Republic of the Congo": "COD",
    "Republic of the Congo": "COG",
    "Palestine": "PSE",
    "Kosovo": "XKX",
    "Svalbard": "SJM"
}

SUPPLEMENTAL_INCOME = {
    "Taiwan": "H"
}

WB_GROUP = {
    "HIC": "H",
    "UMC": "M",
    "LMC": "L+LM",
    "LIC": "L+LM"
}

WB_URL = (
    "https://api.worldbank.org/v2/country"
    "?format=json&per_page=400"
)

MAP_URL = (
    "https://raw.githubusercontent.com/"
    "nvkelso/natural-earth-vector/master/"
    "geojson/ne_110m_admin_0_countries.geojson"
)

# ============================================================
# PUBLICATION FONT
# ============================================================

available_fonts = {
    x.name for x in font_manager.fontManager.ttflist
}
FONT = next(
    (
        f for f in [
            "Arial", "Helvetica",
            "Liberation Sans", "DejaVu Sans"
        ]
        if f in available_fonts
    ),
    "DejaVu Sans"
)

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.linewidth": 0.7,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})

# ============================================================
# COMMON FUNCTIONS
# ============================================================

def iso3(country):
    country = str(country).strip()
    if country in ALIASES:
        return ALIASES[country]
    try:
        return pycountry.countries.lookup(country).alpha_3
    except LookupError:
        return None


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024), b""
        ):
            digest.update(chunk)
    return digest.hexdigest()


def save_600(fig, destination):
    fig.savefig(
        destination,
        dpi=600,
        format="jpeg",
        facecolor="white",
        pil_kwargs={
            "quality": 94,
            "subsampling": 0,
            "optimize": True
        }
    )
    plt.close(fig)


def bubble_size(n):
    return 18 + 24 * np.sqrt(
        np.asarray(n, dtype=float)
    )


def yearly_counts(data, col, categories):
    return (
        data.groupby(["Year", col])
        .size()
        .unstack(fill_value=0)
        .reindex(
            index=range(START_YEAR, END_YEAR + 1),
            columns=categories,
            fill_value=0
        )
    )


def format_axes(ax):
    ax.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#444444")
    ax.spines["bottom"].set_color("#444444")
    ax.tick_params(
        direction="out", length=3.5,
        width=0.7, colors="#333333"
    )
    ax.yaxis.set_major_locator(
        ticker.MaxNLocator(nbins=5, integer=True)
    )
    ax.yaxis.set_major_formatter(
        ticker.StrMethodFormatter("{x:,.0f}")
    )
    ax.set_axisbelow(True)
    ax.set_ylabel("Genomes")


# ============================================================
# 1. VALIDATE THE FINAL MASTER
# ============================================================

print("=" * 90)
print("ECC FIGURE 1A + 1B FINAL NATURE-STYLE REBUILD")
print("=" * 90)

if not MASTER.exists():
    raise FileNotFoundError(MASTER)

df = pd.read_excel(
    MASTER,
    sheet_name="FINAL_MASTER_2000_2026",
    engine="openpyxl",
    dtype={"Assembly_Accession": str}
)

required = [
    "Assembly_Accession",
    "Country_Harmonized",
    "Verified_Collection_Year",
    "Final_Source_Category",
    "Normalized_Species",
    "Original_Species"
]

missing_cols = [
    c for c in required if c not in df.columns
]
if missing_cols:
    raise ValueError(
        "Missing columns: " + str(missing_cols)
    )

N = len(df)

if N != EXPECTED_N:
    raise ValueError(
        f"Expected {EXPECTED_N} genomes, found {N}"
    )

if df["Assembly_Accession"].duplicated().any():
    raise ValueError("Duplicate assembly accessions")

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
df["Country_Harmonized"] = (
    df["Country_Harmonized"]
    .fillna("").astype(str).str.strip()
)
df["Source"] = (
    df["Final_Source_Category"]
    .fillna("").astype(str).str.strip()
)
df["Normalized_Species"] = (
    df["Normalized_Species"]
    .fillna("").astype(str).str.strip()
)

if not df["Source"].isin(SOURCE_ORDER).all():
    raise ValueError("Unexpected One Health sources")

if (
    df["Country_Harmonized"]
    .str.lower()
    .isin(["", "nan", "none", "unknown"])
    .any()
):
    raise ValueError("Missing countries in final master")

print("Verified final genomes:", N)
print("Years:", df.Year.min(), "to", df.Year.max())
print("Normalized taxa:", df.Normalized_Species.nunique())

# ============================================================
# 2. SYNCHRONIZE SUPPLEMENTARY MASTER FIRST
#
# The supplementary folder will contain exactly one file:
# ECC_FINAL_2000_2026_MASTER.xlsx
# ============================================================

print("\nSynchronizing supplementary master...")

if (
    SUPPLEMENT.name != "ECC_SUPPLEMENTARY_AND_SOURCE_DATA"
    or SUPPLEMENT.parent.resolve() != PROJECT.resolve()
):
    raise RuntimeError("Unexpected supplementary folder")

SUPPLEMENT.mkdir(parents=True, exist_ok=True)

temporary_copy = (
    MASTER_DIR / "_ECC_SUPPLEMENT_VERIFIED_COPY.xlsx"
)

shutil.copy2(MASTER, temporary_copy)

if sha256(MASTER) != sha256(temporary_copy):
    raise RuntimeError("Master staging checksum mismatch")

if list(SUPPLEMENT.rglob("~$*.xlsx")):
    raise RuntimeError(
        "Close any open supplementary Excel files first"
    )

for existing in SUPPLEMENT.iterdir():
    if existing.is_dir() and not existing.is_symlink():
        shutil.rmtree(existing)
    else:
        existing.unlink()

supp_master = SUPPLEMENT / MASTER.name
shutil.copy2(temporary_copy, supp_master)
temporary_copy.unlink()

assert len(list(SUPPLEMENT.iterdir())) == 1
assert sha256(MASTER) == sha256(supp_master)

print("Supplementary folder: one verified master only")

# ============================================================
# 3. COUNTRY METADATA AND INCOME CLASSIFICATION
# ============================================================

countries = (
    df.groupby("Country_Harmonized")
    .size()
    .reset_index(name="Genome_Count")
    .sort_values(
        ["Genome_Count", "Country_Harmonized"],
        ascending=[False, True]
    )
    .reset_index(drop=True)
)

countries["ISO3"] = (
    countries["Country_Harmonized"].map(iso3)
)

print("\nDownloading World Bank income metadata...")

session = requests.Session()
session.headers.update({
    "User-Agent": "ECC-Geographical-Figure/1.0"
})

r = session.get(WB_URL, timeout=90)
r.raise_for_status()
payload = r.json()

if not isinstance(payload, list) or len(payload) < 2:
    raise RuntimeError("World Bank API returned invalid data")

wb = {}

for item in payload[1]:
    level = (item.get("incomeLevel") or {}).get("id")
    group = WB_GROUP.get(level)
    if not group:
        continue
    code = str(item.get("id", "")).upper()
    wb[code] = {
        "Income_Group": group,
        "WB_Code": level,
        "Longitude": item.get("longitude"),
        "Latitude": item.get("latitude")
    }

if len(wb) < 180:
    raise RuntimeError("World Bank income metadata incomplete")

def classify_country(row):
    code = row["ISO3"]
    name = row["Country_Harmonized"]

    if code in wb:
        return pd.Series({
            "Income_Group": wb[code]["Income_Group"],
            "Income_Source": "World Bank API",
            "WB_Code": wb[code]["WB_Code"]
        })

    if name in SUPPLEMENTAL_INCOME:
        return pd.Series({
            "Income_Group": SUPPLEMENTAL_INCOME[name],
            "Income_Source": "Supplemental classification",
            "WB_Code": ""
        })

    return pd.Series({
        "Income_Group": "Unclassified",
        "Income_Source": "Unclassified",
        "WB_Code": ""
    })

countries = pd.concat(
    [countries, countries.apply(classify_country, axis=1)],
    axis=1
)

unclassified = countries.loc[
    countries["Income_Group"] == "Unclassified"
]

print("\nIncome classification:")
print(
    countries.groupby("Income_Group")["Genome_Count"]
    .sum()
    .to_string()
)

if (
    unclassified["Genome_Count"].sum()
    > N * 0.10
):
    raise RuntimeError(
        "More than 10% of genomes lack income classification. "
        "Inspect country-income mapping before plotting."
    )

lookup = countries[
    [
        "Country_Harmonized",
        "ISO3",
        "Income_Group",
        "Income_Source",
        "WB_Code"
    ]
]

data = df.merge(
    lookup,
    on="Country_Harmonized",
    how="left",
    validate="many_to_one"
)

if data["Income_Group"].isna().any():
    raise ValueError("Country income merge failed")

income_summary = (
    countries.groupby("Income_Group")
    .agg(
        Countries=("Country_Harmonized", "size"),
        Genomes=("Genome_Count", "sum")
    )
    .reindex(INCOME_ORDER, fill_value=0)
    .reset_index()
)

income_summary["Percent"] = (
    income_summary["Genomes"] / N * 100
).round(3)

# ============================================================
# 4. LOAD WORLD MAP AND BUBBLE CENTERS
# ============================================================

print("\nDownloading Natural Earth map...")

r = session.get(MAP_URL, timeout=120)
r.raise_for_status()

with tempfile.TemporaryDirectory() as temp:
    path = Path(temp) / "world.geojson"
    path.write_bytes(r.content)
    world = gpd.read_file(path)

iso_col = next(
    (
        c for c in [
            "ADM0_A3", "ISO_A3_EH", "ISO_A3"
        ]
        if c in world.columns
    ),
    None
)

if iso_col is None:
    raise RuntimeError("World map ISO code absent")

world["ISO3"] = world[iso_col].astype(str)
world = world[world["ISO3"] != "ATA"].copy()

income_map = dict(
    zip(countries["ISO3"], countries["Income_Group"])
)
world["Income_Group"] = (
    world["ISO3"].map(income_map).fillna("No data")
)

coordinates = {}

for code, record in wb.items():
    try:
        lon = float(record["Longitude"])
        lat = float(record["Latitude"])
        if (
            np.isfinite(lon)
            and np.isfinite(lat)
            and not (lon == 0 and lat == 0)
        ):
            coordinates[code] = (lon, lat)
    except (TypeError, ValueError):
        pass

for _, row in world.iterrows():
    code = row["ISO3"]
    if code not in coordinates:
        point = row.geometry.representative_point()
        coordinates[code] = (point.x, point.y)

coordinates["TWN"] = (121.0, 23.7)
coordinates["HKG"] = (114.17, 22.32)
coordinates["PRI"] = (-66.1, 18.2)
coordinates["REU"] = (55.5, -21.1)
coordinates["NCL"] = (165.6, -21.5)

countries["Longitude"] = countries["ISO3"].map(
    lambda c: coordinates.get(
        c, (np.nan, np.nan)
    )[0]
)
countries["Latitude"] = countries["ISO3"].map(
    lambda c: coordinates.get(
        c, (np.nan, np.nan)
    )[1]
)

countries["Map_Status"] = np.where(
    countries["Longitude"].notna()
    & countries["Latitude"].notna(),
    "Mapped", "No coordinates"
)

# ============================================================
# 5. FIGURE 1A
# ============================================================

def draw_country_list(ax, frame, group, title):
    """
    Render EVERY country assigned to an income group.
    Do not truncate or silently omit countries.
    """

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor("#FBFBFB")

    for spine in ax.spines.values():
        spine.set_color("#D2D2D2")
        spine.set_linewidth(0.7)

    header_color = INCOME_COLORS[group]
    ax.axhspan(
        0.922, 1.0,
        color=header_color
    )

    ax.text(
        0.025, 0.962,
        title,
        fontsize=11.2,
        color="white",
        fontweight="bold",
        va="center"
    )

    n_countries = len(frame)
    n_genomes = int(frame["Genome_Count"].sum())

    ax.text(
        0.025, 0.899,
        f"{n_countries} locations | {n_genomes:,} genomes",
        fontsize=9.5,
        fontweight="bold",
        va="top"
    )

    entries = [
        f"{r.Country_Harmonized} ({int(r.Genome_Count):,})"
        for r in frame.itertuples()
    ]

    if not entries:
        ax.text(
            0.03, 0.80, "No classified countries",
            fontsize=9
        )
        return

    # Dynamic 2-column display, all entries retained.
    rows = math.ceil(len(entries) / 2)
    dy = 0.78 / max(rows, 1)

    font_size = (
        10.0 if rows <= 18
        else 9.3 if rows <= 23
        else 8.5 if rows <= 30
        else 7.6
    )

    for i, entry in enumerate(entries):
        col = i // rows
        row = i % rows

        ax.text(
            0.025 + 0.50 * col,
            0.845 - row * dy,
            entry,
            fontsize=font_size,
            va="top",
            clip_on=True
        )


def figure_1a(output):
    print("\nCreating Figure 1A...")

    fig = plt.figure(
        figsize=(18.5, 18.0),
        facecolor="white"
    )

    fig.text(
        0.035, 0.978,
        "Global geographical distribution of "
        "Enterobacter cloacae complex genomes",
        fontsize=19,
        fontweight="bold",
        va="top"
    )

    fig.text(
        0.035, 0.951,
        f"{N:,} genomes | Collection period: 2000–2026",
        fontsize=12.5,
        color="#555555",
        va="top"
    )

    # MAP
    ax = fig.add_axes(
        [0.027, 0.495, 0.946, 0.435]
    )

    world.plot(
        ax=ax,
        facecolor="#EFEFEF",
        edgecolor="#BBBBBB",
        linewidth=0.38
    )

    for group in INCOME_ORDER:
        block = world[
            world["Income_Group"] == group
        ]

        if len(block):
            block.plot(
                ax=ax,
                facecolor=INCOME_COLORS[group],
                edgecolor="#777777",
                linewidth=0.38,
                alpha=0.68
            )

    mapped = countries.dropna(
        subset=["Longitude", "Latitude"]
    ).copy()

    for group in INCOME_ORDER:
        block = mapped[
            mapped["Income_Group"] == group
        ]

        if len(block):
            ax.scatter(
                block["Longitude"],
                block["Latitude"],
                s=bubble_size(block["Genome_Count"]),
                c=INCOME_COLORS[group],
                edgecolors="white",
                linewidths=0.75,
                alpha=0.95,
                zorder=10
            )

    # Top-country names and counts
    for row in countries.head(10).itertuples():
        if not np.isfinite(row.Longitude):
            continue
        if not np.isfinite(row.Latitude):
            continue

        ax.annotate(
            f"{row.Country_Harmonized}\n"
            f"{row.Genome_Count:,}",
            xy=(row.Longitude, row.Latitude),
            xytext=(0, 12),
            textcoords="offset points",
            ha="center",
            fontsize=9.5,
            fontweight="bold",
            bbox={
                "facecolor": "white",
                "edgecolor": "#CFCFCF",
                "alpha": 0.95,
                "boxstyle": "round,pad=0.2"
            },
            zorder=20
        )

    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 85)
    ax.set_aspect("auto")
    ax.axis("off")

    # Income color legend
    handles = [
        Patch(
            facecolor=INCOME_COLORS[g],
            label=INCOME_NAMES[g].title()
        )
        for g in INCOME_ORDER
    ]

    fig.legend(
        handles=handles,
        loc="center",
        bbox_to_anchor=(0.5, 0.489),
        ncol=4,
        frameon=False,
        fontsize=10.5
    )

    # Income method and threshold
    fig.text(
        0.5, 0.468,
        "World Bank FY2027 income thresholds "
        "(2025 GNI per capita, Atlas method)",
        ha="center",
        fontsize=11.3,
        fontweight="bold"
    )

    fig.text(
        0.5, 0.447,
        "High: > US$14,375  |  "
        "Upper-middle: US$4,636–14,375  |  "
        "Low + lower-middle: ≤ US$4,635",
        ha="center",
        fontsize=10.3
    )

    # All income categories written below the map
    for group, left in [
        ("H", 0.025),
        ("M", 0.350),
        ("L+LM", 0.675)
    ]:
        box = fig.add_axes(
            [left, 0.133, 0.300, 0.290]
        )

        sub = countries[
            countries["Income_Group"] == group
        ].sort_values(
            ["Genome_Count", "Country_Harmonized"],
            ascending=[False, True]
        )

        draw_country_list(
            box, sub, group, INCOME_NAMES[group]
        )

    # Unclassified income locations are explicitly listed
    if not unclassified.empty:
        unresolved_text = "; ".join(
            f"{r.Country_Harmonized} ({r.Genome_Count:,})"
            for r in unclassified.itertuples()
        )
        fig.text(
            0.035, 0.116,
            "Unclassified income: " + unresolved_text,
            fontsize=9.5,
            color="#555555",
            va="top"
        )

    # Bubble legend
    fig.text(
        0.035, 0.089,
        "Genome count (bubble area)",
        fontsize=11,
        fontweight="bold"
    )

    bubble_ax = fig.add_axes(
        [0.31, 0.052, 0.63, 0.075]
    )
    bubble_ax.set_xlim(0, 5.7)
    bubble_ax.set_ylim(0, 2)
    bubble_ax.axis("off")

    for i, count in enumerate([
        10, 100, 500, 1000, 2500
    ]):
        x = 0.50 + i * 1.12
        bubble_ax.scatter(
            x, 1.10,
            s=bubble_size(count),
            color="#737373",
            edgecolor="white",
            linewidth=0.7
        )
        bubble_ax.text(
            x, 0.13, f"{count:,}",
            ha="center", fontsize=10
        )

    fig.text(
        0.035, 0.026,
        "Genome counts reflect publicly available "
        "assemblies, not population prevalence. "
        "Country income and geographic mapping status "
        "are documented in the Excel source data.",
        fontsize=9.2,
        color="#555555"
    )

    save_600(fig, output / FILE_A)

    with pd.ExcelWriter(
        output / DATA_A,
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

        data.groupby([
            "Country_Harmonized", "Source"
        ]).size().reset_index(
            name="Genome_Count"
        ).to_excel(
            writer,
            sheet_name="COUNTRY_SOURCE_COUNTS",
            index=False
        )

        data[[
            "Assembly_Accession",
            "Year",
            "Country_Harmonized",
            "Normalized_Species",
            "Source",
            "Income_Group",
            "Income_Source",
            "WB_Code"
        ]].to_excel(
            writer,
            sheet_name="GENOME_SOURCE_DATA",
            index=False
        )

        pd.DataFrame([
            ("Final genomes", N),
            ("Year window", "2000-2026"),
            ("World Bank API", WB_URL),
            ("Income threshold basis", "FY2027"),
            ("Map boundaries", MAP_URL),
            ("Font", FONT),
            ("Resolution DPI", 600),
            ("Build timestamp", datetime.now().isoformat())
        ], columns=[
            "Parameter", "Value"
        ]).to_excel(
            writer,
            sheet_name="FIGURE_METADATA",
            index=False
        )


# ============================================================
# 6. FIGURE 1B - THINNER ANNUAL BARS
# ============================================================

def figure_1b(output):
    print("\nCreating Figure 1B...")

    years = np.arange(2000, 2027)

    annual = (
        data.groupby("Year")
        .size()
        .reindex(years, fill_value=0)
    )

    species_order = (
        data["Normalized_Species"]
        .value_counts()
        .index.tolist()
    )

    source_tab = yearly_counts(
        data, "Source", SOURCE_ORDER
    )

    species_tab = yearly_counts(
        data, "Normalized_Species", species_order
    )

    income_tab = yearly_counts(
        data, "Income_Group", INCOME_ORDER
    )

    assert int(annual.sum()) == N
    assert int(source_tab.to_numpy().sum()) == N
    assert int(species_tab.to_numpy().sum()) == N
    assert int(income_tab.to_numpy().sum()) == N

    species_palette = {
        name: SPECIES_COLORS[i % len(SPECIES_COLORS)]
        for i, name in enumerate(species_order)
    }

    species_labels = {
        name: name.replace(
            "Enterobacter ", "E. "
        ).replace(
            "E. cloacae complex - unresolved",
            "Unresolved ECC"
        )
        for name in species_order
    }

    fig, axes = plt.subplots(
        4, 1,
        figsize=(15, 15.2),
        sharex=True,
        gridspec_kw={
            "height_ratios": [1.0, 1.28, 1.28, 1.28],
            "hspace": 0.42
        }
    )

    fig.suptitle(
        "Temporal distribution of Enterobacter "
        "cloacae complex genomes",
        fontsize=18,
        fontweight="bold",
        y=0.983
    )

    fig.text(
        0.5, 0.956,
        f"{N:,} genomes collected between 2000 and 2026",
        ha="center",
        fontsize=11.5,
        color="#555555"
    )

    # Overall: thin bars
    axes[0].bar(
        years,
        annual.values,
        width=BAR_WIDTH,
        color="#444444",
        edgecolor="none"
    )

    axes[0].set_title(
        "Annual genome collection",
        loc="left",
        fontweight="bold",
        pad=10
    )

    def stack(ax, frame, color_map, label_map,
              columns=4):

        bottom = np.zeros(len(years))

        for key in frame.columns:
            values = frame[key].to_numpy()

            ax.bar(
                years,
                values,
                width=BAR_WIDTH,
                bottom=bottom,
                color=color_map[key],
                edgecolor="none",
                label=label_map.get(key, key)
            )
            bottom += values

        ax.legend(
            loc="upper left",
            bbox_to_anchor=(0.0, 1.015),
            ncol=columns,
            frameon=False,
            fontsize=9.2,
            handlelength=1.25,
            columnspacing=1.2
        )

    # One Health
    stack(
        axes[1],
        source_tab,
        SOURCE_COLORS,
        {},
        columns=5
    )
    axes[1].set_title(
        "One Health sources",
        loc="left",
        fontweight="bold",
        pad=10
    )

    # Normalized species
    stack(
        axes[2],
        species_tab,
        species_palette,
        species_labels,
        columns=4
    )
    axes[2].set_title(
        "Normalized Enterobacter species",
        loc="left",
        fontweight="bold",
        pad=10
    )

    # Income groups
    stack(
        axes[3],
        income_tab,
        INCOME_COLORS,
        {
            "H": "High",
            "M": "Upper-middle",
            "L+LM": "Low + lower-middle",
            "Unclassified": "Unclassified"
        },
        columns=4
    )
    axes[3].set_title(
        "Country-income groups",
        loc="left",
        fontweight="bold",
        pad=10
    )

    # Publication-style axes
    for ax in axes:
        format_axes(ax)
        ax.set_xlim(1999.4, 2026.6)
        ax.set_ylim(bottom=0)

    axes[-1].set_xticks(
        [2000, 2005, 2010, 2015, 2020, 2026]
    )

    axes[-1].set_xlabel(
        "Collection year",
        fontsize=11.5
    )

    fig.subplots_adjust(
        left=0.105,
        right=0.965,
        top=0.923,
        bottom=0.073
    )

    fig.text(
        0.5, 0.025,
        "Bar height indicates the number of available "
        "genomes and does not represent population "
        "incidence or prevalence.",
        ha="center",
        fontsize=9.4,
        color="#555555"
    )

    save_600(fig, output / FILE_B)

    with pd.ExcelWriter(
        output / DATA_B,
        engine="openpyxl"
    ) as writer:

        annual.rename(
            "Genome_Count"
        ).rename_axis(
            "Year"
        ).reset_index().to_excel(
            writer,
            sheet_name="ANNUAL_COUNTS",
            index=False
        )

        for sheet_name, frame, col in [
            ("YEAR_BY_SOURCE", source_tab, "Source"),
            ("YEAR_BY_SPECIES", species_tab, "Species"),
            ("YEAR_BY_INCOME", income_tab, "Income_Group")
        ]:
            frame.rename_axis(
                "Year"
            ).reset_index().melt(
                id_vars="Year",
                var_name=col,
                value_name="Genome_Count"
            ).to_excel(
                writer,
                sheet_name=sheet_name,
                index=False
            )

        data[[
            "Assembly_Accession",
            "Year",
            "Country_Harmonized",
            "Source",
            "Original_Species",
            "Normalized_Species",
            "Income_Group"
        ]].to_excel(
            writer,
            sheet_name="GENOME_SOURCE_DATA",
            index=False
        )

        pd.DataFrame([
            ("Final genomes", N),
            ("Time window", "2000-2026"),
            ("Bar width", BAR_WIDTH),
            ("Number of species groups", len(species_order)),
            ("DPI", 600),
            ("Font", FONT),
            ("Created", datetime.now().isoformat())
        ], columns=[
            "Parameter", "Value"
        ]).to_excel(
            writer,
            sheet_name="FIGURE_METADATA",
            index=False
        )


# ============================================================
# 7. STAGE AND VALIDATE BOTH NEW FIGURES
# ============================================================

print("\nCreating publication figures in staging folders...")

with tempfile.TemporaryDirectory(
    prefix="ECC_NATURE_FINAL_",
    dir=str(PROJECT)
) as tmp:

    stage = Path(tmp)
    a = stage / "Figure_1A"
    b = stage / "Figure_1B"
    a.mkdir()
    b.mkdir()

    figure_1a(a)
    figure_1b(b)

    expected_a = {FILE_A, DATA_A}
    expected_b = {FILE_B, DATA_B}

    if {p.name for p in a.iterdir()} != expected_a:
        raise RuntimeError("Figure 1A output validation failed")

    if {p.name for p in b.iterdir()} != expected_b:
        raise RuntimeError("Figure 1B output validation failed")

    a_check = pd.read_excel(
        a / DATA_A,
        sheet_name="COUNTRY_INCOME_COUNTS"
    )
    b_check = pd.read_excel(
        b / DATA_B,
        sheet_name="ANNUAL_COUNTS"
    )

    if int(a_check["Genome_Count"].sum()) != N:
        raise RuntimeError("Figure 1A total count mismatch")

    if int(b_check["Genome_Count"].sum()) != N:
        raise RuntimeError("Figure 1B total count mismatch")

    for p in [a / FILE_A, b / FILE_B]:
        if p.stat().st_size < 100000:
            raise RuntimeError(
                "Unexpectedly small image: " + str(p)
            )

    # Only now replace the previous figure folders
    for target, new in [
        (FOLDER_A, a),
        (FOLDER_B, b)
    ]:
        if target.exists():
            shutil.rmtree(target)

        shutil.copytree(new, target)

# ============================================================
# 8. FINAL VERIFICATION AND SUMMARY
# ============================================================

assert len(list(SUPPLEMENT.iterdir())) == 1
assert sha256(MASTER) == sha256(supp_master)

print("\n" + "=" * 90)
print("FINAL ECC FIGURES COMPLETED")
print("=" * 90)

print("\nFINAL MASTER:", MASTER)
print("SUPPLEMENTARY MASTER:", supp_master)
print("GENOMES:", N)
print("YEARS:", START_YEAR, "-", END_YEAR)

print("\nINCOME SUMMARY")
print(income_summary.to_string(index=False))

print("\nSPECIES SUMMARY")
print(
    data["Normalized_Species"]
    .value_counts()
    .to_string()
)

print("\nFIGURE 1A")
for path in sorted(FOLDER_A.iterdir()):
    print(" -", path.name)

print("\nFIGURE 1B")
for path in sorted(FOLDER_B.iterdir()):
    print(" -", path.name)

print("\nSUPPLEMENTARY FOLDER")
for path in sorted(SUPPLEMENT.iterdir()):
    print(" -", path.name)

print("\nOnly one Excel master in supplementary folder.")
print("Figure 1A: income-classified country list included.")
print("Figure 1B: bar width =", BAR_WIDTH)
print("Both figures: 600 DPI JPEG and Excel source data.")
print("=" * 90)
