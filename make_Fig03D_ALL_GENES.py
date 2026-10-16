
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Rectangle
from matplotlib import font_manager
from openpyxl.utils import get_column_letter

# ============================================================
# FIGURE 3D - COMPLETE ECC AMR GENE HEATMAP
#
# ALL SEVEN RESOLVED SPECIES
# ALL GENES FROM AMR_Genes
# ALL FIVE ONE HEALTH SOURCES
# THREE FIXED INCOME GROUPS
#
# ROWS: individual AMR genes, ordered by class
#
# COLUMNS:
# species -> source -> income
#
# COLOR = genome percentage carrying gene
#
# LEFT:
# gene name and distinct AMR-class color strip
#
# Class assignment:
# 1. User-verified GENE_CLASS_MAPPING.xlsx if present
# 2. Provisional recognizable gene-family assignment
# 3. Otherwise "Unassigned - verify"
#
# Provisional assignment is not definitive annotation.
#
# OUTPUT:
# ONE 600-DPI JPEG
# ONE Excel source workbook
#
# NO statistical testing in this gene heatmap.
# The earlier class-level Mann-Whitney results
# must not be reused as gene-level significance.
#
# DOES NOT CHANGE THE MASTER OR PREVIOUS FIGURES.
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

OUT = ROOT / "Fig_03D_ALL_SPECIES_GENE_HEATMAP"
OUT.mkdir(parents=True, exist_ok=True)

MAPPING_FILE = OUT / "GENE_CLASS_MAPPING.xlsx"

JPG = OUT / "Fig_03D_ALL_SPECIES_ALL_AMR_GENES.jpg"
XLSX = OUT / "Fig_03D_ALL_SPECIES_ALL_AMR_GENES_Data.xlsx"

SPECIES = [
    "Enterobacter hormaechei",
    "Enterobacter cloacae",
    "Enterobacter roggenkampii",
    "Enterobacter asburiae",
    "Enterobacter kobei",
    "Enterobacter ludwigii",
    "Enterobacter bugandensis"
]

SOURCES = [
    "Human",
    "Animal",
    "Environment",
    "Food",
    "Plant-associated"
]

INCOMES = ["H", "M", "L+LM"]

INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82"
}

SOURCE_COLORS = {
    "Human": "#E69F00",
    "Animal": "#56B4E9",
    "Environment": "#009E73",
    "Food": "#B07AA1",
    "Plant-associated": "#A6761D"
}

SPECIES_COLORS = {
    "Enterobacter hormaechei": "#0072B2",
    "Enterobacter cloacae": "#D55E00",
    "Enterobacter roggenkampii": "#009E73",
    "Enterobacter asburiae": "#CC79A7",
    "Enterobacter kobei": "#7851A9",
    "Enterobacter ludwigii": "#A6761D",
    "Enterobacter bugandensis": "#56B4E9"
}

CLASS_ORDER = [
    "Beta-lactam",
    "Aminoglycoside",
    "Quinolone",
    "Tetracycline",
    "Sulfonamide",
    "Trimethoprim",
    "Phenicols",
    "Macrolide",
    "Fosfomycin",
    "Colistin",
    "Rifamycin",
    "Bleomycin",
    "Streptothricin",
    "Lincosamide",
    "Other verified class",
    "Unassigned - verify"
]

CLASS_COLORS = {
    "Beta-lactam": "#D55E00",
    "Aminoglycoside": "#0072B2",
    "Quinolone": "#009E73",
    "Tetracycline": "#E69F00",
    "Sulfonamide": "#CC79A7",
    "Trimethoprim": "#7851A9",
    "Phenicols": "#A6761D",
    "Macrolide": "#56B4E9",
    "Fosfomycin": "#D49A6A",
    "Colistin": "#B2182B",
    "Rifamycin": "#4878A8",
    "Bleomycin": "#7A9E48",
    "Streptothricin": "#C17C74",
    "Lincosamide": "#AB6FBE",
    "Other verified class": "#566573",
    "Unassigned - verify": "#B8B8B8"
}

# ============================================================
# FONT
# ============================================================

installed_fonts = {
    f.name for f in font_manager.fontManager.ttflist
}

FONT = next(
    (
        f for f in [
            "Arial",
            "Helvetica",
            "Liberation Sans",
            "DejaVu Sans"
        ] if f in installed_fonts
    ),
    "DejaVu Sans"
)

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 10,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})

# ============================================================
# UTILITY FUNCTIONS
# ============================================================

def parse_gene_set(value):
    if pd.isna(value):
        return set()

    return {
        g.strip()
        for g in str(value).split(";")
        if g.strip()
    }


def provisional_class(gene):
    """
    Transparent family-name classification.
    These are provisional classifications,
    NOT verified AMRFinder gene annotations.

    Ambiguous names remain unassigned.
    """

    g = gene.strip().lower()

    if g.startswith("bla"):
        return "Beta-lactam"

    if re.match(r"^(aac|aph|ant)\(", g):
        return "Aminoglycoside"

    if re.match(r"^(aac|aph|ant)[-_]", g):
        return "Aminoglycoside"

    if re.match(r"^(qnr|oqx)", g):
        return "Quinolone"

    if re.match(r"^tet(\(|[a-z])", g):
        return "Tetracycline"

    if re.match(r"^sul[123]($|[-_])", g):
        return "Sulfonamide"

    if re.match(r"^dfr", g):
        return "Trimethoprim"

    if re.match(r"^(cat|cml|flo)", g):
        return "Phenicols"

    if re.match(r"^(erm|mph|mef|msr)", g):
        return "Macrolide"

    if re.match(r"^fos", g):
        return "Fosfomycin"

    if re.match(r"^mcr[-_]", g):
        return "Colistin"

    if re.match(r"^arr", g):
        return "Rifamycin"

    if re.match(r"^ble", g):
        return "Bleomycin"

    if re.match(r"^sat", g):
        return "Streptothricin"

    if re.match(r"^lnu", g):
        return "Lincosamide"

    return "Unassigned - verify"


def class_sort_key(gene, gene_class, prevalence):
    rank = CLASS_ORDER.index(gene_class[gene])
    return (
        rank,
        -prevalence[gene],
        gene.lower()
    )


# ============================================================
# INPUTS
# ============================================================

print("=" * 100)
print("FIGURE 3D - ALL ECC SPECIES, ALL AMR GENES")
print("=" * 100)

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

required = [
    "Assembly_Accession",
    "Normalized_Species",
    "Final_Source_Category",
    "Country_Harmonized",
    "AMR_Genes",
    "AMRFinder_Status"
]

missing = [c for c in required if c not in df.columns]

if missing:
    raise ValueError(
        "Missing master columns: " + str(missing)
    )

if len(df) != 9157:
    raise ValueError(
        f"Expected 9157 genomes, got {len(df)}"
    )

if (
    df["Assembly_Accession"].isna().any()
    or df["Assembly_Accession"].duplicated().any()
):
    raise ValueError("Duplicate or missing accessions")

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

status = (
    df["AMRFinder_Status"]
    .fillna("").astype(str).str.strip().str.upper()
)

if not status.eq("COMPLETE").all():
    raise ValueError("Incomplete AMRFinder records")

if (
    not df["Source"].isin(SOURCES).all()
    or df["Country"].eq("").any()
):
    raise ValueError("Unexpected source or country")

# ============================================================
# FROZEN INCOME LOOKUP
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
        raise ValueError("Missing income field: " + c)

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
    raise ValueError("Duplicate income country")

check = pd.concat(
    [
        df["Country"].value_counts().rename("Master"),
        income.set_index("Country")[
            "Genome_Count"
        ].rename("Reference")
    ],
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
        "Income database does not match master"
    )

df = df.merge(
    income[["Country", "Income_Group"]],
    on="Country",
    how="left",
    validate="many_to_one"
)

# All seven resolved species with classified income.
data = df[
    df["Species"].isin(SPECIES)
    & df["Income_Group"].isin(INCOMES)
].copy()

if len(data) != 8867:
    raise ValueError(
        f"Expected 8867 included genomes; "
        f"found {len(data)}"
    )

print("Input genomes:", len(df))
print("Analyzed genomes:", len(data))

# ============================================================
# PARSE EVERY AMR GENE
# ============================================================

print("\nParsing all AMR genes...")

data["Gene_Set"] = (
    data["AMR_Genes"].apply(parse_gene_set)
)

gene_occurrences = Counter()

for genes in data["Gene_Set"]:
    gene_occurrences.update(genes)

ALL_GENES = sorted(gene_occurrences)

if not ALL_GENES:
    raise ValueError("No AMR genes found")

print("Unique resistance gene labels:", len(ALL_GENES))

# ============================================================
# BUILD / READ EDITABLE GENE-CLASS MAP
# ============================================================

default_mapping = pd.DataFrame({
    "Gene": ALL_GENES,
    "Assigned_Class": [
        provisional_class(g)
        for g in ALL_GENES
    ],
    "Assignment_Source": [
        (
            "Provisional gene-family rule"
            if provisional_class(g)
            != "Unassigned - verify"
            else "Unassigned"
        )
        for g in ALL_GENES
    ],
    "Occurrences": [
        gene_occurrences[g]
        for g in ALL_GENES
    ]
})

if MAPPING_FILE.is_file():

    print(
        "Using existing editable mapping:",
        MAPPING_FILE
    )

    user_map = pd.read_excel(
        MAPPING_FILE,
        sheet_name="GENE_CLASS_MAPPING",
        engine="openpyxl"
    )

    if not {
        "Gene",
        "Assigned_Class"
    }.issubset(user_map.columns):
        raise ValueError(
            "Mapping requires Gene and Assigned_Class"
        )

    user_map["Gene"] = (
        user_map["Gene"].astype(str).str.strip()
    )

    if user_map["Gene"].duplicated().any():
        raise ValueError(
            "Duplicate gene entries in mapping"
        )

    existing = user_map.set_index("Gene")

    for i, row in default_mapping.iterrows():

        gene = row["Gene"]

        if gene not in existing.index:
            continue

        value = existing.at[
            gene, "Assigned_Class"
        ]

        if pd.isna(value):
            continue

        value = str(value).strip()

        if not value:
            continue

        if value not in CLASS_ORDER:
            raise ValueError(
                f"Unknown class '{value}' "
                f"for gene {gene}"
            )

        default_mapping.at[
            i, "Assigned_Class"
        ] = value

        default_mapping.at[
            i, "Assignment_Source"
        ] = "User mapping file"

else:

    print(
        "\nCreating editable gene-to-class mapping."
    )

    with pd.ExcelWriter(
        MAPPING_FILE,
        engine="openpyxl"
    ) as writer:

        default_mapping.to_excel(
            writer,
            sheet_name="GENE_CLASS_MAPPING",
            index=False
        )

        pd.DataFrame({
            "Valid_Class": CLASS_ORDER
        }).to_excel(
            writer,
            sheet_name="VALID_CLASSES",
            index=False
        )

gene_class = dict(zip(
    default_mapping["Gene"],
    default_mapping["Assigned_Class"]
))

# ============================================================
# CREATE 105 SPECIES / SOURCE / INCOME GROUPS
# ============================================================

groups = [
    (sp, source, inc)
    for sp in SPECIES
    for source in SOURCES
    for inc in INCOMES
]

counts = defaultdict(Counter)
group_n = Counter()

for row in data.itertuples():

    group_key = (
        row.Species,
        row.Source,
        row.Income_Group
    )

    group_n[group_key] += 1

    for gene in row.Gene_Set:
        counts[group_key][gene] += 1

if sum(group_n.values()) != len(data):
    raise RuntimeError("Group counts do not reconcile")

# ============================================================
# MATRIX OF GENE PREVALENCE
# ============================================================

genes_ordered = sorted(
    ALL_GENES,
    key=lambda gene: class_sort_key(
        gene, gene_class, gene_occurrences
    )
)

NROWS = len(genes_ordered)
NCOLS = len(groups)

matrix = np.full(
    (NROWS, NCOLS),
    np.nan,
    dtype=np.float32
)

carriage = np.zeros(
    (NROWS, NCOLS),
    dtype=np.int32
)

group_total = np.zeros(
    NCOLS, dtype=np.int32
)

for col_idx, group in enumerate(groups):

    n = group_n[group]
    group_total[col_idx] = n

    if n == 0:
        continue

    for row_idx, gene in enumerate(genes_ordered):

        positive = counts[group][gene]

        carriage[row_idx, col_idx] = positive

        matrix[row_idx, col_idx] = (
            100.0 * positive / n
        )

# Every gene must occur in at least one included genome.
if not (
    np.nanmax(matrix, axis=1) > 0
).all():
    raise RuntimeError(
        "Some genes disappeared from the matrix"
    )

# ============================================================
# FIGURE DIMENSIONS
#
# Canvas adapts to the number of genes.
# All genes retained.
# Maximum figure height = 90 inches,
# corresponding to 54000 pixels at 600 DPI.
# ============================================================

fig_height = min(
    90.0,
    max(17.0, 6.0 + NROWS * 0.20)
)

fig_width = 30.0

fontsize_gene = (
    9.0 if NROWS <= 200
    else 7.5 if NROWS <= 350
    else 6.2
)

fontsize_header = 7.5

print(
    f"\nHeatmap dimensions: "
    f"{NROWS} genes x {NCOLS} columns"
)

print(
    f"Figure size: "
    f"{fig_width:.1f} x {fig_height:.1f} inches"
)

# ============================================================
# COLOR MAP
# ============================================================

cmap = LinearSegmentedColormap.from_list(
    "AMR_gene_prevalence",
    [
        "#F7FBFF",
        "#D9EAF4",
        "#91C8D7",
        "#45A4AF",
        "#19738F",
        "#16365F"
    ]
)

cmap.set_bad("#DEDEDE")

# ============================================================
# DRAW FIGURE
# ============================================================

print("\nRendering complete gene-level heatmap...")

fig = plt.figure(
    figsize=(fig_width, fig_height),
    facecolor="white"
)

# Separate annotation-strip axes avoid overlap.
#
# Gene text sits at left of the class strip.
# Strip does not sit under gene labels.

ax_gene = fig.add_axes([
    0.012, 0.085, 0.112, 0.835
])

ax_class = fig.add_axes([
    0.128, 0.085, 0.012, 0.835
])

ax_heat = fig.add_axes([
    0.145, 0.085, 0.79, 0.835
])

ax_gene.set_xlim(0, 1)
ax_gene.set_ylim(NROWS - .5, -.5)
ax_gene.axis("off")

ax_class.set_xlim(0, 1)
ax_class.set_ylim(NROWS - .5, -.5)
ax_class.axis("off")

masked = np.ma.masked_invalid(matrix)

im = ax_heat.imshow(
    masked,
    cmap=cmap,
    vmin=0,
    vmax=100,
    aspect="auto",
    interpolation="nearest"
)

ax_heat.set_xlim(-.5, NCOLS - .5)
ax_heat.set_ylim(NROWS - .5, -.5)

ax_heat.set_yticks([])
ax_heat.set_xticks([])

for spine in ax_heat.spines.values():
    spine.set_visible(False)

# ============================================================
# LEFT GENE LABELS AND CLASS COLOR STRIP
# ============================================================

previous_class = None

for i, gene in enumerate(genes_ordered):

    gene_category = gene_class[gene]

    ax_gene.text(
        0.98,
        i,
        gene,
        ha="right",
        va="center",
        fontsize=fontsize_gene,
        color="#222222"
    )

    ax_class.add_patch(
        Rectangle(
            (0, i - 0.5),
            1,
            1,
            facecolor=CLASS_COLORS[gene_category],
            edgecolor="none"
        )
    )

    if (
        previous_class is not None
        and previous_class != gene_category
    ):

        for axis in [
            ax_gene,
            ax_class
        ]:
            axis.axhline(
                i - .5,
                color="#777777",
                linewidth=.45
            )

        ax_heat.axhline(
            i - .5,
            color="#555555",
            linewidth=.9
        )

    previous_class = gene_category

# ============================================================
# TOP COLUMN ANNOTATION STRIPS
# ============================================================

# Strips use heatmap coordinates, preserving
# exact column alignment.

species_top = -4.15
source_top = -3.05
income_top = -1.95

for j, (species, source, inc) in enumerate(groups):

    for y0, color in [
        (species_top, SPECIES_COLORS[species]),
        (source_top, SOURCE_COLORS[source]),
        (income_top, INCOME_COLORS[inc])
    ]:

        ax_heat.add_patch(
            Rectangle(
                (j - .5, y0),
                1,
                .70,
                facecolor=color,
                edgecolor="white",
                linewidth=.15,
                clip_on=False
            )
        )

    # Income letters at top.
    ax_heat.text(
        j,
        -4.75,
        inc,
        rotation=90,
        ha="center",
        va="bottom",
        fontsize=fontsize_header,
        color="#333333",
        clip_on=False
    )

# Source separators every 3 columns.
for j in range(3, NCOLS, 3):
    ax_heat.axvline(
        j - .5,
        color="white",
        linewidth=.75,
        zorder=5
    )

# Species separators every 15 columns.
for j in range(15, NCOLS, 15):
    ax_heat.axvline(
        j - .5,
        color="#262626",
        linewidth=1.6,
        zorder=6
    )

# Species labels, centred above each block.
for species_idx, species in enumerate(SPECIES):

    center = species_idx * 15 + 7

    ax_heat.text(
        center,
        -7.1,
        species.replace("Enterobacter ", "E. "),
        ha="center",
        va="bottom",
        fontsize=12.5,
        fontstyle="italic",
        fontweight="bold",
        clip_on=False
    )

# One Health compartment labels are placed
# BELOW the heatmap, outside the plotting area.
for species_idx, species in enumerate(SPECIES):

    base = species_idx * 15

    for source_idx, source in enumerate(SOURCES):

        center = base + source_idx * 3 + 1

        ax_heat.text(
            center,
            NROWS + 1.2,
            source,
            rotation=90,
            ha="center",
            va="top",
            fontsize=8,
            clip_on=False
        )

fig.text(
    0.50,
    0.985,
    "Global ECC resistance-gene repertoire",
    ha="center",
    fontsize=21,
    fontweight="bold"
)

fig.text(
    0.50,
    0.970,
    "Seven species | Five One Health sources | "
    "Three country-income groups | All observed AMR genes",
    ha="center",
    fontsize=12,
    color="#555555"
)

# ============================================================
# CLASS COLOR LEGEND
# ============================================================

present_classes = [
    cls for cls in CLASS_ORDER
    if cls in set(gene_class.values())
]

legend_handles = [
    Patch(
        facecolor=CLASS_COLORS[cls],
        edgecolor="none",
        label=cls
    )
    for cls in present_classes
]

fig.legend(
    handles=legend_handles,
    loc="lower center",
    bbox_to_anchor=(0.52, 0.027),
    ncol=min(6, max(2, len(present_classes))),
    frameon=False,
    fontsize=9
)

# ============================================================
# HEATMAP COLORBAR
# ============================================================

cax = fig.add_axes([
    0.945,
    0.38,
    0.011,
    0.24
])

cb = fig.colorbar(
    im,
    cax=cax
)

cb.set_label(
    "Genomes carrying AMR gene (%)",
    fontsize=11
)

cb.set_ticks([
    0, 20, 40, 60, 80, 100
])

fig.text(
    0.50,
    0.010,
    "Grey = no genomes in the group. "
    "Gene-class assignments based on provisional "
    "gene-family rules require verification.",
    ha="center",
    fontsize=9,
    color="#555555"
)

# ============================================================
# SAVE 600 DPI JPEG
# ============================================================

fig.savefig(
    JPG,
    dpi=600,
    format="jpeg",
    facecolor="white",
    pil_kwargs={
        "quality": 92,
        "subsampling": 0,
        "optimize": True
    }
)

plt.close(fig)

# ============================================================
# EXCEL SOURCE TABLES
# ============================================================

column_names = [
    (
        species.replace("Enterobacter ", "E."),
        source,
        inc
    )
    for species, source, inc in groups
]

prevalence_wide = pd.DataFrame(
    matrix,
    columns=pd.MultiIndex.from_tuples(
        column_names,
        names=["Species", "Source", "Income"]
    )
)

prevalence_wide.insert(
    0,
    ("Metadata", "AMR_Class", ""),
    [gene_class[g] for g in genes_ordered]
)

prevalence_wide.insert(
    0,
    ("Metadata", "Gene", ""),
    genes_ordered
)

# Flatten for Excel compatibility.
prevalence_wide.columns = [
    " | ".join(
        str(x) for x in col if str(x)
    )
    for col in prevalence_wide.columns
]

carrier_wide = pd.DataFrame(
    carriage,
    columns=[
        f"{species.replace('Enterobacter ', 'E.')} "
        f"| {source} | {inc}"
        for species, source, inc in groups
    ]
)

carrier_wide.insert(
    0, "Gene", genes_ordered
)

carrier_wide.insert(
    1,
    "AMR_Class",
    [gene_class[g] for g in genes_ordered]
)

group_df = pd.DataFrame([
    {
        "Species": sp,
        "Source": source,
        "Income": inc,
        "Genomes": group_n[(sp, source, inc)]
    }
    for sp, source, inc in groups
])

gene_df = pd.DataFrame([
    {
        "Gene": gene,
        "AMR_Class": gene_class[gene],
        "Total_Carrier_Genomes": gene_occurrences[gene],
        "Assignment_Source": default_mapping.loc[
            default_mapping["Gene"] == gene,
            "Assignment_Source"
        ].iloc[0]
    }
    for gene in genes_ordered
])

methods = pd.DataFrame([
    (
        "Input",
        "Final 9157 ECC master (2000-2026)"
    ),
    (
        "Analyzed",
        len(data)
    ),
    (
        "Species",
        len(SPECIES)
    ),
    (
        "One Health sources",
        len(SOURCES)
    ),
    (
        "Income groups",
        len(INCOMES)
    ),
    (
        "Unique AMR genes",
        len(ALL_GENES)
    ),
    (
        "Heatmap measurement",
        "100 x genomes containing gene / "
        "all genomes in species-source-income group"
    ),
    (
        "Statistical significance",
        "Not tested in this heatmap"
    ),
    (
        "Class assignment",
        "Verified mappings when supplied; "
        "otherwise provisional name-family assignment "
        "or unassigned"
    ),
    (
        "Unclassified class treatment",
        "Retained; not omitted"
    ),
    (
        "No genomes",
        "NaN / grey, not zero"
    ),
    (
        "Income colors",
        str(INCOME_COLORS)
    ),
    (
        "Created",
        datetime.now().isoformat()
    )
], columns=["Parameter", "Description"])

with pd.ExcelWriter(
    XLSX,
    engine="openpyxl"
) as writer:

    prevalence_wide.to_excel(
        writer,
        sheet_name="GENE_PREVALENCE_MATRIX",
        index=False
    )

    carrier_wide.to_excel(
        writer,
        sheet_name="GENE_CARRIER_COUNTS",
        index=False
    )

    gene_df.to_excel(
        writer,
        sheet_name="GENE_CLASS_INVENTORY",
        index=False
    )

    default_mapping.to_excel(
        writer,
        sheet_name="GENE_CLASS_MAPPING",
        index=False
    )

    group_df.to_excel(
        writer,
        sheet_name="GROUP_DENOMINATORS",
        index=False
    )

    data[[
        "Assembly_Accession",
        "Species",
        "Source",
        "Country",
        "Income_Group",
        "AMR_Genes"
    ]].to_excel(
        writer,
        sheet_name="GENOME_SOURCE_DATA",
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

print("\n" + "=" * 100)
print("ALL SPECIES AMR GENE HEATMAP COMPLETE")
print("=" * 100)

print("Species:", len(SPECIES))
print("One Health compartments:", len(SOURCES))
print("Income groups:", len(INCOMES))
print("Groups:", len(groups))
print("Genes:", len(ALL_GENES))
print("Genomes:", len(data))

print("\nGENES PER ASSIGNED CLASS")

print(
    gene_df["AMR_Class"]
    .value_counts()
    .to_string()
)

print("\nFigure:", JPG)
print("Excel:", XLSX)
print("Editable gene-class mapping:", MAPPING_FILE)

print("\nJPEG: 600 DPI")
print("Existing figure files unchanged")
print("=" * 100)
