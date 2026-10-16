
from pathlib import Path
from collections import Counter
from datetime import datetime
import re
import gc

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle, Patch
from PIL import Image

# ============================================================
# FIGURE 4 — ONE COMPLETE HEATMAP PER ECC SPECIES
#
# 7 SPECIES = EXACTLY 7 FIGURES
#
# ALL DETECTED AMR GENES
# ALL CLASSES DISPLAYED TOGETHER
# CLASS IDENTIFIED BY LEFT COLOR STRIP
#
# HUMAN / ANIMAL / ENVIRONMENT
# EACH DIVIDED INTO H / M / L+LM
#
# PREVALENCE ONLY
# NO P VALUES
# NO PAIRWISE TESTS
# NO PAGINATION
#
# OUTPUT: 7 JPEGs at 600 DPI + 1 EXCEL
#
# Gene-to-class classification is provisional
# unless a verified mapping is supplied.
# ============================================================

ROOT = Path(r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES")

MASTER_FILE = (
    ROOT
    / "00_FINAL_COUNTRY_REPORTED_MASTER"
    / "ECC_FINAL_2000_2026_MASTER.xlsx"
)

INCOME_FILE = (
    ROOT
    / "Fig_01A_Global_Geographical_Distribution"
    / "Fig_01A_Country_Income_Database.xlsx"
)

OUT = ROOT / "Fig_04_AMR_Gene_Prevalence_Heatmaps"
OUT.mkdir(parents=True, exist_ok=True)

EXCEL = OUT / "Fig04_ALL_GENES_PREVALENCE_MASTER.xlsx"

VERIFIED_MAP = OUT / "VERIFIED_GENE_CLASS_MAPPING.xlsx"

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
    "Environment"
]

INCOMES = [
    "H",
    "M",
    "L+LM"
]

# Preserve established income colors.
INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82"
}

# Source-strip colors, distinct from one another.
SOURCE_COLORS = {
    "Human": "#D86B18",
    "Animal": "#1577BD",
    "Environment": "#008C73"
}

# Each AMR class gets its own class-strip color.
# None exactly duplicates source or income colors.
CLASS_COLORS = {
    "Beta-lactam": "#E85D3F",
    "Aminoglycoside": "#27A8DF",
    "Quinolone": "#55B88A",
    "Tetracycline": "#E7BF32",
    "Sulfonamide": "#DF63B2",
    "Trimethoprim": "#9572CE",
    "Phenicols": "#AA804C",
    "Macrolide": "#45C1CD",
    "Fosfomycin": "#ED9752",
    "Colistin": "#B53B54",
    "Rifamycin": "#698EBB",
    "Bleomycin": "#788B95",
    "Lincosamide": "#BB92C6",
    "Streptothricin": "#A0B74F",
    "Other verified class": "#53687C",
    "Unassigned - verify": "#C8CDD2"
}

CLASS_ORDER = list(CLASS_COLORS.keys())

DPI = 600

# One heatmap per species, whatever the gene count.
# A generous vertical canvas keeps row names readable.
# Maximum image height stays below common JPEG limits.
MAX_FIGURE_HEIGHT = 58

# ============================================================
# HEATMAP COLORS
# ============================================================

CMAP = LinearSegmentedColormap.from_list(
    "AMR_prevalence",
    [
        "#F7FBFF",
        "#DCECF9",
        "#B1D6F3",
        "#70B5E8",
        "#348DD3",
        "#1567B0",
        "#123B78"
    ]
)

CMAP.set_bad("#D9DDE2")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})

# ============================================================
# FUNCTIONS
# ============================================================

def parse_genes(value):
    if pd.isna(value):
        return set()

    return {
        g.strip()
        for g in str(value).split(";")
        if g.strip()
    }


def provisional_class(gene):

    g = gene.lower().strip()

    if g.startswith("bla"):
        return "Beta-lactam"

    if re.match(r"^(aac|aph|ant)(\(|-|_)", g):
        return "Aminoglycoside"

    if g.startswith(("qnr", "oqx")):
        return "Quinolone"

    if re.match(r"^tet(\(|[a-z])", g):
        return "Tetracycline"

    if re.match(r"^sul[123]($|[-_])", g):
        return "Sulfonamide"

    if g.startswith("dfr"):
        return "Trimethoprim"

    if g.startswith(("cat", "cml", "flo")):
        return "Phenicols"

    if g.startswith(("erm", "mph", "mef", "msr")):
        return "Macrolide"

    if g.startswith("fos"):
        return "Fosfomycin"

    if re.match(r"^mcr[-_]", g):
        return "Colistin"

    if g.startswith("arr"):
        return "Rifamycin"

    if g.startswith("ble"):
        return "Bleomycin"

    if g.startswith("lnu"):
        return "Lincosamide"

    if g.startswith("sat"):
        return "Streptothricin"

    return "Unassigned - verify"


def species_tag(name):
    return name.replace(
        "Enterobacter ", "E_"
    ).replace(" ", "_")


def save_image(fig, path):

    fig.savefig(
        path,
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

    with Image.open(path) as im:
        im.verify()

    if path.stat().st_size < 10000:
        raise RuntimeError(
            "Invalid image: " + str(path)
        )


# ============================================================
# LOAD DATA
# ============================================================

print("=" * 100)
print("FIGURE 4 — SINGLE AMR PREVALENCE HEATMAP PER SPECIES")
print("=" * 100)

if not MASTER_FILE.is_file():
    raise FileNotFoundError(MASTER_FILE)

if not INCOME_FILE.is_file():
    raise FileNotFoundError(INCOME_FILE)

master = pd.read_excel(
    MASTER_FILE,
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

for col in required:
    if col not in master.columns:
        raise ValueError(
            "Missing master column: " + col
        )

if len(master) != 9157:
    raise ValueError(
        f"Expected 9157 genomes; found {len(master)}"
    )

if (
    master["Assembly_Accession"].isna().any()
    or master["Assembly_Accession"].duplicated().any()
):
    raise ValueError(
        "Invalid genome accession identifiers"
    )

status = (
    master["AMRFinder_Status"]
    .fillna("")
    .astype(str)
    .str.strip()
    .str.upper()
)

if not status.eq("COMPLETE").all():
    raise ValueError(
        "Incomplete AMRFinder records"
    )

master["Species"] = (
    master["Normalized_Species"]
    .fillna("").astype(str).str.strip()
)

master["Source"] = (
    master["Final_Source_Category"]
    .fillna("").astype(str).str.strip()
)

master["Country"] = (
    master["Country_Harmonized"]
    .fillna("").astype(str).str.strip()
)

income = pd.read_excel(
    INCOME_FILE,
    sheet_name="COUNTRY_INCOME_COUNTS",
    engine="openpyxl"
)

for col in [
    "Country_Harmonized",
    "Income_Group"
]:
    if col not in income.columns:
        raise ValueError(
            "Missing income field: " + col
        )

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
    raise ValueError(
        "Duplicate income country"
    )

master = master.merge(
    income[["Country", "Income_Group"]],
    on="Country",
    how="left",
    validate="many_to_one"
)

data = master[
    master["Species"].isin(SPECIES)
    & master["Source"].isin(SOURCES)
    & master["Income_Group"].isin(INCOMES)
].copy()

if data.empty:
    raise ValueError("No genomes after filtering")

data["Gene_Set"] = (
    data["AMR_Genes"].apply(parse_genes)
)

print("Genomes included:", len(data))

# ============================================================
# GENE-TO-CLASS ANNOTATION
# ============================================================

all_genes = sorted(
    set().union(*data["Gene_Set"].tolist())
)

if not all_genes:
    raise ValueError("No AMR genes identified")

mapping = pd.DataFrame({
    "Gene": all_genes,
    "AMR_Class": [
        provisional_class(g) for g in all_genes
    ],
    "Annotation_Status": [
        (
            "Provisional"
            if provisional_class(g) != "Unassigned - verify"
            else "Unassigned"
        )
        for g in all_genes
    ]
})

if VERIFIED_MAP.exists():

    verified = pd.read_excel(
        VERIFIED_MAP,
        engine="openpyxl"
    )

    if not {
        "Gene", "AMR_Class"
    }.issubset(verified.columns):
        raise ValueError(
            "Verified mapping needs Gene and AMR_Class"
        )

    if verified["Gene"].duplicated().any():
        raise ValueError(
            "Duplicate gene in verified mapping"
        )

    verified_lookup = (
        verified.set_index("Gene")["AMR_Class"]
        .dropna()
        .to_dict()
    )

    for i, gene in enumerate(mapping["Gene"]):

        if gene not in verified_lookup:
            continue

        cls = str(
            verified_lookup[gene]
        ).strip()

        if cls not in CLASS_ORDER:
            raise ValueError(
                f"Unknown verified class: {cls}"
            )

        mapping.at[i, "AMR_Class"] = cls
        mapping.at[
            i, "Annotation_Status"
        ] = "Verified mapping"

gene_class = dict(zip(
    mapping["Gene"],
    mapping["AMR_Class"]
))

class_rank = {
    cls: i
    for i, cls in enumerate(CLASS_ORDER)
}

# ============================================================
# BUILD PREVALENCE MATRICES
# ============================================================

print("\nCalculating gene prevalence...")

species_matrices = {}
long_rows = []
summary_rows = []
denominator_rows = []

for species in SPECIES:

    sub = data.loc[
        data["Species"].eq(species)
    ]

    if sub.empty:
        raise ValueError(
            "No genomes for " + species
        )

    occurrences = Counter()

    for genes in sub["Gene_Set"]:
        occurrences.update(genes)

    genes = sorted(
        occurrences,
        key=lambda gene: (
            class_rank[gene_class[gene]],
            -occurrences[gene],
            gene.lower()
        )
    )

    columns = [
        (source, inc)
        for source in SOURCES
        for inc in INCOMES
    ]

    matrix = np.full(
        (len(genes), 9),
        np.nan,
        dtype=np.float32
    )

    for j, (source, inc) in enumerate(columns):

        group = sub.loc[
            sub["Source"].eq(source)
            & sub["Income_Group"].eq(inc)
        ]

        n = len(group)

        counts = Counter()
        for gene_set in group["Gene_Set"]:
            counts.update(gene_set)

        denominator_rows.append({
            "Species": species,
            "Source": source,
            "Income": inc,
            "Genomes": n
        })

        for i, gene in enumerate(genes):

            positives = counts[gene]

            if n > 0:
                prevalence = 100 * positives / n
                matrix[i, j] = prevalence
            else:
                prevalence = np.nan

            long_rows.append({
                "Species": species,
                "AMR_Class": gene_class[gene],
                "Gene": gene,
                "Source": source,
                "Income": inc,
                "Carrier_Genomes": positives,
                "Total_Genomes": n,
                "Prevalence_Percent": prevalence
            })

    species_matrices[species] = (
        genes,
        matrix
    )

    summary_rows.append({
        "Species": species,
        "Genomes": len(sub),
        "Unique_AMR_Genes": len(genes),
        "Figure_Count": 1
    })

    print(
        species,
        "| genomes:", len(sub),
        "| genes:", len(genes)
    )

prevalence_df = pd.DataFrame(long_rows)
summary_df = pd.DataFrame(summary_rows)
denominator_df = pd.DataFrame(denominator_rows)

# ============================================================
# SAVE EXCEL BEFORE RENDERING
# ============================================================

methods = pd.DataFrame([
    (
        "Figure",
        "Figure 4 - AMR gene prevalence"
    ),
    (
        "Figures per species",
        1
    ),
    (
        "Species",
        7
    ),
    (
        "Source groups",
        "Human; Animal; Environment"
    ),
    (
        "Income",
        "H; M; L+LM"
    ),
    (
        "Outcome",
        "Gene-carrier prevalence (%)"
    ),
    (
        "Calculation",
        "100 x positive genomes / all genomes "
        "in species-source-income group"
    ),
    (
        "Class annotation",
        "Verified mapping if supplied; "
        "otherwise provisional gene-name rules"
    ),
    (
        "Missing groups",
        "Grey cells"
    ),
    (
        "Statistical testing",
        "None; descriptive prevalence"
    ),
    (
        "Created",
        datetime.now().isoformat()
    )
], columns=["Parameter", "Value"])

with pd.ExcelWriter(
    EXCEL,
    engine="openpyxl"
) as writer:

    prevalence_df.to_excel(
        writer,
        sheet_name="GENE_PREVALENCE",
        index=False
    )

    mapping.to_excel(
        writer,
        sheet_name="GENE_CLASS_MAPPING",
        index=False
    )

    summary_df.to_excel(
        writer,
        sheet_name="SPECIES_SUMMARY",
        index=False
    )

    denominator_df.to_excel(
        writer,
        sheet_name="GROUP_DENOMINATORS",
        index=False
    )

    methods.to_excel(
        writer,
        sheet_name="METHODS",
        index=False
    )

print("\nExcel saved:", EXCEL)

# ============================================================
# DRAW ONE COMPLETE HEATMAP PER SPECIES
# ============================================================

print("\nGenerating seven continuous heatmaps...")

generated = []

for number, species in enumerate(SPECIES, start=1):

    genes, matrix = species_matrices[species]

    nrows = len(genes)

    # Dynamic height; NEVER paginate.
    height = min(
        MAX_FIGURE_HEIGHT,
        max(10.5, 4.5 + nrows * 0.087)
    )

    width = 14.5

    fig = plt.figure(
        figsize=(width, height),
        facecolor="white"
    )

    # Separate axes prevent label/strip overlap.
    ax_names = fig.add_axes([
        .025, .13, .195, .715
    ])

    ax_class = fig.add_axes([
        .224, .13, .017, .715
    ])

    ax_heat = fig.add_axes([
        .257, .13, .690, .715
    ])

    for axis in [ax_names, ax_class]:

        axis.set_xlim(0, 1)
        axis.set_ylim(
            nrows - .5, -.5
        )

        axis.axis("off")

    image = ax_heat.imshow(
        np.ma.masked_invalid(matrix),
        cmap=CMAP,
        vmin=0,
        vmax=100,
        interpolation="nearest",
        aspect="auto"
    )

    ax_heat.set_xlim(-.5, 8.5)
    ax_heat.set_ylim(
        nrows - .5, -.5
    )

    ax_heat.set_xticks([])
    ax_heat.set_yticks([])

    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    # Adjust gene font with dataset size.
    if nrows <= 100:
        gene_font = 9.0
    elif nrows <= 250:
        gene_font = 8.1
    elif nrows <= 450:
        gene_font = 7.3
    else:
        gene_font = 6.6

    # ========================================================
    # GENE LABELS AND THIN CLASS STRIP
    # ========================================================

    classes = [
        gene_class[gene]
        for gene in genes
    ]

    for i, gene in enumerate(genes):

        ax_names.text(
            .98,
            i,
            gene,
            ha="right",
            va="center",
            fontsize=gene_font,
            color="#183153"
        )

        ax_class.add_patch(
            Rectangle(
                (0, i-.5),
                1,
                1,
                facecolor=CLASS_COLORS[
                    gene_class[gene]
                ],
                edgecolor="none"
            )
        )

        if i > 0 and classes[i] != classes[i-1]:

            ax_heat.axhline(
                i-.5,
                color="#516070",
                linewidth=.8
            )

            ax_class.axhline(
                i-.5,
                color="white",
                linewidth=.8
            )

    # ========================================================
    # PREVALENCE VALUES
    #
    # For small matrices, show exact cell values.
    # For very large matrices, omit printed values
    # to preserve legibility.
    # Full numerical values are retained in Excel.
    # ========================================================

    if nrows <= 85:

        for i in range(nrows):
            for j in range(9):

                value = matrix[i, j]

                if np.isfinite(value):

                    color = (
                        "white"
                        if value >= 60
                        else "#1D344F"
                    )

                    text_value = f"{value:.1f}"

                else:

                    color = "#697586"
                    text_value = "—"

                ax_heat.text(
                    j,
                    i,
                    text_value,
                    ha="center",
                    va="center",
                    fontsize=8.0,
                    color=color
                )

    # ========================================================
    # TOP SOURCE AND INCOME STRIPS
    # ========================================================

    for source_index, source in enumerate(SOURCES):

        start = source_index * 3

        ax_heat.add_patch(
            Rectangle(
                (start-.5, -2.15),
                3,
                .55,
                facecolor=SOURCE_COLORS[source],
                edgecolor="white",
                linewidth=.8,
                clip_on=False
            )
        )

        ax_heat.text(
            start+1,
            -1.88,
            source,
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
            color="white",
            clip_on=False
        )

        for income_index, inc in enumerate(INCOMES):

            x = start + income_index

            ax_heat.add_patch(
                Rectangle(
                    (x-.5, -1.35),
                    1,
                    .53,
                    facecolor=INCOME_COLORS[inc],
                    edgecolor="white",
                    linewidth=.45,
                    clip_on=False
                )
            )

            ax_heat.text(
                x,
                -1.09,
                inc,
                ha="center",
                va="center",
                fontsize=9.5,
                fontweight="bold",
                color="white",
                clip_on=False
            )

    for boundary in [2.5, 5.5]:

        ax_heat.axvline(
            boundary,
            color="white",
            linewidth=2.1,
            zorder=6
        )

    # ========================================================
    # TITLES
    # ========================================================

    fig.text(
        .5,
        .970,
        species,
        ha="center",
        fontsize=19.5,
        fontstyle="italic",
        fontweight="bold",
        color="#17345C"
    )

    fig.text(
        .5,
        .947,
        "AMR gene carriage across One Health "
        "sources and country-income groups",
        ha="center",
        fontsize=11.0,
        color="#516078"
    )

    fig.text(
        .025,
        .865,
        "AMR genes",
        fontsize=11,
        fontweight="bold",
        color="#17345C"
    )

    fig.text(
        .225,
        .865,
        "Class",
        fontsize=10,
        color="#17345C"
    )

    # ========================================================
    # COLORBAR
    # ========================================================

    cax = fig.add_axes([
        .55, .077, .345, .010
    ])

    cb = fig.colorbar(
        image,
        cax=cax,
        orientation="horizontal"
    )

    cb.set_label(
        "AMR gene-carrier prevalence (%)",
        fontsize=10
    )

    cb.set_ticks(
        [0, 20, 40, 60, 80, 100]
    )

    # ========================================================
    # AMR CLASS LEGEND
    # ========================================================

    present = [
        cls for cls in CLASS_ORDER
        if cls in set(classes)
    ]

    handles = [
        Patch(
            facecolor=CLASS_COLORS[cls],
            edgecolor="none",
            label=cls
        )
        for cls in present
    ]

    fig.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(.025, .017),
        ncol=4,
        frameon=False,
        fontsize=8.2
    )

    fig.text(
        .5,
        .006,
        "Each cell indicates the percentage of genomes "
        "carrying the indicated AMR gene. "
        "Grey = no genomes in that group.",
        ha="center",
        fontsize=8.4,
        color="#647083"
    )

    # ========================================================
    # SAVE FIGURE
    # ========================================================

    filename = OUT / (
        f"Fig04_{number:02d}_"
        f"{species_tag(species)}_"
        "ALL_AMR_GENES.jpg"
    )

    save_image(fig, filename)

    generated.append(filename)

    print(
        f"[SAVED {number}/7] "
        f"{species} | {nrows} genes",
        flush=True
    )

# ============================================================
# FINAL REPORT
# ============================================================

if len(generated) != 7:
    raise RuntimeError(
        "Expected exactly seven figures"
    )

print("\n" + "=" * 100)
print("FIGURE 4 COMPLETE")
print("=" * 100)

print("Species:", len(SPECIES))
print("AMR genes across the dataset:", len(all_genes))
print("JPEG figures generated:", len(generated))
print("Resolution:", DPI, "DPI")

print("\nFIGURES:")

for filename in generated:
    print(filename.name)

print("\nEXCEL:", EXCEL)
print("\nNo gene pagination")
print("No separate heatmaps by AMR class")
print("No statistical recalculation")
print("=" * 100)
