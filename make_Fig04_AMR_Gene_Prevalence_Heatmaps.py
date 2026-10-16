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
# FIGURE 4
# AMR GENE PREVALENCE HEATMAPS
# 7 SPECIES-SPECIFIC FIGURES
#
# USER REQUIREMENTS
# - prevalence only
# - no p-values
# - no pairwise comparison panel
# - one heatmap per species
# - 600 DPI JPEG
# - one Excel workbook
# - different colors for source strips and class strips
#   so that no similar strip colors appear
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

OUTDIR = ROOT / "Fig_04_AMR_Gene_Prevalence_Heatmaps"
OUTDIR.mkdir(parents=True, exist_ok=True)

EXCEL_OUT = OUTDIR / "Fig04_AMR_Gene_Prevalence_Heatmaps_Data.xlsx"

# ------------------------------------------------------------
# SPECIES ORDER
# ------------------------------------------------------------
SPECIES_ORDER = [
    "Enterobacter hormaechei",
    "Enterobacter cloacae",
    "Enterobacter roggenkampii",
    "Enterobacter asburiae",
    "Enterobacter kobei",
    "Enterobacter ludwigii",
    "Enterobacter bugandensis",
]

SOURCES = ["Human", "Animal", "Environment"]
INCOMES = ["H", "M", "L+LM"]

# ------------------------------------------------------------
# DISTINCT STRIP COLORS
# IMPORTANT:
# Source-strip colors are intentionally different from
# AMR class-strip colors.
# ------------------------------------------------------------

# Source strips (top source bars)
SOURCE_STRIP_COLORS = {
    "Human": "#C65D00",         # deep burnt orange
    "Animal": "#0077B6",        # strong blue
    "Environment": "#007F5F",   # deep green
}

# Income strips (top sub-bars)
INCOME_STRIP_COLORS = {
    "H": "#D4A017",      # mustard-gold
    "M": "#7B2CBF",      # purple
    "L+LM": "#D1495B",   # crimson-rose
}

# AMR class strip colors
# Deliberately chosen so they do not visually duplicate
# Human/Animal/Environment strip colors.
CLASS_COLORS = {
    "Beta-lactam": "#F4A261",         # sand orange
    "Aminoglycoside": "#8ECAE6",      # pale cyan
    "Quinolone": "#52B788",           # mint green
    "Tetracycline": "#F6D55C",        # yellow
    "Sulfonamide": "#E78AC3",         # pink
    "Trimethoprim": "#B497E7",        # lavender
    "Phenicols": "#DDB892",           # beige brown
    "Macrolide": "#7FD3E6",           # aqua
    "Fosfomycin": "#F6BD60",          # amber
    "Colistin": "#E76F51",            # coral red
    "Rifamycin": "#84A9D8",           # steel blue
    "Bleomycin": "#BDBDBD",           # grey
    "Lincosamide": "#C6A0D6",         # mauve
    "Streptothricin": "#B7C36B",      # olive
    "Other verified class": "#8C8C8C",
    "Unassigned - verify": "#D9D9D9",
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
    "Lincosamide",
    "Streptothricin",
    "Other verified class",
    "Unassigned - verify",
]

# ------------------------------------------------------------
# HEATMAP SETTINGS
# ------------------------------------------------------------
DPI = 600
MAX_GENES_PER_PAGE = 60  # to keep 600-DPI figures readable and safe
NOT_AVAILABLE_COLOR = "#CFCFCF"

PREV_CMAP = LinearSegmentedColormap.from_list(
    "prev_blue",
    [
        "#F7FBFF",
        "#DCEAF7",
        "#AED4F3",
        "#6FB3E9",
        "#2E86DE",
        "#0B4F9C",
    ]
)
PREV_CMAP.set_bad(NOT_AVAILABLE_COLOR)

# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------
def parse_genes(value):
    if pd.isna(value):
        return set()
    return {
        x.strip()
        for x in str(value).split(";")
        if x.strip()
    }

def classify_gene(gene):
    g = gene.lower().strip()

    if g.startswith("bla"):
        return "Beta-lactam"
    if re.match(r"^(aac|aph|ant)(\(|-|_)", g):
        return "Aminoglycoside"
    if re.match(r"^(qnr|oqx)", g):
        return "Quinolone"
    if re.match(r"^tet(\(|[a-z])", g):
        return "Tetracycline"
    if re.match(r"^sul[123]($|[-_])", g):
        return "Sulfonamide"
    if g.startswith("dfr"):
        return "Trimethoprim"
    if re.match(r"^(cat|cml|flo)", g):
        return "Phenicols"
    if re.match(r"^(erm|mph|mef|msr)", g):
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

def short_tag(name):
    return name.replace("Enterobacter ", "E_").replace(" ", "_")

def save_and_verify(fig, outpath):
    fig.savefig(
        outpath,
        dpi=DPI,
        format="jpeg",
        facecolor="white",
        pil_kwargs={
            "quality": 94,
            "subsampling": 0,
            "optimize": True
        }
    )
    plt.close(fig)
    gc.collect()

    with Image.open(outpath) as im:
        im.verify()

    if not outpath.is_file() or outpath.stat().st_size == 0:
        raise RuntimeError(f"Failed to create output image: {outpath}")

def find_sheet(xls, preferred):
    if preferred in xls.sheet_names:
        return preferred
    return xls.sheet_names[0]

# ------------------------------------------------------------
# LOAD MASTER
# ------------------------------------------------------------
print("=" * 110)
print("FIGURE 4 - AMR GENE PREVALENCE HEATMAPS")
print("=" * 110)

if not MASTER_FILE.is_file():
    raise FileNotFoundError(f"Master file not found:\n{MASTER_FILE}")

if not INCOME_FILE.is_file():
    raise FileNotFoundError(f"Income file not found:\n{INCOME_FILE}")

master_xls = pd.ExcelFile(MASTER_FILE, engine="openpyxl")
master_sheet = find_sheet(master_xls, "FINAL_MASTER_2000_2026")
master = pd.read_excel(master_xls, sheet_name=master_sheet, dtype={"Assembly_Accession": str})

required_master_columns = [
    "Assembly_Accession",
    "Normalized_Species",
    "Final_Source_Category",
    "Country_Harmonized",
    "AMR_Genes",
]

missing_master = [c for c in required_master_columns if c not in master.columns]
if missing_master:
    raise ValueError(f"Missing columns in master file: {missing_master}")

income_xls = pd.ExcelFile(INCOME_FILE, engine="openpyxl")

if "COUNTRY_INCOME_COUNTS" in income_xls.sheet_names:
    income = pd.read_excel(income_xls, sheet_name="COUNTRY_INCOME_COUNTS")
else:
    income = pd.read_excel(income_xls, sheet_name=income_xls.sheet_names[0])

if "Country_Harmonized" not in income.columns or "Income_Group" not in income.columns:
    raise ValueError("Income workbook must contain Country_Harmonized and Income_Group columns")

income = income.copy()
income["Country_Harmonized"] = income["Country_Harmonized"].astype(str).str.strip()
income["Income_Group"] = income["Income_Group"].astype(str).str.strip()

if income["Country_Harmonized"].duplicated().any():
    raise ValueError("Duplicate country in income classification file")

master = master.copy()
master["Species"] = master["Normalized_Species"].astype(str).str.strip()
master["Source"] = master["Final_Source_Category"].astype(str).str.strip()
master["Country"] = master["Country_Harmonized"].astype(str).str.strip()

master = master.merge(
    income[["Country_Harmonized", "Income_Group"]],
    left_on="Country",
    right_on="Country_Harmonized",
    how="left",
    validate="many_to_one"
)

data = master[
    master["Species"].isin(SPECIES_ORDER)
    & master["Source"].isin(SOURCES)
    & master["Income_Group"].isin(INCOMES)
].copy()

if data.empty:
    raise ValueError("No usable records after filtering")

data["Gene_Set"] = data["AMR_Genes"].apply(parse_genes)

print("Genomes retained for Figure 4:", len(data))
print("Species represented:", data["Species"].nunique())

# ------------------------------------------------------------
# BUILD GENE INVENTORY AND PREVALENCE TABLE
# ------------------------------------------------------------
all_gene_rows = []
gene_map_rows = []
figure_summary_rows = []

for species in SPECIES_ORDER:

    species_df = data[data["Species"] == species].copy()

    if species_df.empty:
        print(f"[WARNING] No data for {species}; skipping")
        continue

    occurrences = Counter()
    for genes in species_df["Gene_Set"]:
        occurrences.update(genes)

    gene_list = sorted(
        occurrences.keys(),
        key=lambda g: (
            CLASS_ORDER.index(classify_gene(g)) if classify_gene(g) in CLASS_ORDER else 999,
            -occurrences[g],
            g.lower(),
        )
    )

    for gene in gene_list:
        gene_map_rows.append({
            "Species": species,
            "Gene": gene,
            "AMR_Class": classify_gene(gene),
            "Total_Carriers_In_Species": occurrences[gene],
        })

        row = {
            "Species": species,
            "Gene": gene,
            "AMR_Class": classify_gene(gene),
        }

        for source in SOURCES:
            for income_group in INCOMES:
                subset = species_df[
                    (species_df["Source"] == source)
                    & (species_df["Income_Group"] == income_group)
                ]

                n_total = len(subset)
                n_carrier = sum(gene in gset for gset in subset["Gene_Set"])
                prevalence = (100.0 * n_carrier / n_total) if n_total > 0 else np.nan

                prefix = f"{source}_{income_group}"
                row[f"{prefix}_N"] = n_total
                row[f"{prefix}_Carriers"] = n_carrier
                row[f"{prefix}_Prevalence"] = prevalence

        all_gene_rows.append(row)

    figure_summary_rows.append({
        "Species": species,
        "Genomes": len(species_df),
        "Unique_AMR_Genes": len(gene_list),
        "Pages": int(np.ceil(len(gene_list) / MAX_GENES_PER_PAGE)) if len(gene_list) else 0,
    })

prevalence_df = pd.DataFrame(all_gene_rows)
gene_map_df = pd.DataFrame(gene_map_rows)
figure_summary_df = pd.DataFrame(figure_summary_rows)

if prevalence_df.empty:
    raise ValueError("No gene prevalence table produced")

# ------------------------------------------------------------
# SAVE EXCEL FIRST
# ------------------------------------------------------------
methods_df = pd.DataFrame([
    ("Analysis", "Species-specific AMR gene prevalence heatmaps"),
    ("Rows", "Individual AMR genes grouped by provisional AMR class"),
    ("Columns", "Human, Animal, Environment, each split into H / M / L+LM"),
    ("Displayed value", "Percentage of genomes carrying the gene within each group"),
    ("Formula", "100 x gene-carrier genomes / total genomes in that species-source-income group"),
    ("Unavailable cells", "Grey, when no genome is available in the group"),
    ("Statistics", "Not shown in Figure 4; prevalence only"),
    ("Master file", str(MASTER_FILE)),
    ("Income file", str(INCOME_FILE)),
    ("Generated", datetime.now().isoformat()),
], columns=["Parameter", "Value"])

with pd.ExcelWriter(EXCEL_OUT, engine="openpyxl") as writer:
    prevalence_df.to_excel(writer, sheet_name="GENE_PREVALENCE", index=False)
    gene_map_df.to_excel(writer, sheet_name="GENE_CLASS_MAP", index=False)
    figure_summary_df.to_excel(writer, sheet_name="FIGURE_SUMMARY", index=False)
    methods_df.to_excel(writer, sheet_name="METHODS", index=False)

print("[SAVED] Excel source file")
print(EXCEL_OUT)

# ------------------------------------------------------------
# RENDERING FUNCTION
# ------------------------------------------------------------
def draw_species_heatmap(species, species_df, page_no, page_df):
    genes = page_df["Gene"].tolist()
    classes = page_df["AMR_Class"].tolist()
    n_genes = len(genes)

    # 9 prevalence columns in fixed order
    col_pairs = []
    for source in SOURCES:
        for income_group in INCOMES:
            col_pairs.append((source, income_group))

    matrix = np.full((n_genes, 9), np.nan, dtype=float)

    for i, (_, row) in enumerate(page_df.iterrows()):
        for j, (source, income_group) in enumerate(col_pairs):
            matrix[i, j] = row[f"{source}_{income_group}_Prevalence"]

    fig_height = max(10.5, min(21.5, 4.8 + n_genes * 0.28))

    fig = plt.figure(figsize=(15.6, fig_height), facecolor="white")

    # Axes
    ax_class = fig.add_axes([0.03, 0.12, 0.13, 0.76])
    ax_gene  = fig.add_axes([0.16, 0.12, 0.12, 0.76])
    ax_heat  = fig.add_axes([0.29, 0.12, 0.66, 0.76])

    # Gene/class axes setup
    for ax in [ax_class, ax_gene]:
        ax.set_xlim(0, 1)
        ax.set_ylim(n_genes - 0.5, -0.5)
        ax.axis("off")

    # Heatmap
    im = ax_heat.imshow(
        np.ma.masked_invalid(matrix),
        cmap=PREV_CMAP,
        vmin=0,
        vmax=100,
        interpolation="nearest",
        aspect="auto"
    )

    ax_heat.set_xlim(-0.5, 8.5)
    ax_heat.set_ylim(n_genes - 0.5, -0.5)
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])
    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    # Row labels and class strip
    for i, (gene, amr_class) in enumerate(zip(genes, classes)):
        ax_gene.text(
            0.98, i, gene,
            ha="right", va="center",
            fontsize=9.0,
            color="#17324D",
            fontstyle="italic" if gene.startswith("bla") or gene.startswith("qnr") or gene.startswith("tet") or gene.startswith("sul") or gene.startswith("dfr") else "normal"
        )

        ax_class.add_patch(
            Rectangle(
                (0, i - 0.5), 1, 1,
                facecolor=CLASS_COLORS.get(amr_class, "#D9D9D9"),
                edgecolor="white",
                linewidth=0.4
            )
        )

        # class label only at first row of each class block
        if i == 0 or classes[i - 1] != amr_class:
            # determine span
            end = i
            while end + 1 < n_genes and classes[end + 1] == amr_class:
                end += 1
            center = (i + end) / 2
            ax_class.text(
                0.04, center, amr_class,
                ha="left", va="center",
                fontsize=9.5,
                fontweight="bold",
                color="#0F2440"
            )
            ax_heat.axhline(i - 0.5, color="#425466", linewidth=0.9)

    ax_heat.axhline(n_genes - 0.5, color="#425466", linewidth=0.9)

    # Cell values
    for i in range(n_genes):
        for j in range(9):
            val = matrix[i, j]
            if np.isfinite(val):
                txt_color = "white" if val >= 55 else "#17324D"
                ax_heat.text(
                    j, i, f"{val:.1f}",
                    ha="center", va="center",
                    fontsize=8.2,
                    color=txt_color
                )
            else:
                ax_heat.text(
                    j, i, "–",
                    ha="center", va="center",
                    fontsize=9.5,
                    color="#4B5563",
                    fontweight="bold"
                )

    # Top strips and labels
    # Income strips
    for j, (source, income_group) in enumerate(col_pairs):
        ax_heat.add_patch(
            Rectangle(
                (j - 0.5, -1.00), 1, 0.34,
                facecolor=INCOME_STRIP_COLORS[income_group],
                edgecolor="white",
                linewidth=0.6,
                clip_on=False
            )
        )
        ax_heat.text(
            j, -1.20, income_group,
            ha="center", va="center",
            fontsize=9.5,
            color="white",
            fontweight="bold",
            clip_on=False
        )

    # Source strips
    for block_idx, source in enumerate(SOURCES):
        start = block_idx * 3
        ax_heat.add_patch(
            Rectangle(
                (start - 0.5, -1.90), 3, 0.40,
                facecolor=SOURCE_STRIP_COLORS[source],
                edgecolor="white",
                linewidth=0.8,
                clip_on=False
            )
        )
        ax_heat.text(
            start + 1, -2.13, source,
            ha="center", va="center",
            fontsize=11.5,
            color="white",
            fontweight="bold",
            clip_on=False
        )

    # Left headers
    ax_class.add_patch(
        Rectangle(
            (0, -2.35), 1, 0.85,
            facecolor="#0B2E59",
            edgecolor="white",
            linewidth=0.8,
            clip_on=False
        )
    )
    ax_class.text(
        0.5, -1.92, "AMR class",
        ha="center", va="center",
        fontsize=11.5,
        color="white",
        fontweight="bold",
        clip_on=False
    )

    ax_gene.add_patch(
        Rectangle(
            (0, -2.35), 1, 0.85,
            facecolor="#0B2E59",
            edgecolor="white",
            linewidth=0.8,
            clip_on=False
        )
    )
    ax_gene.text(
        0.5, -1.92, "AMR gene",
        ha="center", va="center",
        fontsize=11.5,
        color="white",
        fontweight="bold",
        clip_on=False
    )

    # Titles
    fig.text(
        0.5, 0.966,
        species,
        ha="center",
        fontsize=22,
        fontweight="bold",
        fontstyle="italic",
        color="#102A43"
    )

    fig.text(
        0.5, 0.935,
        "AMR gene carriage prevalence across One Health sources and country-income groups",
        ha="center",
        fontsize=12.5,
        color="#334E68"
    )

    total_pages = int(np.ceil(len(species_df) / MAX_GENES_PER_PAGE))
    fig.text(
        0.5, 0.912,
        f"Figure 4 | {species} | Page {page_no + 1} of {total_pages}",
        ha="center",
        fontsize=10,
        color="#627D98"
    )

    # Colorbar
    cax = fig.add_axes([0.62, 0.050, 0.22, 0.020])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("Gene-carrier prevalence (%)", fontsize=10.5, color="#17324D")
    cb.ax.tick_params(labelsize=9, colors="#17324D")

    # Not available legend
    fig.patches.append(
        Rectangle(
            (0.86, 0.048), 0.03, 0.024,
            transform=fig.transFigure,
            facecolor=NOT_AVAILABLE_COLOR,
            edgecolor="none"
        )
    )
    fig.text(
        0.895, 0.059,
        "Not available",
        ha="left", va="center",
        fontsize=10,
        color="#17324D"
    )

    # Class legend
    present_classes = [c for c in CLASS_ORDER if c in set(classes)]
    class_handles = [
        Patch(facecolor=CLASS_COLORS[c], edgecolor="none", label=c)
        for c in present_classes
    ]

    fig.legend(
        handles=class_handles,
        loc="lower left",
        bbox_to_anchor=(0.03, 0.01),
        ncol=min(5, len(class_handles)),
        frameon=False,
        fontsize=9.0,
        title="AMR class",
        title_fontsize=10.5
    )

    outname = OUTDIR / f"Fig04_{short_tag(species)}_AMR_Gene_Prevalence_Page_{page_no + 1:02d}.jpg"
    save_and_verify(fig, outname)
    return outname

# ------------------------------------------------------------
# MAKE FIGURES
# ------------------------------------------------------------
print("\nRendering figures...")

created_files = []

for species in SPECIES_ORDER:
    sp_df = prevalence_df[prevalence_df["Species"] == species].copy()

    if sp_df.empty:
        print(f"[SKIPPED] {species} (no data)")
        continue

    # sort rows by class then by max prevalence
    sp_df["Class_Order"] = sp_df["AMR_Class"].apply(
        lambda x: CLASS_ORDER.index(x) if x in CLASS_ORDER else 999
    )

    prev_cols = [c for c in sp_df.columns if c.endswith("_Prevalence")]
    sp_df["Max_Prevalence"] = sp_df[prev_cols].max(axis=1, skipna=True)

    sp_df = sp_df.sort_values(
        by=["Class_Order", "Max_Prevalence", "Gene"],
        ascending=[True, False, True]
    ).reset_index(drop=True)

    pages = [
        sp_df.iloc[i:i + MAX_GENES_PER_PAGE].copy()
        for i in range(0, len(sp_df), MAX_GENES_PER_PAGE)
    ]

    print(f"{species}: {len(sp_df)} genes, {len(pages)} page(s)")

    for page_no, page_df in enumerate(pages):
        img = draw_species_heatmap(species, sp_df, page_no, page_df)
        created_files.append(img)
        print(f"[SAVED] {img.name}")

# ------------------------------------------------------------
# FINAL REPORT
# ------------------------------------------------------------
print("\n" + "=" * 110)
print("FIGURE 4 COMPLETED")
print("=" * 110)
print("Species requested:", len(SPECIES_ORDER))
print("Species with figures:", prevalence_df['Species'].nunique())
print("JPEG files created:", len(created_files))
print("Excel source file:", EXCEL_OUT)
print("Figure folder:", OUTDIR)
print("Resolution:", DPI, "DPI")
print("Note: source-strip colors and class-strip colors were chosen to remain visually distinct.")
print("=" * 110)