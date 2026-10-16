from pathlib import Path
import gc

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle, Patch
from matplotlib import font_manager

# ============================================================
# FIGURE 4
# ONE COMBINED HEATMAP FOR ALL SPECIES
#
# INPUT:
#   Existing source workbook created previously:
#   Fig04_ALL_GENES_PREVALENCE_MASTER.xlsx
#
# OUTPUT:
#   1 JPEG (combined heatmap)
#   1 Excel source workbook
#
# STRUCTURE:
#   Species strip | AMR class strip | gene names | 9 prevalence columns | total carrier bar
#
# NOTE:
#   For this combined figure, cell values are NOT printed,
#   because the figure would become unreadable.
#   Exact values are preserved in the Excel workbook.
# ============================================================

ROOT = Path(r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES\Fig_04_AMR_Gene_Prevalence_Heatmaps")

INPUT = ROOT / "Fig04_ALL_GENES_PREVALENCE_MASTER.xlsx"
OUT_JPG = ROOT / "Fig04_Combined_All_Species_AMR_Gene_Prevalence.jpg"
OUT_XLSX = ROOT / "Fig04_Combined_All_Species_AMR_Gene_Prevalence_SourceData.xlsx"

DPI = 600

SPECIES_ORDER = [
    "Enterobacter hormaechei",
    "Enterobacter cloacae",
    "Enterobacter roggenkampii",
    "Enterobacter asburiae",
    "Enterobacter kobei",
    "Enterobacter ludwigii",
    "Enterobacter bugandensis"
]

SOURCES = ["Human", "Animal", "Environment"]
INCOMES = ["H", "M", "L+LM"]

SOURCE_COLORS = {
    "Human": "#C65D00",
    "Animal": "#0077B6",
    "Environment": "#007F5F"
}

INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82"
}

SPECIES_COLORS = {
    "Enterobacter hormaechei": "#8E6C8A",
    "Enterobacter cloacae": "#6B8E23",
    "Enterobacter roggenkampii": "#6C5B7B",
    "Enterobacter asburiae": "#B56576",
    "Enterobacter kobei": "#5C80BC",
    "Enterobacter ludwigii": "#BC6C25",
    "Enterobacter bugandensis": "#3D5A80"
}

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

# ============================================================
# FONT — Nature Microbiology style: clean sans, compact
# ============================================================

installed = {f.name for f in font_manager.fontManager.ttflist}
FONT = next((x for x in ["Arial", "Liberation Sans", "DejaVu Sans"] if x in installed), "DejaVu Sans")

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 8.5,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42
})

CMAP = LinearSegmentedColormap.from_list(
    "amr_prev",
    [
        "#FFFFFF",
        "#E3F0FA",
        "#B9DAF1",
        "#7CB9E3",
        "#3C91CD",
        "#1768AE",
        "#123C78"
    ]
)
CMAP.set_bad("#D9DDE2")

# ============================================================
# LOAD EXISTING WORKBOOK
# ============================================================

print("=" * 110)
print("FIGURE 4 - ONE COMBINED ALL-SPECIES HEATMAP")
print("=" * 110)

if not INPUT.is_file():
    raise FileNotFoundError(f"Input workbook not found:\n{INPUT}")

df = pd.read_excel(INPUT, sheet_name="GENE_PREVALENCE", engine="openpyxl")

required = [
    "Species", "AMR_Class", "Gene", "Source", "Income",
    "Carrier_Genomes", "Total_Genomes", "Prevalence_Percent"
]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns in input workbook: {missing}")

for col in ["Species", "AMR_Class", "Gene", "Source", "Income"]:
    df[col] = df[col].astype(str).str.strip()

for col in ["Carrier_Genomes", "Total_Genomes", "Prevalence_Percent"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df = df[
    df["Species"].isin(SPECIES_ORDER)
    & df["Source"].isin(SOURCES)
    & df["Income"].isin(INCOMES)
].copy()

if df.empty:
    raise ValueError("No usable rows found")

# One class per Species+Gene only
conflict_check = df.groupby(["Species", "Gene"])["AMR_Class"].nunique()
if (conflict_check > 1).any():
    raise ValueError("Conflicting AMR class assignments detected for some genes")

# ============================================================
# BUILD COMBINED ROW TABLE
# ============================================================

row_records = []

for species in SPECIES_ORDER:
    sub = df[df["Species"] == species].copy()
    if sub.empty:
        continue

    # total absolute carrier value across all 9 groups
    totals = sub.groupby(["Gene", "AMR_Class"], as_index=False).agg(
        Total_Carrier_Genomes=("Carrier_Genomes", "sum"),
        Max_Prevalence=("Prevalence_Percent", "max")
    )

    totals["Species"] = species
    totals["Class_Order"] = totals["AMR_Class"].map(lambda x: CLASS_ORDER.index(x) if x in CLASS_ORDER else 999)

    totals = totals.sort_values(
        by=["Class_Order", "Total_Carrier_Genomes", "Max_Prevalence", "Gene"],
        ascending=[True, False, False, True]
    ).reset_index(drop=True)

    for _, row in totals.iterrows():
        gene = row["Gene"]
        amr_class = row["AMR_Class"]

        record = {
            "Species": species,
            "AMR_Class": amr_class,
            "Gene": gene,
            "Total_Carrier_Genomes": int(row["Total_Carrier_Genomes"]),
        }

        for source in SOURCES:
            for income in INCOMES:
                subset = sub[
                    (sub["Gene"] == gene)
                    & (sub["Source"] == source)
                    & (sub["Income"] == income)
                ]
                if len(subset) != 1:
                    raise ValueError(f"Expected one row for {species} | {gene} | {source} | {income}")
                record[f"{source}_{income}_Prev"] = subset["Prevalence_Percent"].iloc[0]
                record[f"{source}_{income}_Carrier"] = subset["Carrier_Genomes"].iloc[0]
                record[f"{source}_{income}_Total"] = subset["Total_Genomes"].iloc[0]

        row_records.append(record)

rows = pd.DataFrame(row_records)

if rows.empty:
    raise ValueError("Combined row table is empty")

print("Total plotted rows:", len(rows))
print("Species represented:", rows["Species"].nunique())

# ============================================================
# SAVE SOURCE EXCEL
# ============================================================

methods = pd.DataFrame([
    ("Figure", "Figure 4 combined all-species AMR gene prevalence heatmap"),
    ("Rows", "Species + individual AMR genes"),
    ("Heatmap value", "Prevalence (%)"),
    ("Right-side bar", "Absolute total number of carrier genomes across the 9 source-income groups"),
    ("Columns", "Human, Animal, Environment; each subdivided into H, M, L+LM"),
    ("Cell values shown", "No, for readability in the combined figure"),
    ("Reason", "Exact values retained in the Excel source workbook"),
], columns=["Parameter", "Value"])

with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
    rows.to_excel(writer, sheet_name="COMBINED_HEATMAP_DATA", index=False)
    methods.to_excel(writer, sheet_name="METHODS", index=False)

print("[SAVED] Source Excel:")
print(OUT_XLSX)

# ============================================================
# MATRIX
# ============================================================

heat_cols = []
for source in SOURCES:
    for income in INCOMES:
        heat_cols.append(f"{source}_{income}_Prev")

matrix = rows[heat_cols].to_numpy(dtype=float)
bar_values = rows["Total_Carrier_Genomes"].to_numpy(dtype=float)
genes = rows["Gene"].tolist()
classes = rows["AMR_Class"].tolist()
species_vec = rows["Species"].tolist()

n_rows = len(rows)

# ------------------------------------------------------------
# Figure size
# Large digital figure, still one JPEG.
# This is suitable for zooming and supplementary-quality display.
# ------------------------------------------------------------
FIG_WIDTH = 10.0
FIG_HEIGHT = min(74.0, max(28.0, 6.0 + n_rows * 0.035))

print(f"Figure size: {FIG_WIDTH:.2f} x {FIG_HEIGHT:.2f} inches at {DPI} DPI")

fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT), facecolor="white")

bottom = 0.085
height = 0.84

# Axes:
# species strip | class strip | gene names | heatmap | barplot
ax_species = fig.add_axes([0.020, bottom, 0.030, height])
ax_class   = fig.add_axes([0.055, bottom, 0.018, height])
ax_names   = fig.add_axes([0.078, bottom, 0.270, height])
ax_heat    = fig.add_axes([0.355, bottom, 0.340, height])
ax_bar     = fig.add_axes([0.725, bottom, 0.235, height])

for ax in [ax_species, ax_class, ax_names]:
    ax.set_xlim(0, 1)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.axis("off")

# ============================================================
# LEFT SPECIES STRIP
# ============================================================

for i, sp in enumerate(species_vec):
    ax_species.add_patch(
        Rectangle(
            (0, i - 0.5), 1, 1,
            facecolor=SPECIES_COLORS[sp],
            edgecolor="none"
        )
    )

# Species boundaries and labels
species_blocks = []
start = 0
for i in range(1, n_rows + 1):
    if i == n_rows or species_vec[i] != species_vec[start]:
        species_blocks.append((species_vec[start], start, i - 1))
        start = i

for sp, s, e in species_blocks:
    center = (s + e) / 2
    ax_species.text(
        0.5, center,
        sp.replace("Enterobacter ", "E. "),
        ha="center", va="center",
        rotation=90,
        fontsize=8.5,
        color="white",
        fontweight="bold"
    )

# ============================================================
# CLASS STRIP + GENE NAMES
# ============================================================

for i, (cls, gene) in enumerate(zip(classes, genes)):
    ax_class.add_patch(
        Rectangle(
            (0, i - 0.5), 1, 1,
            facecolor=CLASS_COLORS.get(cls, "#D9D9D9"),
            edgecolor="none"
        )
    )
    ax_names.text(
        0.99, i, gene,
        ha="right", va="center",
        fontsize=4.7,
        color="#17324D"
    )

# Species separators
for _, s, e in species_blocks[:-1]:
    y = e + 0.5
    for ax in [ax_species, ax_class, ax_names, ax_heat, ax_bar]:
        ax.axhline(y, color="#2F3E4E", linewidth=0.8, zorder=10)

# Class separators within species
for i in range(1, n_rows):
    if species_vec[i] == species_vec[i - 1] and classes[i] != classes[i - 1]:
        y = i - 0.5
        for ax in [ax_class, ax_heat, ax_bar]:
            ax.axhline(y, color="#A0A9B2", linewidth=0.35, zorder=9)

# ============================================================
# HEATMAP
# ============================================================

im = ax_heat.imshow(
    np.ma.masked_invalid(matrix),
    cmap=CMAP,
    vmin=0,
    vmax=100,
    aspect="auto",
    interpolation="nearest",
    resample=False
)

ax_heat.set_xlim(-0.5, 8.5)
ax_heat.set_ylim(n_rows - 0.5, -0.5)
ax_heat.set_xticks([])
ax_heat.set_yticks([])
for spine in ax_heat.spines.values():
    spine.set_visible(False)

for x in [2.5, 5.5]:
    ax_heat.axvline(x, color="white", linewidth=1.2, zorder=10)

# ============================================================
# RIGHT ABSOLUTE TOTAL BAR
# ============================================================

ypos = np.arange(n_rows)
ax_bar.barh(
    ypos,
    bar_values,
    height=0.72,
    color="#274C77",
    edgecolor="none"
)

ax_bar.set_ylim(n_rows - 0.5, -0.5)
ax_bar.set_yticks([])
ax_bar.set_xlabel("Total carrier genomes (absolute n)", fontsize=8.7, color="#17324D", labelpad=6)
ax_bar.tick_params(axis="x", labelsize=7.8, colors="#17324D")
ax_bar.grid(axis="x", color="#D0D7DE", linewidth=0.5, alpha=0.8)
ax_bar.set_axisbelow(True)
ax_bar.spines["top"].set_visible(False)
ax_bar.spines["right"].set_visible(False)
ax_bar.spines["left"].set_visible(False)
ax_bar.spines["bottom"].set_color("#8091A5")

# ============================================================
# TOP STRIPS
# ============================================================

ax_top = fig.add_axes([0.355, 0.935, 0.340, 0.030])
ax_top.set_xlim(-0.5, 8.5)
ax_top.set_ylim(0, 2)
ax_top.axis("off")

for source_idx, source in enumerate(SOURCES):
    start = source_idx * 3
    ax_top.add_patch(
        Rectangle(
            (start - 0.5, 1.05), 3, 0.90,
            facecolor=SOURCE_COLORS[source],
            edgecolor="white",
            linewidth=0.7
        )
    )
    ax_top.text(
        start + 1, 1.50, source,
        ha="center", va="center",
        fontsize=8.8, color="white", fontweight="bold"
    )

    for income_idx, income in enumerate(INCOMES):
        j = start + income_idx
        ax_top.add_patch(
            Rectangle(
                (j - 0.5, 0.05), 1, 0.85,
                facecolor=INCOME_COLORS[income],
                edgecolor="white",
                linewidth=0.5
            )
        )
        ax_top.text(
            j, 0.47, income,
            ha="center", va="center",
            fontsize=7.8, color="white", fontweight="bold"
        )

# ============================================================
# TITLES
# ============================================================

fig.text(
    0.50, 0.986,
    "Combined AMR gene prevalence landscape across Enterobacter species",
    ha="center", va="top",
    fontsize=12.8, fontweight="bold", color="#102A43"
)

fig.text(
    0.50, 0.972,
    "Rows are species-specific AMR genes; heatmap shows prevalence (%) and the right-side bar shows the absolute total number of carrier genomes",
    ha="center", va="top",
    fontsize=8.6, color="#486581"
)

fig.text(0.035, 0.918, "Species", ha="center", fontsize=8.4, fontweight="bold", color="#102A43")
fig.text(0.064, 0.918, "Class",   ha="center", fontsize=8.4, fontweight="bold", color="#102A43")
fig.text(0.205, 0.918, "AMR gene", ha="center", fontsize=8.4, fontweight="bold", color="#102A43")
fig.text(0.525, 0.918, "Prevalence (%)", ha="center", fontsize=8.8, fontweight="bold", color="#102A43")
fig.text(0.840, 0.918, "Absolute abundance", ha="center", fontsize=8.8, fontweight="bold", color="#102A43")

# ============================================================
# COLORBAR
# ============================================================

cax = fig.add_axes([0.355, 0.045, 0.200, 0.010])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_ticks([0, 20, 40, 60, 80, 100])
cb.set_label("AMR gene prevalence (%)", fontsize=8.3, color="#17324D", labelpad=3)
cb.ax.tick_params(labelsize=7.4, colors="#17324D")

# ============================================================
# LEGENDS
# ============================================================

# Class legend
present_classes = [c for c in CLASS_ORDER if c in set(classes)]
handles = [Patch(facecolor=CLASS_COLORS[c], edgecolor="none", label=c) for c in present_classes]

fig.legend(
    handles=handles,
    loc="lower left",
    bbox_to_anchor=(0.015, 0.002),
    ncol=4,
    frameon=False,
    fontsize=7.0,
    title="AMR class",
    title_fontsize=7.8,
    columnspacing=0.8,
    handlelength=1.1
)

fig.text(
    0.73, 0.012,
    "Grey heatmap cells indicate that no genome was available\nfor the corresponding species-source-income group.",
    ha="left", va="bottom",
    fontsize=7.2, color="#52606D"
)

# ============================================================
# SAVE
# ============================================================

fig.savefig(
    OUT_JPG,
    dpi=DPI,
    format="jpeg",
    facecolor="white",
    pil_kwargs={
        "quality": 90,
        "subsampling": 2,
        "optimize": False
    }
)

plt.close(fig)
gc.collect()

if not OUT_JPG.is_file():
    raise RuntimeError("Combined figure was not created")

if OUT_JPG.stat().st_size < 10000:
    raise RuntimeError("Combined figure file is too small and likely invalid")

print("\n[SAVED] Combined heatmap:")
print(OUT_JPG)

print("\n" + "=" * 110)
print("FIGURE 4 COMBINED HEATMAP COMPLETED")
print("=" * 110)
print("Rows plotted:", n_rows)
print("Species plotted:", rows['Species'].nunique())
print("Figure:", OUT_JPG)
print("Source data:", OUT_XLSX)
print("Right bar = absolute total carrier genomes (n)")
print("Heatmap = prevalence (%)")
print("Note: exact prevalence values are retained in Excel for readability.")
print("=" * 110)