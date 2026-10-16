
#!/usr/bin/env python3

# ============================================================
# ECC: SEVEN LARGE SPHERICAL NETWORKS IN ONE FRAME
# A4 PORTRAIT PROPORTIONS - ENLARGED MANUSCRIPT FIGURE
#
# Outputs:
#   1 curved-edge JPEG (600 DPI)
#   1 straight-edge JPEG (600 DPI)
#   1 complete Excel workbook
#
# ALL replicons, AMR genes and observed associations
# ALL full replicon and AMR gene names at RIGHT
# TIGHTER right-side legend spacing
# MUCH LARGER spheres, arranged 3 + 2 + 2
# NO top figure title or top panel headings
# ============================================================

from pathlib import Path
from collections import Counter, defaultdict
import math
import colorsys
import gc

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle, Wedge
from matplotlib.colors import to_hex
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ============================================================
# PATHS
# ============================================================

ROOT = Path(
    r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES"
    r"\ECC_SUPPLEMENTARY_AND_SOURCE_DATA"
)

INPUT = ROOT / "ECC_FINAL_2000_2026_MASTER.xlsx"

OUT = ROOT / "FIG_SEVEN_SPHERES_MAXIMUM_SIZE"
OUT.mkdir(parents=True, exist_ok=True)

CURVED_JPG = OUT / "ECC_7_SPHERES_MAXIMUM_CURVED_600DPI.jpg"
STRAIGHT_JPG = OUT / "ECC_7_SPHERES_MAXIMUM_STRAIGHT_600DPI.jpg"
EXCEL = OUT / "ECC_7_SPHERES_MAXIMUM_SOURCE_DATA.xlsx"

# ============================================================
# FIGURE SETTINGS
# ============================================================

# A4 portrait proportions: width / height ~ 0.707
# Enlarged print canvas, not literal A4 paper size.
FIG_W = 30.0
FIG_H = 42.4

DPI = 600

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": [
        "Arial", "Helvetica", "Liberation Sans",
        "DejaVu Sans"
    ],
    "font.size": 10,
    "text.color": "#252525",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.facecolor": "white",
    "pdf.fonttype": 42
})

SPECIES = [
    "Enterobacter hormaechei",
    "Enterobacter cloacae",
    "Enterobacter roggenkampii",
    "Enterobacter asburiae",
    "Enterobacter kobei",
    "Enterobacter ludwigii",
    "Enterobacter bugandensis"
]

SPECIES_COLORS = {
    "Enterobacter hormaechei": "#0072B2",
    "Enterobacter cloacae": "#D55E00",
    "Enterobacter roggenkampii": "#009E73",
    "Enterobacter asburiae": "#CC79A7",
    "Enterobacter kobei": "#7C3AED",
    "Enterobacter ludwigii": "#A16207",
    "Enterobacter bugandensis": "#5B6770"
}

SOURCE_COLORS = {
    "Human": "#E69F00",
    "Animal": "#56B4E9",
    "Environment": "#009E73"
}

AMR_COLORS = {
    "Beta-lactam": "#1D4ED8",
    "Aminoglycoside": "#EA580C",
    "Sulfonamide": "#059669",
    "Trimethoprim": "#9333EA",
    "Tetracycline": "#DB2777",
    "Quinolone": "#0891B2",
    "Phenicol": "#92400E",
    "Macrolide": "#BE185D",
    "Fosfomycin": "#0D9488",
    "Rifamycin": "#A16207",
    "Colistin": "#DC2626",
    "Bleomycin": "#64748B",
    "Other": "#94A3B8"
}

EXCLUDED_COUNTRIES = {
    "not applicable", "puerto rico", "svalbard"
}

# ============================================================
# NETWORK GEOMETRY
# ============================================================

CENTRE_RADIUS = 0.24
REPLICON_RING_RADIUS = 0.72
AMR_INNER_RADIUS = 1.18
AMR_OUTER_RADIUS = 2.07

CURVED_ALPHA = 0.16
STRAIGHT_ALPHA = 0.14

# In-sphere labels; all full labels are at right.
SPHERE_REPLICON_LABELS = 10
SPHERE_AMR_LABELS = 18

# Full right-side legends
REPLICON_COLUMNS = 2
AMR_COLUMNS = 4

REPLICON_FONT = 8.3
AMR_FONT = 8.0

# Layout sizes
LEFT_START = 0.018
LEFT_END = 0.737
RIGHT_START = 0.755
RIGHT_END = 0.988

# ============================================================
# HELPERS
# ============================================================

def parse_set(value):
    if pd.isna(value):
        return set()

    return {
        x.strip()
        for x in str(value).split(";")
        if x.strip()
        and x.strip().lower()
        not in {"nan", "none", "null", ""}
    }


def short_species(value):
    return value.replace("Enterobacter ", "E. ")


def classify_amr(gene):
    g = str(gene).lower()

    if g.startswith("bla"):
        return "Beta-lactam"
    if g.startswith(("aac", "aad", "aph", "ant")):
        return "Aminoglycoside"
    if g.startswith("sul"):
        return "Sulfonamide"
    if g.startswith("dfr"):
        return "Trimethoprim"
    if g.startswith("tet"):
        return "Tetracycline"
    if g.startswith(("qnr", "qep")):
        return "Quinolone"
    if g.startswith(("cat", "flo", "cml")):
        return "Phenicol"
    if g.startswith(("mph", "erm", "ere", "mef")):
        return "Macrolide"
    if g.startswith("fos"):
        return "Fosfomycin"
    if g.startswith(("arr", "rph")):
        return "Rifamycin"
    if g.startswith("mcr"):
        return "Colistin"
    if g.startswith("ble"):
        return "Bleomycin"

    return "Other"


def polar(radius, angle):
    return np.array([
        radius * math.cos(angle),
        radius * math.sin(angle)
    ])


def unit(vector):
    length = np.linalg.norm(vector)
    if length == 0:
        return np.array([1.0, 0.0])
    return vector / length


def make_colors(names):
    result = {}
    golden = 0.618033988749895

    for i, name in enumerate(names):
        h = (0.07 + golden * i) % 1.0
        s = [0.78, 0.62, 0.88][i % 3]
        v = [0.85, 0.72, 0.96][i % 3]

        result[name] = to_hex(
            colorsys.hsv_to_rgb(h, s, v)
        )

    return result


def bezier_points(p0, p1, curvature=0.18, count=14):
    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)

    delta = p1 - p0
    distance = np.linalg.norm(delta)

    if distance < 1e-10:
        return np.vstack([p0, p1])

    normal = np.array([
        -delta[1], delta[0]
    ]) / distance

    sign = (
        1 if p0[0] + p1[1] - p0[1] >= 0
        else -1
    )

    control = (
        (p0 + p1) / 2
        + sign * curvature * distance * normal
    )

    t = np.linspace(0, 1, count)

    return (
        ((1 - t) ** 2)[:, None] * p0
        + (2 * (1 - t) * t)[:, None] * control
        + (t ** 2)[:, None] * p1
    )


def source_ring(ax, x, y, radius, counts):
    total = sum(counts.values())

    if total <= 0:
        return

    start = 90

    for source, colour in SOURCE_COLORS.items():
        n = counts.get(source, 0)

        if n <= 0:
            continue

        angle = 360 * n / total

        ax.add_patch(
            Wedge(
                (x, y), radius,
                start, start + angle,
                width=radius * 0.24,
                facecolor=colour,
                edgecolor="white",
                linewidth=0.25,
                zorder=8
            )
        )

        start += angle


# ============================================================
# LOAD MASTER
# ============================================================

print("=" * 78)
print("ECC LARGE A4-PROPORTION SEVEN-SPHERE NETWORK")
print("=" * 78)

if not INPUT.exists():
    raise FileNotFoundError(INPUT)

required = [
    "Assembly_Accession",
    "Normalized_Species",
    "Country_Harmonized",
    "Final_Source_Category",
    "Plasmid_Replicons",
    "AMR_Genes",
    "BioProject",
    "Verified_Collection_Year"
]

master = pd.read_excel(
    INPUT,
    sheet_name="FINAL_MASTER_2000_2026",
    usecols=required
)

country = (
    master["Country_Harmonized"]
    .fillna("")
    .astype(str)
    .str.strip()
    .str.lower()
)

master = master.loc[
    master["Normalized_Species"].isin(SPECIES)
    & master["Final_Source_Category"].isin(SOURCE_COLORS)
    & ~country.isin(EXCLUDED_COUNTRIES)
].copy()

master["REP_SET"] = (
    master["Plasmid_Replicons"].apply(parse_set)
)

master["GENE_SET"] = (
    master["AMR_Genes"].apply(parse_set)
)

if master.empty:
    raise RuntimeError("No eligible genomes.")

print("Eligible genomes:", len(master))

# ============================================================
# GLOBAL COLOUR KEYS
# ============================================================

global_rep = Counter()
global_gene = Counter()

for values in master["REP_SET"]:
    global_rep.update(values)

for values in master["GENE_SET"]:
    global_gene.update(values)

REP_ORDER = sorted(
    global_rep,
    key=lambda x: (-global_rep[x], x)
)

GENE_ORDER = sorted(
    global_gene,
    key=lambda x: (-global_gene[x], x)
)

REP_COLOR = make_colors(REP_ORDER)

REP_ID = {
    name: f"R{i+1:03d}"
    for i, name in enumerate(REP_ORDER)
}

GENE_ID = {
    name: f"G{i+1:03d}"
    for i, name in enumerate(GENE_ORDER)
}

print("All replicons:", len(REP_ORDER))
print("All AMR genes:", len(GENE_ORDER))

# ============================================================
# SPECIES ANALYSIS
# ============================================================

species_data = {}
summary_rows = []
rep_rows = []
gene_rows = []
association_rows = []

for species in SPECIES:

    df = master.loc[
        master["Normalized_Species"] == species
    ].copy()

    N = len(df)

    if N == 0:
        continue

    reps = Counter()
    genes = Counter()
    pairs = Counter()

    rep_sources = defaultdict(Counter)
    gene_sources = defaultdict(Counter)
    pair_sources = defaultdict(Counter)

    for row in df.itertuples(index=False):
        rr = row.REP_SET
        gg = row.GENE_SET
        src = row.Final_Source_Category

        reps.update(rr)
        genes.update(gg)

        for gene in gg:
            gene_sources[gene][src] += 1

        for rep in rr:
            rep_sources[rep][src] += 1

            for gene in gg:
                pairs[(rep, gene)] += 1
                pair_sources[(rep, gene)][src] += 1

    ordered_reps = sorted(
        reps,
        key=lambda x: (-reps[x], x)
    )

    # Inner ring
    rp = {}

    if ordered_reps:
        angles = np.linspace(
            np.pi / 2,
            5 * np.pi / 2,
            len(ordered_reps),
            endpoint=False
        )

        for i, (rep, angle) in enumerate(
            zip(ordered_reps, angles)
        ):
            r = (
                REPLICON_RING_RADIUS
                + 0.10 * (i % 3) / 2
            )

            rp[rep] = polar(r, angle)

    # Outer AMR class crown
    groups = defaultdict(list)

    for gene in genes:
        groups[classify_amr(gene)].append(gene)

    classes = [
        cls for cls in AMR_COLORS
        if groups[cls]
    ]

    gp = {}

    if classes:
        weights = np.array([
            math.sqrt(len(groups[cls]))
            for cls in classes
        ], dtype=float)

        weights /= weights.sum()

        current = np.pi / 2
        sectors = {}

        for cls, frac in zip(classes, weights):
            width = 2 * np.pi * frac
            sectors[cls] = (current, current + width)
            current += width

        for cls in classes:
            names = sorted(
                groups[cls],
                key=lambda g: (-genes[g], g)
            )

            n = len(names)
            a0, a1 = sectors[cls]
            pad = min(0.05, (a1 - a0) * 0.08)

            if n == 1:
                angles = [(a0 + a1) / 2]
            else:
                angles = np.linspace(
                    a0 + pad,
                    a1 - pad,
                    n
                )

            for i, (gene, angle) in enumerate(
                zip(names, angles)
            ):
                frac = (i + 0.5) / n

                r = (
                    AMR_INNER_RADIUS
                    + (AMR_OUTER_RADIUS-AMR_INNER_RADIUS)
                    * math.sqrt(frac)
                )

                r += 0.09 * math.sin(i * 2.4)
                r += 0.045 * (i % 4)

                gp[gene] = polar(r, angle)

    # Statistics
    local_stats = []

    for (rep, gene), a in pairs.items():
        b = reps[rep] - a
        c = genes[gene] - a
        d = N - a - b - c

        odds, p = fisher_exact(
            [[a, b], [c, d]],
            alternative="two-sided"
        )

        expected = reps[rep] * genes[gene] / N

        item = {
            "Species": species,
            "Replicon_ID": REP_ID[rep],
            "Replicon": rep,
            "Replicon_Color": REP_COLOR[rep],
            "Gene_ID": GENE_ID[gene],
            "AMR_Gene": gene,
            "AMR_Class": classify_amr(gene),
            "Cooccurrence_Genomes": a,
            "Species_Genomes": N,
            "Replicon_Genomes": reps[rep],
            "AMR_Gene_Genomes": genes[gene],
            "Expected_Cooccurrence": expected,
            "Lift": a / expected if expected else np.nan,
            "Odds_Ratio": odds,
            "P_Value": p
        }

        for src in SOURCE_COLORS:
            item[src + "_Cooccurrence"] = (
                pair_sources[(rep, gene)][src]
            )

        local_stats.append(item)

    if local_stats:
        q = multipletests(
            [r["P_Value"] for r in local_stats],
            method="fdr_bh"
        )[1]

        for record, value in zip(local_stats, q):
            record["FDR_q"] = value
            record["FDR_Significant"] = value < 0.05

    association_rows.extend(local_stats)

    for rep in ordered_reps:
        item = {
            "Species": species,
            "Replicon_ID": REP_ID[rep],
            "Replicon": rep,
            "Color_HEX": REP_COLOR[rep],
            "Positive_Genomes": reps[rep],
            "Species_Genomes": N,
            "Prevalence_Percent": 100 * reps[rep] / N,
            "X": rp[rep][0],
            "Y": rp[rep][1]
        }

        for src in SOURCE_COLORS:
            item[src + "_Genomes"] = (
                rep_sources[rep][src]
            )

        rep_rows.append(item)

    for gene in genes:
        item = {
            "Species": species,
            "Gene_ID": GENE_ID[gene],
            "AMR_Gene": gene,
            "AMR_Class": classify_amr(gene),
            "Color_HEX": AMR_COLORS[classify_amr(gene)],
            "Positive_Genomes": genes[gene],
            "Species_Genomes": N,
            "Prevalence_Percent": 100 * genes[gene] / N,
            "X": gp[gene][0],
            "Y": gp[gene][1]
        }

        for src in SOURCE_COLORS:
            item[src + "_Genomes"] = (
                gene_sources[gene][src]
            )

        gene_rows.append(item)

    species_data[species] = {
        "N": N,
        "reps": reps,
        "genes": genes,
        "pairs": pairs,
        "rp": rp,
        "gp": gp,
        "ordered_reps": ordered_reps,
        "rep_sources": rep_sources
    }

    summary_rows.append({
        "Species": species,
        "Genomes": N,
        "Replicons": len(reps),
        "AMR_Genes": len(genes),
        "Observed_Associations": len(pairs)
    })

    print(
        short_species(species),
        "| genomes:", N,
        "| replicons:", len(reps),
        "| AMR genes:", len(genes),
        "| edges:", len(pairs)
    )

# ============================================================
# DRAW ONE SPHERE
# ============================================================

def draw_sphere(ax, species, curved):

    d = species_data[species]

    reps = d["reps"]
    genes = d["genes"]
    pairs = d["pairs"]

    rp = d["rp"]
    gp = d["gp"]

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-2.62, 2.62)
    ax.set_ylim(-2.62, 2.62)

    paths = []
    edge_colors = []
    edge_widths = []

    maximum = max(pairs.values(), default=1)

    for (rep, gene), count in pairs.items():
        p0 = rp[rep]
        p1 = gp[gene]

        if curved:
            points = bezier_points(p0, p1)
        else:
            points = np.vstack([p0, p1])

        paths.append(points)
        edge_colors.append(REP_COLOR[rep])

        edge_widths.append(
            0.13 + 1.25 * math.sqrt(
                count / maximum
            )
        )

    if paths:
        ax.add_collection(
            LineCollection(
                paths,
                colors=edge_colors,
                linewidths=edge_widths,
                alpha=(
                    CURVED_ALPHA
                    if curved else STRAIGHT_ALPHA
                ),
                rasterized=True,
                zorder=1
            )
        )

    # Inner spokes
    for rep in d["ordered_reps"]:
        x, y = rp[rep]

        ax.plot(
            [0, x], [0, y],
            color=REP_COLOR[rep],
            linewidth=0.85,
            alpha=0.58,
            zorder=2
        )

    # Species hub
    ax.add_patch(
        Circle(
            (0, 0),
            CENTRE_RADIUS,
            facecolor=SPECIES_COLORS[species],
            edgecolor="white",
            linewidth=1.4,
            zorder=10
        )
    )

    ax.text(
        0, 0.04,
        short_species(species),
        fontsize=11,
        fontweight="bold",
        color="white",
        ha="center",
        va="center",
        zorder=11
    )

    ax.text(
        0, -0.06,
        f"n={d['N']:,}",
        fontsize=8.5,
        color="white",
        ha="center",
        va="center",
        zorder=11
    )

    # Replicons
    maximum_rep = max(reps.values(), default=1)

    for rep in d["ordered_reps"]:
        x, y = rp[rep]

        r = (
            0.022
            + 0.035 * math.sqrt(
                reps[rep] / maximum_rep
            )
        )

        ax.add_patch(
            Circle(
                (x, y),
                r,
                facecolor=REP_COLOR[rep],
                edgecolor="white",
                linewidth=0.45,
                zorder=5
            )
        )

        source_ring(
            ax, x, y,
            r * 1.25,
            d["rep_sources"][rep]
        )

    # More visible R-IDs
    for rep in d["ordered_reps"][:SPHERE_REPLICON_LABELS]:
        p = rp[rep]
        label = p + 0.075 * unit(p)

        ax.text(
            label[0], label[1],
            REP_ID[rep],
            fontsize=8.0,
            fontweight="bold",
            color=REP_COLOR[rep],
            ha="left" if label[0] >= 0 else "right",
            va="center",
            zorder=12
        )

    # AMR gene nodes
    max_gene = max(genes.values(), default=1)

    for cls, colour in AMR_COLORS.items():
        matching = [
            gene for gene in genes
            if classify_amr(gene) == cls
        ]

        if not matching:
            continue

        coords = np.array([
            gp[gene] for gene in matching
        ])

        sizes = [
            12 + 160 * math.sqrt(
                genes[gene] / max_gene
            )
            for gene in matching
        ]

        ax.scatter(
            coords[:, 0],
            coords[:, 1],
            s=sizes,
            c=colour,
            marker="s",
            edgecolors="white",
            linewidths=0.35,
            alpha=0.96,
            rasterized=True,
            zorder=6
        )

    # G-IDs
    for gene, count in genes.most_common(
        SPHERE_AMR_LABELS
    ):
        p = gp[gene]
        label = p + 0.08 * unit(p)

        ax.text(
            label[0], label[1],
            GENE_ID[gene],
            fontsize=7.4,
            color="#303030",
            ha="left" if label[0] >= 0 else "right",
            va="center",
            zorder=12
        )


# ============================================================
# DENSE, COMPLETE RIGHT-SIDE KEYS
# ============================================================

def draw_full_key(
    ax, entries, columns,
    x0, x1, top, bottom,
    heading, fontsize
):

    count = len(entries)

    ax.text(
        x0, top,
        f"{heading} (n = {count})",
        fontsize=13,
        fontweight="bold",
        va="top",
        color="#202020"
    )

    text_top = top - 0.018

    rows = math.ceil(count / columns)
    col_width = (x1 - x0) / columns
    row_step = (text_top - bottom) / max(rows, 1)

    # Explicit check; never silently omit labels.
    row_height_pt = (
        row_step * FIG_H * 72
    )

    if row_height_pt < fontsize * 1.22:
        raise RuntimeError(
            f"Not enough height for {heading}. "
            f"Available row spacing = {row_height_pt:.1f} pt, "
            f"required = {fontsize * 1.22:.1f} pt."
        )

    for i, (colour, label) in enumerate(entries):
        col = i // rows
        row = i % rows

        x = x0 + col * col_width
        y = text_top - (row + 0.5) * row_step

        ax.plot(
            x + 0.0025,
            y,
            marker="s",
            markersize=3.8,
            color=colour,
            markeredgewidth=0,
            linestyle="None"
        )

        ax.text(
            x + 0.0085,
            y,
            label,
            ha="left",
            va="center",
            fontsize=fontsize,
            color="#252525"
        )


# ============================================================
# FULL FIGURE
# ============================================================

def create_figure(filename, curved):

    mode = "CURVED" if curved else "STRAIGHT"
    print(f"\nCreating large {mode} figure...")

    fig = plt.figure(
        figsize=(FIG_W, FIG_H),
        facecolor="white"
    )

    # Larger spheres in three rows:
    #
    # Row 1: 3 species
    # Row 2: 2 species
    # Row 3: 2 species
    #
    # This provides more physical space per sphere.


    # ========================================================
    # MAXIMUM-SIZE SPHERICAL NETWORK ARRANGEMENT
    #
    # 2 + 2 + 2 + 1
    # Four closely spaced rows.
    # Last species centred.
    # ========================================================

    row_layout = [
        SPECIES[0:2],
        SPECIES[2:4],
        SPECIES[4:6],
        SPECIES[6:7]
    ]

    # Optimised vertical distribution
    row_bottoms = [
        0.775,
        0.545,
        0.315,
        0.085
    ]

    row_height = 0.220

    for row_index, row_species in enumerate(row_layout):

        n = len(row_species)

        total_width = LEFT_END - LEFT_START

        # Two larger network panels per row.
        sphere_width = 0.350
        gap = 0.012

        group_width = (
            n * sphere_width
            + (n - 1) * gap
        )

        first_x = (
            LEFT_START
            + (total_width - group_width) / 2
        )

        for j, species in enumerate(row_species):

            if species not in species_data:
                continue

            x = first_x + j * (sphere_width + gap)

            ax = fig.add_axes([
                x,
                row_bottoms[row_index],
                sphere_width,
                row_height
            ])

            draw_sphere(ax, species, curved)

    # --------------------------------------------------------
    # All right-hand names
    # --------------------------------------------------------

    legend = fig.add_axes([0, 0, 1, 1])

    legend.axis("off")
    legend.set_xlim(0, 1)
    legend.set_ylim(0, 1)

    legend.plot(
        [0.746, 0.746],
        [0.045, 0.980],
        linewidth=0.7,
        color="#CBD5E1"
    )

    rep_entries = [
        (
            REP_COLOR[rep],
            f"{REP_ID[rep]}  {rep}"
        )
        for rep in REP_ORDER
    ]

    draw_full_key(
        legend,
        rep_entries,
        columns=REPLICON_COLUMNS,
        x0=RIGHT_START,
        x1=RIGHT_END,
        top=0.975,
        bottom=0.760,
        heading="Plasmid replicons",
        fontsize=REPLICON_FONT
    )

    legend.plot(
        [RIGHT_START, RIGHT_END],
        [0.747, 0.747],
        linewidth=0.7,
        color="#CBD5E1"
    )

    gene_entries = [
        (
            AMR_COLORS[classify_amr(gene)],
            f"{GENE_ID[gene]}  {gene}"
        )
        for gene in GENE_ORDER
    ]

    draw_full_key(
        legend,
        gene_entries,
        columns=AMR_COLUMNS,
        x0=RIGHT_START,
        x1=RIGHT_END,
        top=0.735,
        bottom=0.058,
        heading="Antimicrobial resistance genes",
        fontsize=AMR_FONT
    )

    # --------------------------------------------------------
    # Bottom information: left network section
    # --------------------------------------------------------

    foot = fig.add_axes([
        0.035, 0.007,
        0.685, 0.057
    ])

    foot.axis("off")
    foot.set_xlim(0, 1)
    foot.set_ylim(0, 1)

    total_edges = sum(
        len(d["pairs"])
        for d in species_data.values()
    )

    total_rep_nodes = sum(
        len(d["reps"])
        for d in species_data.values()
    )

    mode_name = (
        "Curved edges" if curved else "Straight edges"
    )

    foot.text(
        0.0, 0.92,
        f"Seven spherical species networks | {mode_name}",
        fontsize=13,
        fontweight="bold",
        va="top"
    )

    foot.text(
        0.0, 0.70,
        (
            f"{len(master):,} genomes | "
            f"{len(species_data)} species | "
            f"{len(REP_ORDER)} replicon types | "
            f"{total_rep_nodes} replicon nodes | "
            f"{len(GENE_ORDER)} AMR genes | "
            f"{total_edges:,} observed associations"
        ),
        fontsize=9.5,
        va="top"
    )

    for i, species in enumerate(SPECIES):
        if species not in species_data:
            continue

        col = i % 4
        row = i // 4

        x = col * 0.245
        y = 0.49 - row * 0.22

        foot.scatter(
            [x + 0.007],
            [y],
            s=65,
            color=SPECIES_COLORS[species]
        )

        foot.text(
            x + 0.022,
            y,
            (
                f"{short_species(species)} "
                f"(n={species_data[species]['N']:,})"
            ),
            va="center",
            fontsize=9
        )

    foot.text(
        0.0, 0.095,
        (
            "Replicon circles and their connections use "
            "identical replicon-specific colours. "
            "AMR squares indicate resistance class. "
            "Source rings: Human (orange), Animal (blue), "
            "Environment (green)."
        ),
        fontsize=8.5,
        va="center"
    )

    foot.text(
        0.0, 0.005,
        (
            "Edges represent genome-level co-occurrence, "
            "not confirmed physical linkage of AMR genes "
            "to plasmids. Complete identities and "
            "statistics are supplied in Excel."
        ),
        fontsize=8.0,
        va="bottom",
        color="#475569"
    )

    # High-memory step: process each image separately.
    fig.savefig(
        filename,
        dpi=DPI,
        format="jpg",
        facecolor="white",
        pil_kwargs={"quality": 90}
    )

    plt.close(fig)
    del fig
    gc.collect()

    print("Saved:", filename)


# ============================================================
# SAVE EXCEL FIRST
# ============================================================

print("\nSaving workbook...")

rep_key = pd.DataFrame([
    {
        "Replicon_ID": REP_ID[rep],
        "Replicon_Name": rep,
        "Colour_HEX": REP_COLOR[rep],
        "Positive_Genomes": global_rep[rep],
        "Overall_Prevalence_Percent":
            100 * global_rep[rep] / len(master)
    }
    for rep in REP_ORDER
])

gene_key = pd.DataFrame([
    {
        "Gene_ID": GENE_ID[gene],
        "AMR_Gene": gene,
        "AMR_Class": classify_amr(gene),
        "Colour_HEX": AMR_COLORS[classify_amr(gene)],
        "Positive_Genomes": global_gene[gene],
        "Overall_Prevalence_Percent":
            100 * global_gene[gene] / len(master)
    }
    for gene in GENE_ORDER
])

with pd.ExcelWriter(EXCEL, engine="openpyxl") as writer:

    pd.DataFrame(summary_rows).to_excel(
        writer,
        sheet_name="Species_Summary",
        index=False
    )

    rep_key.to_excel(
        writer,
        sheet_name="All_Replicon_Labels",
        index=False
    )

    gene_key.to_excel(
        writer,
        sheet_name="All_AMR_Gene_Labels",
        index=False
    )

    pd.DataFrame(rep_rows).to_excel(
        writer,
        sheet_name="Replicon_Prevalence",
        index=False
    )

    pd.DataFrame(gene_rows).to_excel(
        writer,
        sheet_name="AMR_Gene_Prevalence",
        index=False
    )

    pd.DataFrame(association_rows).to_excel(
        writer,
        sheet_name="All_Associations",
        index=False
    )

    master[[
        "Assembly_Accession",
        "Normalized_Species",
        "Country_Harmonized",
        "Final_Source_Category",
        "BioProject",
        "Verified_Collection_Year",
        "Plasmid_Replicons",
        "AMR_Genes"
    ]].to_excel(
        writer,
        sheet_name="Genome_Source_Data",
        index=False
    )

print("Excel saved:", EXCEL)

# ============================================================
# BOTH NETWORK VERSIONS
# ============================================================

create_figure(CURVED_JPG, curved=True)
create_figure(STRAIGHT_JPG, curved=False)

print("\n" + "=" * 78)
print("COMPLETED")
print("=" * 78)
print("Curved :", CURVED_JPG)
print("Straight:", STRAIGHT_JPG)
print("Excel  :", EXCEL)
print("=" * 78)
