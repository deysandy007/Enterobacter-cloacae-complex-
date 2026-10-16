
#!/usr/bin/env python3

# ============================================================
# FIGURE 8B — COMPLETE MLST TREEMAP WITH INCOME STRIPS
#
# Every valid ST receives one independent rectangle.
# Every ST rectangle receives a proportional income strip.
#
# Species = background color and block
# ST abundance = rectangle area
# Income = H / M / L+LM / Unclassified strip
#
# No "Other STs" grouping.
# One landscape figure, 600 DPI JPEG.
# Detailed Excel with all ST-income counts.
# ============================================================

from pathlib import Path
import re
import gc

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.colors import to_rgb

import squarify

from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter


# ============================================================
# PATHS
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

OUT = (
    ROOT / "Fig_08B_Combined_MLST_Treemap"
    / "ALL_ST_INCOME_STRIPS_EVERY_TILE"
)

OUT.mkdir(parents=True, exist_ok=True)

JPG = OUT / "Fig08B_ALL_ST_Income_Strips_Every_Tile_600DPI.jpg"

EXCEL = OUT / "Fig08B_ALL_ST_Income_Strips_Data.xlsx"

DPI = 600

# Large landscape format.
FIG_W = 25.0
FIG_H = 17.2

CANVAS_W = 1000.0
CANVAS_H = 720.0

# Income-strip fraction of each tile's height.
# No minimum tile-width or tile-height exclusion.
STRIP_FRACTION = 0.24

# Every ST tile is drawn, even if too small to label.
DRAW_ALL_STS = True


# ============================================================
# SPECIES AND COLORS
# ============================================================

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
    "Enterobacter kobei": "#7851A9",
    "Enterobacter ludwigii": "#A6761D",
    "Enterobacter bugandensis": "#56B4E9"
}

SPECIES_ROWS = [
    ["Enterobacter hormaechei"],
    [
        "Enterobacter cloacae",
        "Enterobacter roggenkampii",
        "Enterobacter asburiae"
    ],
    [
        "Enterobacter kobei",
        "Enterobacter ludwigii",
        "Enterobacter bugandensis"
    ]
]

INCOME_ORDER = [
    "H",
    "M",
    "L+LM",
    "Unclassified"
]

INCOME_COLORS = {
    "H": "#C98A00",
    "M": "#673AB7",
    "L+LM": "#C94F82",
    "Unclassified": "#929DA8"
}

INCOME_LABELS = {
    "H": "High income",
    "M": "Upper-middle income",
    "L+LM": "Lower-middle + low income",
    "Unclassified": "Income unclassified"
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "figure.facecolor": "white",
    "savefig.facecolor": "white"
})


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def clean_st(value):

    if pd.isna(value):
        return ""

    s = str(value).strip().strip("'\"").strip()

    if s.upper().startswith("ST"):
        s = s[2:].strip()

    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".")[0]

    if not re.fullmatch(r"\d+", s):
        return ""

    number = int(s)

    return f"ST{number}" if number > 0 else ""


def country_key(value):

    if pd.isna(value):
        return ""

    return re.sub(
        r"\s+", " ", str(value).strip()
    ).casefold()


def normalize_income(value):

    if pd.isna(value):
        return "Unclassified"

    s = re.sub(
        r"[\s_\-]+",
        "",
        str(value).strip().upper()
    )

    mapping = {
        "H": "H",
        "HIGH": "H",
        "HIGHINCOME": "H",
        "HIGHINCOMECOUNTRIES": "H",

        "M": "M",
        "UM": "M",
        "UMIC": "M",
        "UPPERMIDDLE": "M",
        "UPPERMIDDLEINCOME": "M",
        "UPPERMIDDLEINCOMECOUNTRIES": "M",

        "L+LM": "L+LM",
        "LM+L": "L+LM",
        "LOWERMIDDLE+LOW": "L+LM",
        "LOW+LOWERMIDDLE": "L+LM",
        "LOWERMIDDLEANDLOW": "L+LM",
        "LOWERMIDDLEANDLOWINCOME": "L+LM",
        "LOWERMIDDLE": "L+LM",
        "LOWERMIDDLEINCOME": "L+LM",
        "LOW": "L+LM",
        "LOWINCOME": "L+LM"
    }

    return mapping.get(s, "Unclassified")


def lighten(color, amount):

    rgb = np.array(to_rgb(color))

    return tuple(
        np.clip(
            rgb + (1-rgb)*amount,
            0, 1
        )
    )


def choose_text_color(color):

    r, g, b = to_rgb(color)

    luminance = (
        0.2126*r
        + 0.7152*g
        + 0.0722*b
    )

    return "white" if luminance < 0.47 else "#173148"


def style_excel(writer):

    for ws in writer.book.worksheets:

        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions

        for cell in ws[1]:

            cell.fill = PatternFill(
                "solid",
                fgColor="173352"
            )

            cell.font = Font(
                color="FFFFFF",
                bold=True
            )

        for j, cell in enumerate(ws[1], 1):

            ws.column_dimensions[
                get_column_letter(j)
            ].width = min(
                47,
                max(
                    17,
                    len(str(cell.value)) + 4
                )
            )


# ============================================================
# LOAD MASTER
# ============================================================

print("=" * 100)
print("FIGURE 8B — INCOME STRIPS ON EVERY ST TILE")
print("=" * 100)

for file in [MASTER, INCOME_FILE]:

    if not file.is_file():
        raise FileNotFoundError(file)

required = [
    "Assembly_Accession",
    "Normalized_Species",
    "MLST_Scheme",
    "MLST_ST",
    "Country_Harmonized"
]

df = pd.read_excel(
    MASTER,
    sheet_name="FINAL_MASTER_2000_2026",
    engine="openpyxl",
    usecols=required,
    dtype=str
)

if len(df) != 9157:
    raise ValueError(
        f"Expected 9157 genomes, found {len(df)}"
    )

if (
    df["Assembly_Accession"].isna().any()
    or df["Assembly_Accession"].duplicated().any()
):
    raise ValueError(
        "Missing or duplicate assembly accessions."
    )

df["Species"] = (
    df["Normalized_Species"]
    .fillna("")
    .str.strip()
)

df["Scheme"] = (
    df["MLST_Scheme"]
    .fillna("")
    .str.strip()
    .str.lower()
)

df["ST"] = (
    df["MLST_ST"].map(clean_st)
)

df["Country_Key"] = (
    df["Country_Harmonized"].map(country_key)
)

df["Resolved"] = (
    df["Species"].isin(SPECIES)
)

df["Valid_ST"] = (
    df["Scheme"].eq("ecloacae")
    & df["ST"].ne("")
)

df["Eligible"] = (
    df["Resolved"]
    & df["Valid_ST"]
)


# ============================================================
# LOAD FROZEN COUNTRY-INCOME DATABASE
# ============================================================

income = pd.read_excel(
    INCOME_FILE,
    sheet_name="COUNTRY_INCOME_COUNTS",
    engine="openpyxl",
    dtype=str
)

for required_column in [
    "Country_Harmonized",
    "Income_Group"
]:
    if required_column not in income.columns:
        raise ValueError(
            f"Missing income column: {required_column}"
        )

income["Country_Key"] = (
    income["Country_Harmonized"].map(country_key)
)

income["Income"] = (
    income["Income_Group"].map(normalize_income)
)

income_nonempty = income[
    income["Country_Key"].ne("")
].copy()

conflicts = (
    income_nonempty
    .groupby("Country_Key")["Income"]
    .nunique()
)

if (conflicts > 1).any():

    bad = conflicts[conflicts > 1]

    raise ValueError(
        "Conflicting country-income labels: "
        + ", ".join(bad.index[:15])
    )

income_map = (
    income_nonempty
    .drop_duplicates("Country_Key")
    .set_index("Country_Key")["Income"]
)

df["Income"] = (
    df["Country_Key"]
    .map(income_map)
    .fillna("Unclassified")
)

resolved = df[df["Resolved"]].copy()
typed = df[df["Eligible"]].copy()

if len(resolved) != 8882:
    raise ValueError(
        f"Expected 8882 resolved genomes, got {len(resolved)}"
    )

if typed.empty:
    raise ValueError("No valid ST-assigned genomes")

TOTAL_TYPED = len(typed)

print("Master genomes:", len(df))
print("Resolved genomes:", len(resolved))
print("ST-assigned genomes:", TOTAL_TYPED)

print("\nIncome distribution:")
print(typed["Income"].value_counts().to_string())


# ============================================================
# COUNT EVERY SPECIES-SPECIFIC ST
# ============================================================

species_total = (
    resolved.groupby("Species")
    .size()
    .reindex(SPECIES, fill_value=0)
)

species_typed = (
    typed.groupby("Species")
    .size()
    .reindex(SPECIES, fill_value=0)
)

st_counts = (
    typed.groupby(["Species", "ST"])
    .agg(
        Genome_Count=(
            "Assembly_Accession",
            "nunique"
        )
    )
    .reset_index()
)

income_counts = (
    typed.groupby(
        ["Species", "ST", "Income"]
    )
    .size()
    .unstack("Income", fill_value=0)
    .reindex(
        columns=INCOME_ORDER,
        fill_value=0
    )
    .reset_index()
)

income_counts = income_counts.rename(
    columns={
        "H": "N_H",
        "M": "N_M",
        "L+LM": "N_LLM",
        "Unclassified": "N_Unclassified"
    }
)

all_st = st_counts.merge(
    income_counts,
    on=["Species", "ST"],
    validate="one_to_one"
)

all_st["ST_Number"] = (
    all_st["ST"].str[2:].astype(int)
)

all_st["Species_Total"] = (
    all_st["Species"].map(species_total)
)

all_st["Species_Typed"] = (
    all_st["Species"].map(species_typed)
)

all_st["Percent_Within_Species"] = (
    100 * all_st["Genome_Count"]
    / all_st["Species_Total"]
)

all_st["Percent_Among_Typed"] = (
    100 * all_st["Genome_Count"]
    / all_st["Species_Typed"]
)

count_columns = [
    "N_H", "N_M", "N_LLM", "N_Unclassified"
]

for name, col in [
    ("H", "N_H"),
    ("M", "N_M"),
    ("LLM", "N_LLM"),
    ("Unclassified", "N_Unclassified")
]:
    all_st[f"Percent_{name}"] = (
        100 * all_st[col]
        / all_st["Genome_Count"]
    )

all_st["Income_Breadth"] = (
    (all_st[["N_H", "N_M", "N_LLM"]] > 0)
    .sum(axis=1)
)

all_st["Cross_Income"] = (
    all_st["Income_Breadth"] >= 2
)

all_st = all_st.sort_values(
    ["Species", "Genome_Count", "ST_Number"],
    ascending=[True, False, True]
).reset_index(drop=True)

all_st["Rank"] = (
    all_st.groupby("Species")
    .cumcount() + 1
)

if not (
    all_st[count_columns].sum(axis=1)
    == all_st["Genome_Count"]
).all():
    raise ValueError(
        "ST income counts do not sum to genome counts."
    )

if int(all_st["Genome_Count"].sum()) != TOTAL_TYPED:
    raise ValueError(
        "Abundance does not sum to typed genomes."
    )


# ============================================================
# FIGURE CANVAS
# ============================================================

fig = plt.figure(
    figsize=(FIG_W, FIG_H),
    facecolor="white"
)

AX_LEFT = 0.026
AX_BOTTOM = 0.137
AX_WIDTH = 0.948
AX_HEIGHT = 0.775

ax = fig.add_axes([
    AX_LEFT,
    AX_BOTTOM,
    AX_WIDTH,
    AX_HEIGHT
])

AX_W_IN = FIG_W * AX_WIDTH
AX_H_IN = FIG_H * AX_HEIGHT

ax.set_xlim(0, CANVAS_W)
ax.set_ylim(0, CANVAS_H)
ax.invert_yaxis()
ax.set_aspect("auto")
ax.axis("off")

tile_rows = []
block_rows = []


# ============================================================
# INCOME STRIP — APPLIED TO EVERY TILE
# ============================================================

def draw_every_income_strip(
    x, y, w, h, row
):
    """
    Always draw a proportional income strip.

    No threshold on tile width or tile height.
    No omitted tiles.

    In extremely small tiles, individual segments may
    be subpixel at final JPEG resolution; this is an
    unavoidable consequence of proportional area.
    """

    if w <= 0 or h <= 0:
        raise ValueError("Invalid ST tile geometry")

    strip_h = h * STRIP_FRACTION

    # Position strip flush with top edge.
    # Do not inset it, because very small tiles
    # would otherwise lose their available width.
    strip_x = x
    strip_y = y
    strip_w = w

    counts = [
        int(row.N_H),
        int(row.N_M),
        int(row.N_LLM),
        int(row.N_Unclassified)
    ]

    total = int(row.Genome_Count)

    if total <= 0:
        raise ValueError(
            "Invalid zero-genome ST encountered."
        )

    position = strip_x

    for group, count in zip(
        INCOME_ORDER, counts
    ):

        if count == 0:
            continue

        segment_width = (
            strip_w * count/total
        )

        rect = Rectangle(
            (position, strip_y),
            segment_width,
            strip_h,
            facecolor=INCOME_COLORS[group],
            edgecolor="none",
            linewidth=0,
            antialiased=False,
            zorder=8
        )

        ax.add_patch(rect)

        position += segment_width

    return strip_h


# ============================================================
# FIT LABEL BELOW STRIP
# ============================================================

def draw_label_below_strip(
    patch,
    x, y, w, h,
    strip_h,
    row,
    face
):

    usable_h = h - strip_h

    width_in = w/CANVAS_W * AX_W_IN
    height_in = usable_h/CANVAS_H * AX_H_IN

    if width_in >= 1.0 and height_in >= 0.70:

        label = (
            f"{row.ST}\n"
            f"n={int(row.Genome_Count):,}\n"
            f"{float(row.Percent_Within_Species):.1f}%"
        )

        fs = min(
            12.5,
            width_in*10.5,
            height_in*13
        )

    elif width_in >= 0.60 and height_in >= 0.34:

        label = (
            f"{row.ST}\n"
            f"n={int(row.Genome_Count):,}"
        )

        fs = min(
            10.0,
            width_in*11,
            height_in*16
        )

    elif width_in >= 0.38 and height_in >= 0.20:

        label = row.ST

        fs = min(
            8.0,
            width_in*11.5,
            height_in*19
        )

    else:
        return "Too small to label"

    if fs < 5.7:
        return "Too small to label"

    artist = ax.text(
        x + w/2,
        y + strip_h + usable_h/2,
        label,
        ha="center",
        va="center",
        fontsize=fs,
        fontweight="bold",
        linespacing=1.12,
        color=choose_text_color(face),
        zorder=12
    )

    artist.set_clip_path(patch)

    return "Labeled"


# ============================================================
# COMBINED SEVEN-SPECIES TREEMAP
# ============================================================

current_y = 0.0

for group in SPECIES_ROWS:

    row_total = sum(
        int(species_typed[s])
        for s in group
    )

    if row_total <= 0:
        continue

    row_h = (
        CANVAS_H * row_total / TOTAL_TYPED
    )

    current_x = 0.0

    for species in group:

        n_typed = int(species_typed[species])
        n_all = int(species_total[species])

        if n_typed <= 0:
            continue

        block_w = (
            CANVAS_W * n_typed / row_total
        )

        block_h = row_h

        x = current_x
        y = current_y

        base = SPECIES_COLORS[species]

        # Species block.
        ax.add_patch(
            Rectangle(
                (x+1, y+1),
                block_w-2,
                block_h-2,
                facecolor=base,
                edgecolor="white",
                linewidth=2.3,
                zorder=1
            )
        )

        header_h = min(
            30.0,
            block_h*0.22
        )

        species_pct = (
            100 * n_typed / TOTAL_TYPED
        )

        if block_w >= 340:

            heading = (
                f"{species} "
                f"(n={n_typed:,}; {species_pct:.1f}%)"
            )

            font = 13.5

        else:

            heading = (
                species.replace("Enterobacter ", "E. ")
                + f" (n={n_typed:,})"
            )

            font = (
                10.0 if block_w > 180 else 8.0
            )

        ax.text(
            x+8,
            y+header_h*0.55,
            heading,
            va="center",
            ha="left",
            fontsize=font,
            fontweight="bold",
            fontstyle="italic",
            color="white",
            zorder=15
        )

        # The total tile region is exactly proportional
        # to the number of genomes in this species.
        px = x+6
        py = y+header_h+2

        pw = max(0.1, block_w-12)
        ph = max(0.1, block_h-header_h-8)

        species_st = all_st[
            all_st["Species"].eq(species)
        ].sort_values(
            ["Genome_Count", "ST_Number"],
            ascending=[False, True]
        ).reset_index(drop=True)

        sizes = species_st[
            "Genome_Count"
        ].astype(float).tolist()

        layout = squarify.squarify(
            squarify.normalize_sizes(
                sizes, pw, ph
            ),
            px, py, pw, ph
        )

        if len(layout) != len(species_st):
            raise ValueError(
                "Not all STs received rectangles "
                f"for {species}"
            )

        for i, (geometry, row) in enumerate(
            zip(
                layout,
                species_st.itertuples(index=False)
            )
        ):

            rx = float(geometry["x"])
            ry = float(geometry["y"])
            rw = float(geometry["dx"])
            rh = float(geometry["dy"])

            shade = (
                0.15
                + 0.62*i/max(1, len(species_st)-1)
            )

            face = lighten(base, shade)

            # Thin borders even for cross-income STs,
            # so they do not obscure small color strips.
            if int(row.Income_Breadth) >= 2:
                border = "#24384A"
                linewidth = 0.48
            else:
                border = "white"
                linewidth = 0.25

            if rw < 1.0 or rh < 1.0:
                linewidth = 0.08

            tile_patch = Rectangle(
                (rx, ry),
                rw, rh,
                facecolor=face,
                edgecolor=border,
                linewidth=linewidth,
                antialiased=False,
                zorder=3
            )

            ax.add_patch(tile_patch)

            # Income strip applied UNCONDITIONALLY.
            strip_h = draw_every_income_strip(
                rx, ry, rw, rh, row
            )

            label_status = draw_label_below_strip(
                tile_patch,
                rx, ry, rw, rh,
                strip_h,
                row,
                face
            )

            tile_rows.append({
                "Species": species,
                "ST": row.ST,
                "Genome_Count": int(row.Genome_Count),
                "Percent_Within_Species": float(
                    row.Percent_Within_Species
                ),
                "N_H": int(row.N_H),
                "N_M": int(row.N_M),
                "N_LLM": int(row.N_LLM),
                "N_Unclassified": int(
                    row.N_Unclassified
                ),
                "Percent_H": float(row.Percent_H),
                "Percent_M": float(row.Percent_M),
                "Percent_LLM": float(
                    row.Percent_LLM
                ),
                "Percent_Unclassified": float(
                    row.Percent_Unclassified
                ),
                "Income_Breadth": int(
                    row.Income_Breadth
                ),
                "Income_Strip_Drawn": True,
                "Label_Status": label_status,
                "Tile_X": rx,
                "Tile_Y": ry,
                "Tile_Width": rw,
                "Tile_Height": rh,
                "Strip_Height": strip_h
            })

        block_rows.append({
            "Species": species,
            "Total_Genomes": n_all,
            "ST_Assigned_Genomes": n_typed,
            "Unassigned_ST_Genomes": n_all-n_typed,
            "Unique_STs": len(species_st),
            "Block_X": x,
            "Block_Y": y,
            "Block_Width": block_w,
            "Block_Height": block_h
        })

        current_x += block_w

    current_y += row_h


# ============================================================
# TITLES AND LEGEND
# ============================================================

fig.text(
    0.5, 0.970,
    "Combined MLST abundance and country-income composition",
    ha="center",
    fontsize=22,
    fontweight="bold",
    color="#173352"
)

fig.text(
    0.5, 0.941,
    (
        f"Enterobacter cloacae complex | "
        f"{TOTAL_TYPED:,} ST-assigned genomes"
    ),
    ha="center",
    fontsize=12,
    color="#52687B"
)

handles = [
    Patch(
        facecolor=INCOME_COLORS[g],
        edgecolor="none",
        label=INCOME_LABELS[g]
    )
    for g in INCOME_ORDER
]

handles.append(
    Patch(
        facecolor="white",
        edgecolor="#24384A",
        linewidth=1.2,
        label="ST across 2+ income groups"
    )
)

fig.legend(
    handles=handles,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.099),
    ncol=5,
    frameon=False,
    fontsize=10,
    handlelength=1.6,
    columnspacing=1.5
)

fig.text(
    0.5, 0.071,
    "Every individual ST tile contains a proportional "
    "country-income color strip along its top edge.",
    ha="center",
    fontsize=10.0,
    color="#344B60"
)

fig.text(
    0.5, 0.050,
    "Species and ST areas are proportional to genome "
    "counts; no STs are merged or excluded from the "
    "treemap because of low abundance.",
    ha="center",
    fontsize=9.2,
    color="#52687B"
)

fig.text(
    0.5, 0.030,
    "Very rare STs may be too small for readable "
    "labels or individually visible income segments; "
    "exact values are retained in Excel.",
    ha="center",
    fontsize=8.8,
    color="#65798B"
)

fig.savefig(
    JPG,
    dpi=DPI,
    format="jpeg",
    facecolor="white",
    pil_kwargs={
        "quality": 94,
        "subsampling": 1
    }
)

plt.close(fig)
gc.collect()

print("[SAVED FIGURE]", JPG)


# ============================================================
# VALIDATION AND EXCEL
# ============================================================

tiles = pd.DataFrame(tile_rows)
blocks = pd.DataFrame(block_rows)

if len(tiles) != len(all_st):
    raise ValueError(
        "Every ST must have its own rectangle."
    )

if not tiles["Income_Strip_Drawn"].all():
    raise ValueError(
        "Some ST tiles are missing income strips."
    )

if tiles.duplicated(["Species", "ST"]).any():
    raise ValueError(
        "Duplicate ST tile encountered."
    )

if int(tiles["Genome_Count"].sum()) != TOTAL_TYPED:
    raise ValueError(
        "Tile totals do not match typed genomes."
    )

species_income = (
    typed.groupby(["Species", "Income"])
    .size()
    .unstack("Income", fill_value=0)
    .reindex(
        index=SPECIES,
        columns=INCOME_ORDER,
        fill_value=0
    )
    .reset_index()
)

qc = pd.DataFrame([
    ["Master genomes", len(df)],
    ["Resolved seven-species genomes", len(resolved)],
    ["ST-assigned genomes", TOTAL_TYPED],
    ["Unassigned ST genomes", len(resolved)-TOTAL_TYPED],
    ["Distinct species-ST combinations", len(all_st)],
    ["Individual rectangles created", len(tiles)],
    ["Income strips created", int(tiles["Income_Strip_Drawn"].sum())],
    ["Tiles without strips", 0],
    ["Other ST groups", 0],
    ["Income-strip placement", "Top 24% of every ST rectangle"],
    ["Income-strip metric", "Income genomes / total genomes within ST"],
    ["Species/ST tile area", "ST-assigned genome count"],
    ["Income classification", "Frozen Figure 1A country database"],
    ["Figure", "One combined landscape 600 DPI JPEG"],
    ["Limitation", "Subpixel segments in very rare ST tiles may not be visually distinguishable"],
], columns=["Parameter", "Value"])

genome_export = df[[
    "Assembly_Accession",
    "Species",
    "ST",
    "Scheme",
    "Country_Harmonized",
    "Income",
    "Resolved",
    "Valid_ST",
    "Eligible"
]].copy()

with pd.ExcelWriter(
    EXCEL,
    engine="openpyxl"
) as writer:

    genome_export.to_excel(
        writer,
        sheet_name="ALL_9157_GENOMES",
        index=False
    )

    all_st.to_excel(
        writer,
        sheet_name="ALL_ST_INCOME_COUNTS",
        index=False
    )

    tiles.to_excel(
        writer,
        sheet_name="EVERY_ST_TILE",
        index=False
    )

    species_income.to_excel(
        writer,
        sheet_name="SPECIES_INCOME",
        index=False
    )

    blocks.to_excel(
        writer,
        sheet_name="SPECIES_BLOCKS",
        index=False
    )

    qc.to_excel(
        writer,
        sheet_name="QC_METHODS",
        index=False
    )

    style_excel(writer)

print()
print("=" * 100)
print("FIGURE 8B FINISHED")
print("=" * 100)
print("Total typed genomes:", TOTAL_TYPED)
print("Individual ST rectangles:", len(tiles))
print("Income strips drawn:", len(tiles))
print("STs grouped:", 0)
print("JPEG:", JPG)
print("Excel:", EXCEL)
print("=" * 100)
