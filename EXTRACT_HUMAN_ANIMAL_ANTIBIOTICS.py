
from pathlib import Path
import re
import warnings
import pandas as pd
import pdfplumber

warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path(r"D:\Enterobacter")
BASE = ROOT / "ORIGINAL_ANTIBIOTIC_CONSUMPTION"
OUT = BASE / "HUMAN_ANIMAL_CLASS_SPECIFIC"
OUT.mkdir(parents=True, exist_ok=True)

ECDC = BASE / "ECDC_ESACNET_2021_OFFICIAL.xlsx"
EMA = BASE / "EMA_ESVAC_2010_2022_OFFICIAL.pdf"

ECC = (
    ROOT / "ECC_NPJ_AMR_FIGURES"
    / "08_ANTIBIOTIC_CONSUMPTION"
    / "ECC_ANTIBIOTIC_CONSUMPTION_COUNTRY_YEAR_MASTER.xlsx"
)

ECDC_SHEETS = {
    "D11_J01DH_HC": ("Carbapenems", "J01DH", "Hospital"),
    "D16_J01XB_HC": ("Polymyxins", "J01XB", "Hospital"),
    "D1_J01A_AC": ("Tetracyclines", "J01A", "Community"),
    "D8_J01A_HC": ("Tetracyclines", "J01A", "Hospital"),
    "D6_J01M_AC": ("Quinolones", "J01M", "Community"),
    "D14_J01M_HC": ("Quinolones", "J01M", "Hospital"),
    "D3_J01D_AC": ("Other beta-lactams", "J01D", "Community"),
    "D10_J01D_HC": ("Other beta-lactams", "J01D", "Hospital"),
}

COUNTRY_FIX = {
    "Czech Republic": "Czechia",
    "The Netherlands": "Netherlands",
    "UK": "United Kingdom",
    "Slovak Republic": "Slovakia",
}

GENES = [
    "blaNDM", "blaIMP", "blaKPC", "blaVIM",
    "blaOXA48like_selected", "mcr_any", "mcr10", "tetX4"
]

def country_clean(x):
    if pd.isna(x):
        return None
    s = str(x).strip()
    return COUNTRY_FIX.get(s, s)

def get_year(x):
    try:
        y = int(float(x))
        return y if 2012 <= y <= 2021 else None
    except (ValueError, TypeError):
        return None

def extract_human():
    print("\n[1] EXTRACTING ECDC HUMAN ANTIBIOTICS")

    records = []

    for sheet, (drug, atc, sector) in ECDC_SHEETS.items():
        raw = pd.read_excel(ECDC, sheet_name=sheet, header=None)

        headers = raw.iloc[1].tolist()
        years = {
            j: get_year(h)
            for j, h in enumerate(headers)
            if get_year(h) is not None
        }

        count = 0

        for i in range(2, len(raw)):
            country = country_clean(raw.iat[i, 0])

            if country is None:
                continue

            if (
                country.startswith(("Note", "Source", "Table"))
                or country in ["EU/EEA", "EU/EEA average", "Total"]
            ):
                continue

            for j, year in years.items():
                val = pd.to_numeric(raw.iat[i, j], errors="coerce")

                if pd.isna(val) or val < 0:
                    continue

                records.append({
                    "Country": country,
                    "Year": year,
                    "Sector": "Human",
                    "Care_Setting": sector,
                    "Antibiotic_Class": drug,
                    "ATC_Code": atc,
                    "Consumption": float(val),
                    "Unit": "DDD/1000 inhabitants/day",
                    "Source_Table": sheet
                })
                count += 1

        print(f"  {sheet}: {count} measurements")

    data = pd.DataFrame(records)

    data.to_csv(
        OUT / "HUMAN_ANTIBIOTIC_CLASS_2012_2021.csv",
        index=False
    )

    print("Human records:", len(data))
    print("Human countries:", data.Country.nunique())

    return data

def extract_ema_source_tables():
    print("\n[2] EXTRACTING ORIGINAL EMA VETERINARY TABLES")

    keywords = [
        "tetracyclines",
        "polymyxins",
        "fluoroquinolones",
        "cephalosporins",
        "aminoglycosides"
    ]

    audit = []
    raw_tables = []

    with pdfplumber.open(EMA) as pdf:
        print("EMA PDF pages:", len(pdf.pages))

        for number, page in enumerate(pdf.pages, start=1):
            text = page.extract_text() or ""
            lower = text.lower()
            matches = [x for x in keywords if x in lower]

            if not matches:
                continue

            tables = page.extract_tables() or []

            audit.append({
                "PDF_Page": number,
                "Antibiotic_Classes_Found": "; ".join(matches),
                "Extracted_Tables": len(tables),
                "Page_Contains_mg_PCU": "mg/pcu" in lower,
                "Text_Preview": re.sub(r"\s+", " ", text[:350])
            })

            for table_number, table in enumerate(tables, start=1):
                for row_number, row in enumerate(table, start=1):
                    if row is None:
                        continue

                    record = {
                        "PDF_Page": number,
                        "Table_Number": table_number,
                        "Row_Number": row_number
                    }

                    for col_number, value in enumerate(row, start=1):
                        record[f"Column_{col_number}"] = value

                    raw_tables.append(record)

    audit_df = pd.DataFrame(audit)
    table_df = pd.DataFrame(raw_tables)

    audit_df.to_csv(
        OUT / "EMA_VETERINARY_TABLE_PAGE_AUDIT.csv",
        index=False
    )

    table_df.to_csv(
        OUT / "EMA_ORIGINAL_EXTRACTED_TABLE_ROWS.csv",
        index=False
    )

    print("Candidate EMA pages:", len(audit_df))
    print("Raw EMA table rows:", len(table_df))

    print(
        "EMA numerical tables remain unverified; "
        "do not treat raw PDF cells as validated measurements."
    )

    return audit_df, table_df

def prepare_ecc():
    print("\n[3] EXTRACTING ECC COUNTRY-YEAR GENE COUNTS")

    df = pd.read_excel(
        ECC,
        sheet_name="COUNTRY_YEAR_GENE_LINK"
    )

    df["Country"] = df["ECC_Country"].map(country_clean)
    df["Year"] = pd.to_numeric(df["Year"], errors="coerce")
    df["Genome_Denominator"] = pd.to_numeric(
        df["Genome_Denominator"], errors="coerce"
    )

    df = df[
        df["Year"].between(2012, 2022)
        & df["Country"].notna()
        & (df["Genome_Denominator"] > 0)
        & df["Source"].isin(["Human", "Animal", "Environment"])
    ].copy()

    df["Year"] = df["Year"].astype(int)

    positive_cols = [g + "_Positive" for g in GENES]

    for col in positive_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    grouped = df.groupby(
        ["Country", "Year", "Source"],
        as_index=False
    )[["Genome_Denominator"] + positive_cols].sum(min_count=1)

    for gene in GENES:
        grouped[gene + "_Prevalence_Pct"] = (
            grouped[gene + "_Positive"]
            / grouped["Genome_Denominator"] * 100
        )

    grouped.to_csv(
        OUT / "ECC_COUNTRY_YEAR_ONEHEALTH_GENE_PREVALENCE.csv",
        index=False
    )

    return grouped

def match_human(human, ecc):
    print("\n[4] MATCHING HUMAN CONSUMPTION TO HUMAN ECC")

    h = ecc[ecc["Source"] == "Human"].copy()

    matched = human.merge(
        h,
        on=["Country", "Year"],
        how="inner",
        validate="many_to_one"
    )

    matched["Meets_Minimum_5_Genomes"] = (
        matched["Genome_Denominator"] >= 5
    )

    matched.to_csv(
        OUT / "HUMAN_CLASS_USE_ECC_MATCHED.csv",
        index=False
    )

    print("Matched consumption rows:", len(matched))
    print("Matched countries:", matched.Country.nunique())
    print(
        "Independent country-year pairs:",
        len(matched[["Country", "Year"]].drop_duplicates())
    )
    print(
        "Country-year pairs with >=5 genomes:",
        len(
            matched[
                matched["Meets_Minimum_5_Genomes"]
            ][["Country", "Year"]].drop_duplicates()
        )
    )

    return matched

def animal_coverage(ecc):
    print("\n[5] CHECKING ANIMAL ECC COVERAGE")

    animal = ecc[ecc["Source"] == "Animal"].copy()

    animal.to_csv(
        OUT / "ANIMAL_ECC_COUNTRY_YEAR_FOR_EMA_MATCHING.csv",
        index=False
    )

    print("Animal countries:", animal.Country.nunique())
    print("Animal country-year pairs:", len(animal))
    print(
        "Animal pairs with >=5 genomes:",
        (animal.Genome_Denominator >= 5).sum()
    )

    return animal

def main():
    for path in [ECDC, EMA, ECC]:
        if not path.is_file():
            raise FileNotFoundError(str(path))

    human = extract_human()
    ema_pages, ema_tables = extract_ema_source_tables()
    ecc = prepare_ecc()
    matched = match_human(human, ecc)
    animal = animal_coverage(ecc)

    lines = [
        "ECC ONE HEALTH CLASS-SPECIFIC ANTIBIOTIC DATA",
        "=" * 65,
        "Human data: ECDC 2012-2021",
        "Animal source: EMA ESVAC 2010-2022",
        "",
        f"ECDC human records: {len(human)}",
        f"ECDC countries: {human.Country.nunique()}",
        f"Human consumption/ECC matched rows: {len(matched)}",
        f"Human matched countries: {matched.Country.nunique()}",
        "Human matched country-year pairs: "
        + str(len(matched[["Country", "Year"]].drop_duplicates())),
        "Human matched country-year pairs with >=5 genomes: "
        + str(len(
            matched.loc[
                matched.Meets_Minimum_5_Genomes,
                ["Country", "Year"]
            ].drop_duplicates()
        )),
        "",
        f"Animal ECC countries: {animal.Country.nunique()}",
        f"Animal ECC country-year pairs: {len(animal)}",
        "Animal country-year pairs with >=5 genomes: "
        + str((animal.Genome_Denominator >= 5).sum()),
        f"EMA candidate pages: {len(ema_pages)}",
        f"EMA raw extracted table rows: {len(ema_tables)}",
        "",
        "STATUS",
        "Human consumption and gene prevalence merged.",
        "Animal country-year AMR prevalence prepared.",
        "EMA raw PDF tables extracted for verification.",
        "Longitudinal animal class values NOT YET VALIDATED.",
        "No antibiotic-AMR correlation computed.",
        "No missing consumption values replaced with zero.",
        "ECDC J01D is a broad class, not specific ESBL exposure.",
        "Veterinary sales are not direct antibiotic administration."
    ]

    report = "\n".join(lines)

    (OUT / "ONE_HEALTH_MATCHING_SUMMARY.txt").write_text(
        report, encoding="utf-8"
    )

    print("\n" + report)
    print("\nOUTPUT DIRECTORY:")
    print(OUT)

if __name__ == "__main__":
    main()
