
from pathlib import Path
import re
import warnings
import pandas as pd

ROOT = Path(r"D:\Enterobacter")
OFFICIAL = ROOT / "ORIGINAL_ANTIBIOTIC_CONSUMPTION"

ECDC = OFFICIAL / "ECDC_ESACNET_2021_OFFICIAL.xlsx"
EMA = OFFICIAL / "EMA_ESVAC_2010_2022_OFFICIAL.pdf"

ECC = (
    ROOT / "ECC_NPJ_AMR_FIGURES"
    / "08_ANTIBIOTIC_CONSUMPTION"
    / "ECC_ANTIBIOTIC_CONSUMPTION_COUNTRY_YEAR_MASTER.xlsx"
)

OUT = OFFICIAL / "ECC_OVERLAP_AUDIT"
OUT.mkdir(parents=True, exist_ok=True)

CODES = {
    "D11_J01DH_HC": ("Carbapenems", "J01DH", "Hospital"),
    "D16_J01XB_HC": ("Polymyxins", "J01XB", "Hospital"),
    "D1_J01A_AC": ("Tetracyclines", "J01A", "Community"),
    "D8_J01A_HC": ("Tetracyclines", "J01A", "Hospital"),
    "D6_J01M_AC": ("Quinolones", "J01M", "Community"),
    "D14_J01M_HC": ("Quinolones", "J01M", "Hospital"),
    "D3_J01D_AC": ("Other beta-lactams (broad)", "J01D", "Community"),
    "D10_J01D_HC": ("Other beta-lactams (broad)", "J01D", "Hospital"),
}

COUNTRY_FIX = {
    "Czech Republic": "Czechia",
    "UK": "United Kingdom",
    "The Netherlands": "Netherlands",
    "Slovak Republic": "Slovakia",
}

def norm_country(x):
    if pd.isna(x):
        return None
    s = str(x).strip()
    return COUNTRY_FIX.get(s, s)

def year_from_column(x):
    try:
        if isinstance(x, (int, float)):
            n = float(x)
            if n.is_integer() and 2012 <= int(n) <= 2021:
                return int(n)
    except (TypeError, ValueError):
        pass

    s = str(x).strip()
    m = re.fullmatch(r"20(?:1[2-9]|2[01])(?:\.0+)?", s)

    return int(float(s)) if m else None

def extract_ecdc():
    print("\n========== ECDC EXTRACTION ==========")

    rows = []
    layout = []

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        for sheet, (drug_class, atc, sector) in CODES.items():
            raw = pd.read_excel(
                ECDC,
                sheet_name=sheet,
                header=None,
                dtype=object
            )

            heading = str(raw.iat[0, 0])
            headers = raw.iloc[1].tolist()

            year_columns = {
                i: year_from_column(h)
                for i, h in enumerate(headers)
                if year_from_column(h) is not None
            }

            layout.append({
                "Sheet": sheet,
                "Antibiotic_Class": drug_class,
                "ATC": atc,
                "Sector": sector,
                "Years": ";".join(
                    str(y) for y in year_columns.values()
                ),
                "Original_Title": heading
            })

            if not year_columns:
                print("NO YEAR COLUMNS:", sheet)
                continue

            for i in range(2, len(raw)):
                country = norm_country(raw.iat[i, 0])

                if not country:
                    continue

                if country.lower() in [
                    "total",
                    "eu/eea",
                    "eu/eea average"
                ]:
                    continue

                if country.lower().startswith(
                    ("note", "source", "table")
                ):
                    continue

                for j, year in year_columns.items():
                    value = pd.to_numeric(
                        raw.iat[i, j],
                        errors="coerce"
                    )

                    if pd.isna(value) or value < 0:
                        continue

                    rows.append({
                        "Country": country,
                        "Year": year,
                        "Sector": sector,
                        "Antibiotic_Class": drug_class,
                        "ATC": atc,
                        "Consumption_DDD_1000_Day": float(value),
                        "Original_Sheet": sheet
                    })

            print(
                sheet,
                "| Class:", drug_class,
                "| Years:", sorted(set(year_columns.values()))
            )

    result = pd.DataFrame(rows)

    result.to_csv(
        OUT / "ECDC_ACTUAL_CLASS_USE_2012_2021.csv",
        index=False
    )

    pd.DataFrame(layout).to_csv(
        OUT / "ECDC_WORKSHEET_AUDIT.csv",
        index=False
    )

    print("ECDC RECORDS:", len(result))
    print("ECDC COUNTRIES:", result["Country"].nunique())

    return result

def extract_ecc_overlap(ecdc):
    print("\n========== ECC GENOME MATCHING ==========")

    ecc = pd.read_excel(
        ECC,
        sheet_name="ECC_9157_EXACT_MATCH",
        usecols=[
            "Assembly_Accession",
            "ECC_Country",
            "Year",
            "Source",
            "Country_ISO3"
        ]
    )

    ecc["Year"] = pd.to_numeric(
        ecc["Year"],
        errors="coerce"
    )

    ecc["Country"] = ecc["ECC_Country"].map(norm_country)

    ecc = ecc[
        ecc["Year"].between(2012, 2021)
    ].copy()

    ecc["Year"] = ecc["Year"].astype(int)

    denominator = (
        ecc.groupby(
            ["Country", "Year", "Source"],
            dropna=False
        )["Assembly_Accession"]
        .nunique()
        .rename("ECC_Genomes")
        .reset_index()
    )

    human = denominator[
        denominator["Source"] == "Human"
    ].copy()

    matched = human.merge(
        ecdc,
        on=["Country", "Year"],
        how="inner"
    )

    matched.to_csv(
        OUT / "HUMAN_ECDC_ECC_COUNTRY_YEAR_MATCH.csv",
        index=False
    )

    denominator.to_csv(
        OUT / "ECC_COUNTRY_YEAR_SOURCE_DENOMINATORS.csv",
        index=False
    )

    unique_pairs = matched[
        ["Country", "Year"]
    ].drop_duplicates()

    print("MATCHED RECORDS:", len(matched))
    print("MATCHED COUNTRY-YEAR PAIRS:", len(unique_pairs))
    print("MATCHED COUNTRIES:", matched["Country"].nunique())

    return denominator, matched

def inspect_ema_pdf():
    import pdfplumber

    print("\n========== EMA PDF AUDIT ==========")

    keywords = [
        "tetracyclines",
        "polymyxins",
        "fluoroquinolones",
        "third- and fourth-generation cephalosporins",
        "mg/pcu",
        "2010–2022"
    ]

    pages = []

    with pdfplumber.open(EMA) as pdf:
        total_pages = len(pdf.pages)

        for page_number, page in enumerate(pdf.pages, 1):
            text = page.extract_text() or ""

            found = [
                k for k in keywords
                if k.lower() in text.lower()
            ]

            if len(found) >= 2:
                pages.append({
                    "PDF_Page": page_number,
                    "Matched_Keywords": "; ".join(found),
                    "Preview": re.sub(
                        r"\s+",
                        " ",
                        text[:230]
                    )
                })

    pd.DataFrame(
        pages,
        columns=["PDF_Page", "Matched_Keywords", "Preview"]
    ).to_csv(
        OUT / "EMA_CANDIDATE_PAGES_FOR_VERIFICATION.csv",
        index=False
    )

    print("EMA TOTAL PAGES:", total_pages)
    print("EMA CANDIDATE PAGES:", len(pages))

    return total_pages, len(pages)

def main():
    for f in [ECDC, EMA, ECC]:
        if not f.exists():
            raise SystemExit(
                "MISSING REQUIRED FILE: " + str(f)
            )

    ecdc = extract_ecdc()
    denominators, matched = extract_ecc_overlap(ecdc)
    ema_pages, candidates = inspect_ema_pdf()

    summary = [
        "ECC ANTIBIOTIC CLASS CONSUMPTION AUDIT",
        "=" * 60,
        f"ECDC records: {len(ecdc)}",
        f"ECDC countries: {ecdc['Country'].nunique()}",
        f"ECDC years: {sorted(ecdc['Year'].unique().tolist())}",
        f"ECC source-specific country-year groups: {len(denominators)}",
        f"Human ECC / ECDC matched records: {len(matched)}",
        "Independent matched country-year pairs: "
        + str(len(matched[['Country', 'Year']].drop_duplicates())),
        f"Matched countries: {matched['Country'].nunique()}",
        f"EMA PDF pages: {ema_pages}",
        f"EMA candidate pages: {candidates}",
        "",
        "IMPORTANT:",
        "EMA numerical tables have NOT yet been extracted.",
        "This workflow assesses overlap only.",
        "No AMR-consumption correlation has been performed.",
        "Broad J01D is not specific third-generation cephalosporin use.",
        "Missing values are not treated as zero."
    ]

    (OUT / "AUDIT_SUMMARY.txt").write_text(
        "\n".join(summary),
        encoding="utf-8"
    )

    print("\n" + "\n".join(summary))
    print("\nOUTPUT:", OUT)

if __name__ == "__main__":
    main()
