
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path(r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES")

ECDC = (
    Path(r"D:\Enterobacter\ORIGINAL_ANTIBIOTIC_CONSUMPTION")
    / "ECDC_ESACNET_2021_OFFICIAL.xlsx"
)

ECC = (
    ROOT / "08_ANTIBIOTIC_CONSUMPTION"
    / "ECC_ANTIBIOTIC_CONSUMPTION_COUNTRY_YEAR_MASTER.xlsx"
)

OUT = ROOT / "FINAL_HUMAN_CARBAPENEM_CORRELATION"
OUT.mkdir(parents=True, exist_ok=True)

GENES = ["blaNDM", "blaIMP", "blaKPC", "blaVIM"]
MIN_GENOMES = 5
MIN_COUNTRIES = 5
RANDOM_SEED = 42
PERMUTATIONS = 10000

COUNTRY_ALIASES = {
    "czech republic": "czechia",
    "the netherlands": "netherlands",
    "slovak republic": "slovakia",
    "uk": "united kingdom",
}

def clean_country(x):
    if pd.isna(x):
        return None
    s = str(x).strip().lower()
    return COUNTRY_ALIASES.get(s, s)

def extract_ecdc():
    print("\n[1] Reading actual ECDC carbapenem use")

    raw = pd.read_excel(
        ECDC,
        sheet_name="D11_J01DH_HC",
        header=None
    )

    years = {}
    for j, val in enumerate(raw.iloc[1]):
        try:
            y = int(float(val))
            if 2012 <= y <= 2021:
                years[j] = y
        except (TypeError, ValueError):
            pass

    records = []

    for i in range(2, len(raw)):
        country = clean_country(raw.iat[i, 0])
        if not country or country.startswith(
            ("note", "source", "table", "eu/eea")
        ):
            continue

        for j, year in years.items():
            value = pd.to_numeric(
                raw.iat[i, j], errors="coerce"
            )
            if pd.notna(value) and value >= 0:
                records.append({
                    "Country": country,
                    "Year": year,
                    "Carbapenem_DDD_1000_Day": float(value)
                })

    result = pd.DataFrame(records)

    result = result.drop_duplicates(
        ["Country", "Year"]
    )

    print("ECDC countries:", result.Country.nunique())
    print("ECDC measurements:", len(result))

    return result

def prepare_ecc():
    print("\n[2] Preparing human ECC carbapenemase prevalence")

    df = pd.read_excel(
        ECC, sheet_name="COUNTRY_YEAR_GENE_LINK"
    )

    df = df[
        df["Source"].astype(str).str.strip().str.lower()
        == "human"
    ].copy()

    df["Country"] = df["ECC_Country"].map(clean_country)
    df["Year"] = pd.to_numeric(df["Year"], errors="coerce")
    df["Genome_Denominator"] = pd.to_numeric(
        df["Genome_Denominator"], errors="coerce"
    )

    df = df[
        df["Year"].between(2012, 2021)
        & df["Country"].notna()
        & (df["Genome_Denominator"] > 0)
    ].copy()

    df["Year"] = df["Year"].astype(int)

    for gene in GENES:
        col = gene + "_Positive"
        df[col] = pd.to_numeric(
            df[col], errors="coerce"
        )

    cols = (
        ["Genome_Denominator"]
        + [g + "_Positive" for g in GENES]
    )

    # Sum the species-stratified genome counts per
    # country and year. These are not independent
    # statistical replicates.
    grouped = df.groupby(
        ["Country", "Year"], as_index=False
    )[cols].sum(min_count=1)

    for gene in GENES:
        grouped[gene + "_Prevalence"] = (
            100 * grouped[gene + "_Positive"]
            / grouped["Genome_Denominator"]
        )

    return grouped

def permutation_pvalue(x, y):
    """
    Two-sided permutation test on country-level observations.
    """
    rng = np.random.default_rng(RANDOM_SEED)

    observed = spearmanr(x, y).statistic

    if not np.isfinite(observed):
        return np.nan

    extreme = 0

    for _ in range(PERMUTATIONS):
        permuted = rng.permutation(y)
        r = spearmanr(x, permuted).statistic

        if np.isfinite(r) and abs(r) >= abs(observed):
            extreme += 1

    return (extreme + 1) / (PERMUTATIONS + 1)

def main():
    if not ECDC.exists():
        raise FileNotFoundError(ECDC)

    if not ECC.exists():
        raise FileNotFoundError(ECC)

    consumption = extract_ecdc()
    genomes = prepare_ecc()

    print("\n[3] Matching consumption and ECC")

    matched = genomes.merge(
        consumption,
        on=["Country", "Year"],
        how="inner",
        validate="one_to_one"
    )

    matched = matched[
        matched["Genome_Denominator"] >= MIN_GENOMES
    ].copy()

    print("Eligible matched country-years:", len(matched))
    print("Eligible matched countries:", matched.Country.nunique())

    # Country is the independent unit in the primary analysis.
    # Consumption is averaged over matched years.
    # Gene prevalence is calculated from pooled counts.
    #
    # A country with many years must not be treated as
    # many statistically independent countries.

    country_rows = []

    for country, sub in matched.groupby("Country"):
        total = sub["Genome_Denominator"].sum()

        record = {
            "Country": country,
            "Matched_Years": len(sub),
            "Total_Human_ECC_Genomes": int(total),
            "Mean_Carbapenem_Use": sub[
                "Carbapenem_DDD_1000_Day"
            ].mean()
        }

        for gene in GENES:
            pos = sub[gene + "_Positive"].sum()
            record[gene + "_Positive"] = pos
            record[gene + "_Prevalence"] = (
                100 * pos / total
            )

        country_rows.append(record)

    country_df = pd.DataFrame(country_rows)

    # Combined presence requires genome-level co-occurrence
    # information to avoid double counting.
    # Therefore, do not sum gene-positive counts to invent
    # a combined carbapenemase prevalence.

    print("\n[4] Running country-level Spearman correlations")

    results = []

    for gene in GENES:
        ycol = gene + "_Prevalence"

        sub = country_df.dropna(
            subset=["Mean_Carbapenem_Use", ycol]
        ).copy()

        n = len(sub)

        result = {
            "Gene": gene,
            "N_Countries": n,
            "N_Country_Years": int(
                sub["Matched_Years"].sum()
            ) if n else 0,
            "Spearman_Rho": np.nan,
            "P_Value": np.nan,
            "Permutation_P": np.nan,
            "Status": "NOT_TESTABLE"
        }

        if (
            n >= MIN_COUNTRIES
            and sub[ycol].nunique() >= 2
            and sub["Mean_Carbapenem_Use"].nunique() >= 2
        ):
            x = sub["Mean_Carbapenem_Use"].to_numpy()
            y = sub[ycol].to_numpy()

            rho, p = spearmanr(x, y)
            perm_p = permutation_pvalue(x, y)

            result.update({
                "Spearman_Rho": rho,
                "P_Value": p,
                "Permutation_P": perm_p,
                "Status": "TESTED"
            })

        results.append(result)

    stats = pd.DataFrame(results)
    stats["FDR_P"] = np.nan

    valid = stats["Permutation_P"].notna()

    if valid.any():
        stats.loc[valid, "FDR_P"] = multipletests(
            stats.loc[valid, "Permutation_P"],
            method="fdr_bh"
        )[1]

    stats["Significance"] = "NT"
    stats.loc[valid, "Significance"] = "ns"
    stats.loc[
        stats["FDR_P"] < 0.05, "Significance"
    ] = "*"
    stats.loc[
        stats["FDR_P"] < 0.01, "Significance"
    ] = "**"
    stats.loc[
        stats["FDR_P"] < 0.001, "Significance"
    ] = "***"

    stats["Decision"] = np.where(
        stats["FDR_P"] < 0.05,
        "STATISTICAL_ASSOCIATION",
        "NO_FDR_SIGNIFICANT_ASSOCIATION"
    )

    stats.loc[
        stats["Status"] == "NOT_TESTABLE", "Decision"
    ] = "NOT_TESTABLE"

    print("\n" + "=" * 90)
    print("FINAL RESULTS")
    print("=" * 90)
    print(stats.to_string(index=False))

    print("\n" + "=" * 90)
    print("MANUSCRIPT DECISION")
    print("=" * 90)

    significant = stats[
        stats["FDR_P"] < 0.05
    ]

    if significant.empty:
        decision = (
            "DROP ANTIBIOTIC CONSUMPTION CORRELATIONS "
            "FROM MAIN MANUSCRIPT: no tested association "
            "survived FDR correction."
        )
    else:
        decision = (
            "REVIEW SIGNIFICANT EXPLORATORY ASSOCIATIONS: "
            + ", ".join(significant["Gene"])
            + ". Association does not establish causality."
        )

    print(decision)

    stats.to_csv(
        OUT / "FINAL_CORRELATION_RESULTS.csv",
        index=False
    )

    matched.to_csv(
        OUT / "MATCHED_COUNTRY_YEAR_DATA.csv",
        index=False
    )

    country_df.to_csv(
        OUT / "COUNTRY_LEVEL_DATA.csv",
        index=False
    )

    (OUT / "FINAL_DECISION.txt").write_text(
        "Human hospital carbapenem consumption versus "
        "ECC carbapenemase gene prevalence\n"
        "ECDC 2012-2021\n"
        "Country is primary statistical unit.\n"
        "Gene prevalence pooled across matched eligible years.\n"
        "Each country-year required at least 5 ECC genomes.\n"
        "Spearman correlation with permutation P-values.\n"
        "Benjamini-Hochberg FDR correction.\n"
        "Exploratory ecological association only.\n\n"
        + stats.to_string(index=False)
        + "\n\n"
        + decision,
        encoding="utf-8"
    )

    print("\nResults saved to:", OUT)

if __name__ == "__main__":
    main()
