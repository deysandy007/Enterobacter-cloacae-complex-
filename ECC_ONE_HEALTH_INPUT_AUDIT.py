
from pathlib import Path
import csv
from collections import Counter
from datetime import datetime

ROOT = Path(r"D:\Enterobacter")
OUTPUT = Path(r"D:\Enterobacter\ECC_NPJ_AMR_FIGURES")
OUTPUT.mkdir(parents=True, exist_ok=True)

LOG = OUTPUT / "ECC_ONE_HEALTH_INPUT_AUDIT.txt"

EXTENSIONS = {".tsv", ".csv", ".xlsx", ".xls"}

KEYWORDS = [
    "metadata", "master", "amrfinder", "amr",
    "consumption", "antibiotic", "plasmid",
    "country", "one_health", "onehealth",
    "filtered", "genome"
]

def print_both(message="", fh=None):
    print(message)
    if fh:
        fh.write(str(message) + "\n")

def inspect_header(path):
    suffix = path.suffix.lower()

    try:
        if suffix in {".csv", ".tsv"}:
            delimiter = "\t" if suffix == ".tsv" else ","
            with open(path, "r", encoding="utf-8-sig",
                      errors="replace", newline="") as f:
                reader = csv.reader(f, delimiter=delimiter)
                return next(reader, [])

        if suffix in {".xlsx", ".xls"}:
            import pandas as pd
            df = pd.read_excel(path, nrows=0)
            return list(df.columns)

    except Exception as e:
        return [f"ERROR: {e}"]

    return []

def main():
    with open(LOG, "w", encoding="utf-8") as log:

        print_both("=" * 80, log)
        print_both("ECC ONE HEALTH AMR INPUT DATABASE AUDIT", log)
        print_both("=" * 80, log)
        print_both(f"Date: {datetime.now()}", log)
        print_both(f"Root directory: {ROOT}", log)

        if not ROOT.exists():
            print_both("ERROR: ROOT DIRECTORY NOT FOUND", log)
            return

        all_files = []
        for path in ROOT.rglob("*"):
            if path.is_file() and path.suffix.lower() in EXTENSIONS:
                all_files.append(path)

        print_both(f"\nTotal candidate data files: {len(all_files)}", log)

        groups = {
            "GENOME METADATA": [
                "metadata", "master", "filtered", "genome"
            ],
            "AMRFINDERPLUS RESULTS": [
                "amrfinder", "amr"
            ],
            "ANTIBIOTIC CONSUMPTION": [
                "consumption", "antibiotic", "aware"
            ],
            "ONE HEALTH / COUNTRY": [
                "one_health", "onehealth", "country", "source"
            ]
        }

        for group_name, words in groups.items():

            print_both("\n" + "=" * 80, log)
            print_both(group_name, log)
            print_both("=" * 80, log)

            matched = [
                f for f in all_files
                if any(
                    word in f.name.lower()
                    or word in f.parent.name.lower()
                    for word in words
                )
            ]

            print_both(f"Matching files: {len(matched)}", log)

            for path in sorted(matched)[:50]:
                try:
                    size_mb = path.stat().st_size / 1048576
                    print_both(
                        f"{path}  |  {size_mb:.2f} MB", log
                    )
                except OSError:
                    pass

        print_both("\n" + "=" * 80, log)
        print_both("IMPORTANT MASTER FILE COLUMN HEADERS", log)
        print_both("=" * 80, log)

        candidates = [
            f for f in all_files
            if any(k in f.name.lower() for k in [
                "master", "metadata", "consumption",
                "country", "one_health"
            ])
        ]

        candidates.sort(
            key=lambda p: (
                "master" not in p.name.lower(),
                str(p).lower()
            )
        )

        for path in candidates[:30]:
            print_both(f"\nFILE: {path}", log)
            columns = inspect_header(path)

            for col in columns[:100]:
                print_both(f"  {col}", log)

        print_both("\n" + "=" * 80, log)
        print_both("AUDIT COMPLETED", log)
        print_both("=" * 80, log)

    print(f"\nAudit report saved to:\n{LOG}")

if __name__ == "__main__":
    main()
