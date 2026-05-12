#!/usr/bin/env python3
"""
combine_datasets.py

Performs a left join between two tabular datasets (CSV or TSV),
allowing enrichment of a base dataset with additional columns
from a secondary dataset.

Features:
- Automatic delimiter detection (CSV/TSV)
- Flexible key matching (same or different column names)
- Collision-safe column handling
- Optional column filtering (--keep-base-cols, --keep-enrich-cols)
- Column inspection mode (--show-cols)

Typical Use Case:
- Merge pipeline outputs with external summaries
  (e.g., genome FASTA summary + protein pipeline results)

Usage:
  python combine_datasets.py \
    -b <base.tsv> \
    -e <enrich.tsv> \
    -k BaseKey:EnrichKey \
    -o output.tsv

Example:
  python combine_datasets.py \
    -b master_summary_report.tsv \
    -e genome_summary.tsv \
    -k Species:basename \
    -o merged.tsv

Author: Adekola Owoyemi (Casola Lab, ECCB, Texas A&M University)
Version: 1.0.0
"""

import argparse
import csv
import os
import sys

__version__ = "1.0.0"
__author__ = "Adekola Owoyemi (Casola Lab, ECCB, Texas A&M University)"

def detect_delimiter(filepath: str) -> str:
    """Robustly detect the delimiter for a text table file."""
    try:
        with open(filepath, 'r', newline='', encoding='utf-8') as f:
            sample = f.read(8192)
    except FileNotFoundError:
        print(f" Warning: file not found during detection: {filepath}. Defaulting to ','", file=sys.stderr)
        return ','
    except Exception as e:
        print(f" Warning: could not read sample from {filepath} ({e}). Defaulting to ','", file=sys.stderr)
        return ','

    # Heuristic: tabs present and commas absent -> TSV
    if '\t' in sample and ',' not in sample:
        print(f" Detected TSV (tab-delimited) for: {os.path.basename(filepath)}", file=sys.stderr)
        return '\t'

    # Try csv.Sniffer
    try:
        dialect = csv.Sniffer().sniff(sample)
        delim = getattr(dialect, 'delimiter', ',')
        if delim in ('\t', ','):
            print(f" Detected delimiter '{delim}' for: {os.path.basename(filepath)}", file=sys.stderr)
            return delim
        else:
            print(f" Sniffer suggested delimiter '{delim}' for: {os.path.basename(filepath)}; using fallback.", file=sys.stderr)
    except Exception:
        pass

    # Fallback by simple counts
    tab_count = sample.count('\t')
    comma_count = sample.count(',')
    if tab_count > comma_count:
        print(f" Fallback: chose '\\t' for: {os.path.basename(filepath)} (tabs={tab_count}, commas={comma_count})", file=sys.stderr)
        return '\t'
    else:
        print(f" Fallback: chose ',' for: {os.path.basename(filepath)} (tabs={tab_count}, commas={comma_count})", file=sys.stderr)
        return ','


def load_enrichment_data(filepath: str, key_column: str, delimiter: str):
    """Load enrichment (secondary) file into a dict keyed by key_column."""
    data_dict = {}
    headers = []
    print(f"--- Loading enrichment data from: {filepath} ---", file=sys.stderr)
    try:
        with open(filepath, 'r', newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            headers = [h.strip() for h in (reader.fieldnames or [])]
            
            if key_column not in headers:
                print(f"Error: Key column '{key_column}' not found in enrichment file.", file=sys.stderr)
                print(f" Available columns in enrichment file: {headers}", file=sys.stderr)
                sys.exit(1)
            
            reader.fieldnames = headers

            count = 0
            for row in reader:
                key_value = row.get(key_column)
                if key_value is not None:
                    key_value = key_value.strip()
                if key_value:
                    data_dict[key_value] = row
                    count += 1
        print(f" Loaded {count} enrichment records (Key: {key_column}).", file=sys.stderr)
        return data_dict, headers
    except FileNotFoundError:
        print(f"Error: Enrichment file not found: {filepath}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while loading enrichment file: {e}", file=sys.stderr)
        sys.exit(1)

def show_columns(filepath, delimiter, label):
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter=delimiter)
            headers = next(reader)
            print(f"\n[{label} FILE COLUMNS]")
            for col in headers:
                print(col.strip())
    except Exception as e:
        print(f"[ERROR] Could not read {filepath}: {e}", file=sys.stderr)
        sys.exit(1)

def main(args):
    base_delim = args.base_delim or detect_delimiter(args.base_file)
    enrich_delim = args.enrich_delim or detect_delimiter(args.enrich_file)

    if args.show_cols:
        show_columns(args.base_file, base_delim, "BASE")
        show_columns(args.enrich_file, enrich_delim, "ENRICH")
        sys.exit(0)

    if not args.show_cols and not args.output_file:
        print("[ERROR] --output-file is required unless --show-cols is used", file=sys.stderr)
        sys.exit(1)
        
    # --- Parse key ---
    if ':' in args.key:
        base_key, enrich_key = args.key.split(':', 1)
        base_key = base_key.strip()
        enrich_key = enrich_key.strip()
        print(f" Mapping keys: Base['{base_key}'] <--> Enrich['{enrich_key}']", file=sys.stderr)
    else:
        base_key = enrich_key = args.key.strip()
        print(f" Using single key for both files: '{base_key}'", file=sys.stderr)

    # --- Handle keep columns ---
    keep_base_cols = None
    keep_enrich_cols = None

    if args.keep_base_cols:
        keep_base_cols = [c.strip() for c in args.keep_base_cols.split(",")]
        print(f" Keeping base columns: {keep_base_cols}", file=sys.stderr)

    if args.keep_enrich_cols:
        keep_enrich_cols = [c.strip() for c in args.keep_enrich_cols.split(",")]
        print(f" Keeping enrichment columns: {keep_enrich_cols}", file=sys.stderr)

    # --- Load enrichment ---
    enrichment_data, enrich_headers = load_enrichment_data(
        args.enrich_file, enrich_key, enrich_delim
    )

    print(f"--- Processing base file: {args.base_file} ---", file=sys.stderr)

    try:
        with open(args.base_file, 'r', newline='', encoding='utf-8') as infile, \
            open(args.output_file, 'w', newline='', encoding='utf-8') as outfile:

            # ✅ robust parsing
            reader = csv.DictReader(infile, delimiter=base_delim, skipinitialspace=True)
            base_headers = [h.strip() for h in (reader.fieldnames or [])]

            # ✅ filter base columns if requested
            if keep_base_cols:
                base_headers = [h for h in base_headers if h in keep_base_cols or h == base_key]

            if base_key not in base_headers:
                print(f"Error: Key column '{base_key}' not found in base file.", file=sys.stderr)
                print(f" Available columns in base file: {base_headers}", file=sys.stderr)
                sys.exit(1)

            # --- Column handling ---
            enrich_col_map = {}
            columns_to_add = []

            for h in enrich_headers:

                if h == enrich_key:
                    continue

                if keep_enrich_cols and h not in keep_enrich_cols:
                    continue

                if h in base_headers:
                    new_name = f"{h}_enrich"
                    enrich_col_map[h] = new_name
                    print(f" Note: Column '{h}' exists in base file. Renaming to '{new_name}'.", file=sys.stderr)
                else:
                    enrich_col_map[h] = h

                columns_to_add.append(enrich_col_map[h])

            final_headers = base_headers + columns_to_add

            # ✅ FIX: match output delimiter to base file
            writer = csv.DictWriter(
                outfile,
                fieldnames=final_headers,
                delimiter=base_delim,
                extrasaction='ignore'
            )
            writer.writeheader()

            processed_count = 0

            for base_row in reader:
                key_value = base_row.get(base_key)
                if key_value is not None:
                    key_value = key_value.strip()

                # ✅ only keep filtered base columns
                if keep_base_cols:
                    final_row = {k: base_row.get(k) for k in base_headers}
                else:
                    final_row = base_row.copy()

                if key_value and key_value in enrichment_data:
                    enrichment_row = enrichment_data[key_value]

                    for orig_col, target_col in enrich_col_map.items():
                        value = enrichment_row.get(orig_col)

                        if value is None and args.fill_value is not None:
                            final_row[target_col] = args.fill_value
                        elif value is not None:
                            final_row[target_col] = value

                else:
                    if args.fill_value is not None:
                        for target_col in columns_to_add:
                            final_row[target_col] = args.fill_value

                writer.writerow(final_row)
                processed_count += 1
                
    except FileNotFoundError:
        print(f"Error: Base file not found: {args.base_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

    print("\n[INFO] Combination complete")
    print(f"[INFO] Processed {processed_count} records.")
    print(f"[INFO] Output saved to: {args.output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Combines two datasets (CSV/TSV) by performing a left join.\n"
            "Supports different key column names via -k BASE:ENRICH."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-b", "--base-file",
        required=True,
        dest="base_file",
        help="The primary file (e.g., genes.csv)",
    )
    parser.add_argument(
        "-e", "--enrich-file",
        required=True,
        dest="enrich_file",
        help="The enrichment file (e.g., proteins.csv)",
    )
    parser.add_argument(
        "-o", "--output-file",
        required=False,
        dest="output_file",
        help="The path for the final output file",
    )
    parser.add_argument(
        "-k", "--key",
        default="Gene_ID",
        dest="key",
        help=(
            "The key column to join on.\n"
            "If same in both files: -k Gene_ID\n"
            "If different: -k BaseColumn:EnrichColumn (e.g. -k Gene_ID:Protein_ID)"
        ),
    )
    parser.add_argument(
        "-f", "--fill-value",
        default=None,
        dest="fill_value",
        help="Value to fill for missing matches (e.g., 'NA', '0')",
    )
    parser.add_argument(
        "--base-delim",
        choices=[",", "\t"],
        default=None,
        dest="base_delim",
        help="Force delimiter for base file (overrides auto-detection).",
    )

    parser.add_argument(
        "--show-cols",
        action="store_true",
        help="Print column names from both files and exit"
    )

    parser.add_argument(
        "--keep-base-cols",
        help="Comma-separated columns to retain from base file"
    )

    parser.add_argument(
        "--keep-enrich-cols",
        help="Comma-separated columns to retain from enrichment file"
    )

    parser.add_argument(
        "--enrich-delim",
        choices=[",", "\t"],
        default=None,
        dest="enrich_delim",
        help="Force delimiter for enrichment file (overrides auto-detection).",
    )


    args = parser.parse_args()
    main(args)
    