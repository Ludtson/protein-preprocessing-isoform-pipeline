#!/usr/bin/env python3
"""
fasta_summary.py

Summarizes FASTA files by reporting:
- Basename (no extension)
- File size (MB)
- Number of sequences

Supports both single-file and directory modes, with optional parallel execution.

Designed as a lightweight companion to the LINGUA preprocessing pipeline.

Usage:
------
# Single file
python fasta_summary.py -i genome.fa

# Directory (sequential)
python fasta_summary.py -d genome_dir/

# Directory (parallel)
python fasta_summary.py -d genome_dir -t 8

# Save output
python fasta_summary.py -d genome_dir -t 8 -o summary.tsv

Output:
-------
basename    size_MB    n_seqs

Author: Adekola Owoyemi (Casola Lab, ECCB, Texas A&M University)
Version: 1.0.0
"""

import argparse
import os
import sys
import subprocess
from concurrent.futures import ProcessPoolExecutor
import time

__version__ = "1.0.0"
__author__ = "Adekola Owoyemi (Casola Lab, ECCB, Texas A&M University)"


def count_sequences(filepath):
    """Fast sequence counting using grep."""
    try:
        result = subprocess.run(
            ["grep", "-c", "^>", filepath],
            capture_output=True,
            text=True
        )
        return int(result.stdout.strip())
    except Exception as e:
        print(f"[WARNING] Could not count sequences in {filepath}: {e}", file=sys.stderr)
        return 0


def summarize_file(filepath):
    """Return (basename, size_MB, n_seqs)."""
    try:
        size_mb = os.path.getsize(filepath) / (1024 * 1024)
    except Exception as e:
        print(f"[WARNING] Could not get size for {filepath}: {e}", file=sys.stderr)
        size_mb = 0

    n_seqs = count_sequences(filepath)

    basename = os.path.basename(filepath)
    basename = os.path.splitext(basename)[0]

    return basename, round(size_mb, 2), n_seqs


def process_directory(directory):
    """Sequential processing."""
    results = []

    for file in os.listdir(directory):
        if file.startswith("."):
            continue

        path = os.path.join(directory, file)

        if os.path.isfile(path):
            results.append(summarize_file(path))

    return results


def process_directory_parallel(directory, threads):
    """Parallel processing."""
    files = [
        os.path.join(directory, f)
        for f in os.listdir(directory)
        if not f.startswith(".") and os.path.isfile(os.path.join(directory, f))
    ]

    results = []

    with ProcessPoolExecutor(max_workers=threads) as executor:
        for res in executor.map(summarize_file, files):
            results.append(res)

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Summarize FASTA files (basename, size, sequence count)."
    )

    parser.add_argument("-i", "--input", help="Single FASTA file")
    parser.add_argument("-d", "--dir", help="Directory containing FASTA files")
    parser.add_argument("-o", "--output", help="Output TSV file")

    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1)"
    )

    args = parser.parse_args()

    if not args.input and not args.dir:
        print("[ERROR] Provide either -i (file) or -d (directory)", file=sys.stderr)
        sys.exit(1)
    
    start_time = time.time()

    results = []

    # --- Single file ---
    if args.input:
        if not os.path.isfile(args.input):
            print(f"[ERROR] File not found: {args.input}", file=sys.stderr)
            sys.exit(1)

        results.append(summarize_file(args.input))

    if args.dir:
        files = [
            f for f in os.listdir(args.dir)
            if not f.startswith(".") and os.path.isfile(os.path.join(args.dir, f))
        ]
        
        print(f"[INFO] Found {len(files)} FASTA files", file=sys.stderr)

    # --- Directory ---
    if args.dir:
        if not os.path.isdir(args.dir):
            print(f"[ERROR] Directory not found: {args.dir}", file=sys.stderr)
            sys.exit(1)

        if args.threads > 1:
            print(f"[INFO] Processing {len(files)} files using {args.threads} threads", file=sys.stderr)
            results.extend(process_directory_parallel(args.dir, args.threads))
        else:
            results.extend(process_directory(args.dir))

    # --- Sort output (nice for reproducibility) ---
    results.sort(key=lambda x: x[0])

    # --- Output ---
    output_lines = ["basename\tsize_MB\tn_seqs"]

    for name, size, nseq in results:
        output_lines.append(f"{name}\t{size}\t{nseq}")

    if args.output:
        try:
            with open(args.output, "w") as out:
                out.write("\n".join(output_lines) + "\n")
        except Exception as e:
            print(f"[ERROR] Could not write output: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print("\n".join(output_lines))
    
    end_time = time.time()
    elapsed = round(end_time - start_time, 2)

    print(f"[INFO] Completed: {len(results)} files in {elapsed} seconds", file=sys.stderr)

if __name__ == "__main__":
    main()
