#!/usr/bin/env python3
"""
fasta_header_tool.py

A flexible FASTA header processing tool.

Features:
- Simplifies headers to primary identifier (first token)
- Optional addition of species name as suffix
- Optional removal of species suffix
- Supports single file and batch (directory) processing
- Optional recursive directory traversal
- Optional sequence linearization
- Supports gzipped FASTA files

Usage:

1. Clean headers:
   python fasta_header_tool.py -i input.faa -o output.faa

2. Add species manually:
   python fasta_header_tool.py -i input.faa -o output.faa --add-species "Arabidopsis thaliana"

3. Add species from map:
   python fasta_header_tool.py -i input_dir -o output_dir --add-species-map species.map

4. Remove species suffix:
   python fasta_header_tool.py -i input.faa -o output.faa --remove-species

5. Process directory recursively:
   python fasta_header_tool.py -i input_dir -o output_dir --recursive

Output:
- Writes FASTA files with modified headers
- Sequence content unchanged (unless --linearize)

Notes:
- Species names appended as "__Species_Name"
- Map fallback = filename basename

Author: Adekola Owoyemi (Casola Lab, ECCB, Texas A&M University)
Version: 1.0.0
"""

import argparse
import os
import gzip
import sys
import datetime

__version__ = "1.0.0"
__author__ = "Adekola Owoyemi (Casola Lab, ECCB, Texas A&M University)"

VALID_EXT = (".fa", ".fasta", ".faa", ".fna")

# ------------------ Utilities ------------------

def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def clean_header(header):
    return header.split()[0]

def add_species(header, species):
    return f"{header}__{species.replace(' ', '_')}"

def remove_species(header):
    return header.split("__")[0]

def load_map(path):
    m = {}
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    m[parts[0]] = parts[1]
    return m

def log(msg):
    t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{t}] [INFO] {msg}", file=sys.stderr)
# ------------------ Core ------------------

def process_file(infile, outfile, args, species_map=None):
    base = os.path.basename(infile).split(".")[0]

    # Determine species
    species = None
    if species_map:
        species = species_map.get(base, base)
    elif args.add_species:
        species = args.add_species

    log(f"Processing: {infile}")

    with open_file(infile, "rt") as fin, open_file(outfile, "wt") as fout:
        seq_buffer = []

        for line in fin:
            if line.startswith(">"):
                # flush buffer if linearizing
                if args.linearize and seq_buffer:
                    fout.write("".join(seq_buffer) + "\n")
                    seq_buffer = []

                hdr = clean_header(line[1:].strip())

                if args.remove_species:
                    hdr = remove_species(hdr)

                if species:
                    hdr = add_species(hdr, species)

                fout.write(f">{hdr}\n")

            else:
                if args.linearize:
                    seq_buffer.append(line.strip())
                else:
                    fout.write(line)

        # flush last sequence
        if args.linearize and seq_buffer:
            fout.write("".join(seq_buffer) + "\n")
    
    log(f"Output written to: {outfile}")

# ------------------ Main ------------------

def main():
    parser = argparse.ArgumentParser(description="Flexible FASTA header tool")

    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)

    parser.add_argument("--add-species", help="Add species name")
    parser.add_argument("--add-species-map", help="Map file (basename -> species)")
    parser.add_argument("--remove-species", action="store_true")

    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--linearize", action="store_true")

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}"
    )

    args = parser.parse_args()

    # ✅ Print version
    log(f"Running fasta_header_tool.py v{__version__}")

    # ✅ conflict check
    if args.add_species and args.add_species_map:
        print("[ERROR] Choose either --add-species or --add-species-map", file=sys.stderr)
        sys.exit(1)

    species_map = load_map(args.add_species_map) if args.add_species_map else None

    # ✅ directory mode
    if os.path.isdir(args.input):

        if not os.path.exists(args.output):
            os.makedirs(args.output)

        # ✅ choose traversal mode
        if args.recursive:
            dirs = os.walk(args.input)
        else:
            dirs = [(args.input, [], os.listdir(args.input))]

        for root, _, files in dirs:
            for f in files:
                if f.endswith(VALID_EXT):
                    infile = os.path.join(root, f)
                    rel_path = os.path.relpath(infile, args.input)
                    outfile = os.path.join(args.output, rel_path)

                    os.makedirs(os.path.dirname(outfile), exist_ok=True)

                    process_file(infile, outfile, args, species_map)

    # ✅ single file mode
    else:
        process_file(args.input, args.output, args, species_map)

if __name__ == "__main__":
    main()
