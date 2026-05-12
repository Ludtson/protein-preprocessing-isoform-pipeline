#!/usr/bin/env bash

# This script is a wrapper to run qc_protein.v.0.2.py on all FASTA files in a directory.
# It handles mapping filenames to species names and organizes all outputs.

# --- Safety checks ---
set -e
set -u
set -o pipefail

# --- Help Message ---
usage() {
    echo "Usage: $0 -i <input_dir> -m <map_file> -o <output_dir> [-l <min_length>]"
    echo "  -i <input_dir>   : Directory containing input FASTA files (*.faa)."
    echo "  -m <map_file>    : Tab-delimited file mapping FASTA basenames to full species names."
    echo "  -o <output_dir>  : Directory where all outputs will be saved."
    echo "  -l <min_length>  : (Optional) Minimum protein length to pass QC (default: 30)."
    exit 1
}

# --- Default values ---
MIN_LENGTH=30

# --- Parse Command-Line Arguments ---
while getopts ":i:m:o:l:" opt; do
    case ${opt} in
        i) INPUT_DIR=$OPTARG ;;
        m) MAP_FILE=$OPTARG ;;
        o) OUTPUT_DIR=$OPTARG ;;
        l) MIN_LENGTH=$OPTARG ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Check if mandatory arguments were provided
if [[ -z "${INPUT_DIR-}" || -z "${MAP_FILE-}" || -z "${OUTPUT_DIR-}" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# --- Get the absolute path of the directory containing this script ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Setup ---
echo "--- Setting up output directories ---"
PROTEIN_OUT_DIR="$OUTPUT_DIR/protein"
ISSUES_OUT_DIR="$OUTPUT_DIR/issues"
mkdir -p "$PROTEIN_OUT_DIR"
mkdir -p "$ISSUES_OUT_DIR"
echo "  Cleaned proteins will be saved in: $PROTEIN_OUT_DIR"
echo "  Detailed issue logs will be saved in: $ISSUES_OUT_DIR"

# --- Load Species Map into an associative array for fast lookups ---
echo "--- Loading species map from: $MAP_FILE ---"
declare -A species_map
while IFS=$'\t' read -r fname species; do
    [[ "$fname" == "FName" ]] && continue
    if [[ -n "$fname" && -n "$species" ]]; then
        species_map["$fname"]="$species"
    fi
done < "$MAP_FILE"
echo "  Loaded ${#species_map[@]} mappings."

# --- Main Processing Loop ---
echo -e "\n--- Starting protein QC for all FASTA files in: $INPUT_DIR ---"

for fasta_file in "$INPUT_DIR"/*.faa; do
    [ -e "$fasta_file" ] || { echo "Error: No .faa files found in '$INPUT_DIR'."; exit 1; }
    basename=$(basename "$fasta_file" .faa)
    species_name=${species_map[$basename]:-"$basename"}
    
    echo -e "\n  Processing file: $fasta_file"
    echo "    -> Species: $species_name"

    safe_species_name=$(echo "$species_name" | tr ' ' '_')
    output_faa="$PROTEIN_OUT_DIR/${safe_species_name}.faa"
    detailed_tsv="$ISSUES_OUT_DIR/${safe_species_name}.tsv"

    # Run the qc_protein script. Its summary will print to the screen,
    # and its detailed TSV will be saved.
    "$SCRIPT_DIR/qc_protein.v.0.2.py" \
        -i "$fasta_file" \
        -o "$output_faa" \
        -s "$species_name" \
        -d "$detailed_tsv" \
        -l "$MIN_LENGTH"
done

# --- NEW: Consolidate all detailed issue logs into a single file ---
echo -e "\n--- Consolidating all detailed issue logs ---"
FINAL_ISSUES_TABLE="$OUTPUT_DIR/all_species_qc_issues.tsv"
is_first_file=true

# Find all tsv files in the issues directory
shopt -s nullglob
issue_files=("$ISSUES_OUT_DIR"/*.tsv)
shopt -u nullglob

if [ ${#issue_files[@]} -gt 0 ]; then
    for tsv_file in "${issue_files[@]}"; do
        if [ "$is_first_file" = true ]; then
            # For the first file, copy the whole thing (including the header)
            cat "$tsv_file" > "$FINAL_ISSUES_TABLE"
            is_first_file=false
        else
            # For subsequent files, append everything *except* the header
            tail -n +2 "$tsv_file" >> "$FINAL_ISSUES_TABLE"
        fi
    done
    echo "  Successfully combined ${#issue_files[@]} issue logs."
else
    echo "  Warning: No issue logs were found to consolidate."
fi

echo -e "\n--- Pipeline Finished ---"
echo "All FASTA files have been processed."
echo "Final consolidated issues table saved to: $FINAL_ISSUES_TABLE"

