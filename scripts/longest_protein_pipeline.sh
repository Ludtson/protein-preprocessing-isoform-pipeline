#!/usr/bin/env bash
set -euo pipefail

# ==========================================
# LINGUA Protein Pipeline
# Version: 1.0.0
# Author: Adekola Owoyemi (Casola Lab, Texas A&M University)
# ==========================================

# --- SCRIPT DESCRIPTION ---
#
# longest_protein_pipeline.sh
#
# This script automates a three-step protein analysis pipeline for multiple species.
#
# The pipeline for each species consists of:
# 1. Protein QC: Cleans input FASTA files (using `qc_proteins.py`)
# 2. Isoform Selection: Selects the longest protein isoform per gene (using `select_longest_proteins.py`)
# 3. Optional FASTA header formatting (using `fasta_header_tool.py`)
#
# The script processes a list of species from a map file, locates the corresponding
# input files, and orchestrates the entire workflow. It generates individual output
# files for each species and compiles two master summary reports for all species.
#
# Main outputs:
# 1. `master_summary_report.tsv`: A high-level summary of processed files and gene/protein counts.
# 2. `master_longest_isoforms.tsv`: A combined detailed report from all species.
# 3. Individual files for each species saved in dedicated subdirectories.
#
# Usage:
#   bash longest_protein_pipeline.sh -m <map_file> -g <gff_dir> -f <fasta_dir> [options]
#
# Example:
#   bash longest_protein_pipeline.sh -m species.map -g gff_files -f fasta_files -o pipeline_output
#
# --- Script Directory (makes the script portable) ---
# This command finds the directory where the current script is located.
# This allows the script to find the Python files even if they are not in the PATH.
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

# --- DEFAULT VALUES ---
OUTPUT_DIR=""
GFF_EXT=".gff"
FASTA_EXT=".faa"
WITH_TRANSCRIPTS=false

SPECIES_MAP_FILE=""
CDS_DIR=""
CDNA_DIR=""


# Now, the Python scripts are referenced by their path relative to the script's location.
QC_PYTHON_CMD="$SCRIPT_DIR/../bin/qc_proteins.py"
ISO_SEL_PYTHON_CMD="$SCRIPT_DIR/../bin/select_longest_proteins.py"
HEADER_TOOL_CMD="$SCRIPT_DIR/../bin/fasta_header_tool.py"

# --- Function for Help Message ---
show_help() {
    echo "Usage: $0 -m <map_file> -g <gff_dir> -f <fasta_dir> [options]"
    echo ""
    echo "Runs the full protein isoform selection pipeline for all species defined in a map file."
    echo ""
    echo "Required arguments:"
    echo "  -m, --map-file <file>     Tab-separated file with columns: species_name and file_basename"
    echo "  -g, --gff-dir <dir>       Directory containing GFF files."
    echo "  -f, --fasta-dir <dir>     Directory containing original protein FASTA files."
    echo ""
    echo "Optional arguments:"
    echo "  -o, --out-dir <dir>       Main output directory for all results. (Default: ./analysis_results)"
    echo "  -k, --keep_noM_only       Keep sequences with only 'noM' defect"    
    echo "  -l, --min-len <int>       Minimum protein length cutoff for QC (default: 30)"
    echo " --format-headers           Enable FASTA header formatting"
    echo " --add-species <name>       Add species name to FASTA headers"
    echo " --remove-species           Remove species suffix from headers"
    echo ""
    echo "Optional transcript extraction:"
    echo "  --with-transcripts        Enable CDS/cDNA extraction"
    echo "  --cds-dir <dir>           Directory containing CDS FASTA files"
    echo "  --cdna-dir <dir>          Directory containing cDNA FASTA files"
    echo ""
    echo " -h, --help                 Show this help message and exit."
}

# --- Parse Command-Line Options ---
MIN_LEN=30  # Default value
FORMAT_HEADERS=false
ADD_SPECIES=""
REMOVE_SPECIES=false
KEEP_NOM_FLAG=""  # Default: don't pass -k
THREADS=1   # default (safe, sequential)

while [[ "$#" -gt 0 ]]; do
  case $1 in
    -m|--map-file) SPECIES_MAP_FILE="$2"; shift ;;
    -g|--gff-dir) GFF_DIR="$2"; shift ;;
    -f|--fasta-dir) FASTA_DIR="$2"; shift ;;
    -o|--out-dir) OUTPUT_DIR="$2"; shift ;;
    -l|--min-len) MIN_LEN="$2"; shift ;;
    -k|--keep_noM_only) KEEP_NOM_FLAG="-k" ;;  # <-- Add this line
    -h|--help) show_help; exit 0 ;;
    --format-headers) FORMAT_HEADERS=true ;;
    --add-species) ADD_SPECIES="$2"; shift ;;
    --remove-species) REMOVE_SPECIES=true ;;
    --with-transcripts) WITH_TRANSCRIPTS=true ;;
    --cds-dir) CDS_DIR="$2"; shift ;;
    --cdna-dir) CDNA_DIR="$2"; shift ;;
    --threads) THREADS="$2"; shift ;;
    --version) echo "LINGUA Protein Pipeline v1.0.0"; exit 0 ;;
    *) echo "Unknown parameter passed: $1"; show_help; exit 1 ;;
  esac
  shift
done

# --- Validate Required Arguments ---
if [ -z "$GFF_DIR" ] || [ -z "$FASTA_DIR" ]; then

    echo "Error: Missing required arguments: --gff-dir and --fasta-dir"
    show_help
    exit 1
fi

if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="analysis_results"
    echo "[INFO] No output directory specified. Using default: $OUTPUT_DIR"
else
    echo "[INFO] Using user-specified output directory: $OUTPUT_DIR"
fi

if [ "$WITH_TRANSCRIPTS" = true ]; then
    if [ -z "$CDS_DIR" ] && [ -z "$CDNA_DIR" ]; then
        echo "[ERROR] --with-transcripts requires at least --cds-dir or --cdna-dir"
        exit 1
    fi
fi

# --- Configuration ---
GFF_DIR=${GFF_DIR%/}
FASTA_DIR=${FASTA_DIR%/}
OUTPUT_DIR=${OUTPUT_DIR%/}
SUMMARY_FILE="${OUTPUT_DIR}/master_summary_report.tsv"
MASTER_DETAILED_TSV="${OUTPUT_DIR}/master_longest_isoforms.tsv"
LOG_FILE="${OUTPUT_DIR}/pipeline.log"
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "$LOG_DIR"

if [ "$WITH_TRANSCRIPTS" = true ]; then
    TRANSCRIPTS_SUBDIR="${OUTPUT_DIR}/optional_longest_transcripts"
    mkdir -p "$TRANSCRIPTS_SUBDIR"
fi

process_one_species () {

    local species_name="$1"
    local base_name="$2"

    input_fasta="${FASTA_DIR}/${base_name}${FASTA_EXT}"
    input_gff="${GFF_DIR}/${base_name}.gff"
    if [ ! -f "$input_gff" ]; then
        input_gff="${GFF_DIR}/${base_name}.gff3"
    fi

    if [ ! -f "$input_fasta" ] || [ ! -f "$input_gff" ]; then
        echo "[WARNING] Missing input files for ${base_name}, skipping" >> "$LOG_FILE"
        return
    fi

    # -------- Step 1: QC --------
    output_cleaned_fasta="${CLEANED_FASTA_SUBDIR}/${base_name}_cleaned${FASTA_EXT}"
    output_qc_tsv="${QC_DETAILS_SUBDIR}/${base_name}_qc_details.tsv"

    qc_summary=$("$QC_PYTHON_CMD" \
        -i "$input_fasta" \
        -o "$output_cleaned_fasta" \
        -d "$output_qc_tsv" \
        -s "$species_name" \
        -l "$MIN_LEN" \
        $KEEP_NOM_FLAG 2>>"$LOG_FILE")

    # ✅ FIX: preserve original delimiter (NO whitespace collapsing)
    qc_line=$(echo "$qc_summary" | grep "^${species_name}[[:space:]]")

    # -------- Step 2: Isoform --------
    output_iso_tsv="${ISOFORM_DETAILS_SUBDIR}/${base_name}_isoform_details.tsv"
    output_primary_fasta="${PRIMARY_FASTA_SUBDIR}/${base_name}_longest_proteins${FASTA_EXT}"

    iso_summary=$("$ISO_SEL_PYTHON_CMD" \
        -g "$input_gff" \
        -f "$output_cleaned_fasta" \
        -v "$output_iso_tsv" \
        -o "$output_primary_fasta" \
        -s "$species_name" 2>>"$LOG_FILE")

    iso_line=$(echo "$iso_summary" | grep "^${species_name}[[:space:]]")

    # ✅ guard (unchanged)
    if [ -z "$qc_line" ] || [ -z "$iso_line" ]; then
        echo "[WARNING] Missing summary for ${base_name}, skipping" >> "$LOG_FILE"
        return
    fi

    # ✅ FIX: enforce clean TSV extraction
    iso_trimmed=$(echo "$iso_line" | awk -v OFS='\t' '{print $2, $3, $4}')

    # ✅ FIX: safe, consistent TSV write
    printf "%s\t%s\n" "$qc_line" "$iso_trimmed" > "$OUTPUT_DIR/tmp_${base_name}_summary.tsv"

    if [ ! -s "$output_primary_fasta" ]; then
        echo "[ERROR] Isoform output empty for ${base_name}, skipping" >> "$LOG_FILE"
        return
    fi

    # -------- Step 3: Header (optional) --------
    output_final_fasta="${FINAL_FASTA_SUBDIR}/${base_name}_final${FASTA_EXT}"

    if [ "$FORMAT_HEADERS" = true ]; then
        "$HEADER_TOOL_CMD" -i "$output_primary_fasta" -o "$output_final_fasta" >> "$LOG_FILE" 2>&1
    else
        cp "$output_primary_fasta" "$output_final_fasta"
    fi

    # -------- Step 4: Transcript extraction -------- ✅ KEEP THIS
    if [ "$WITH_TRANSCRIPTS" = true ]; then

        cds_input=()
        cdna_input=()
        cds_out=()
        cdna_out=()

        if [ -n "$CDS_DIR" ]; then
            for ext in fna fa fasta; do
                if [ -f "${CDS_DIR}/${base_name}.${ext}" ]; then
                    cds_input=(--cds "${CDS_DIR}/${base_name}.${ext}")
                    cds_out=(--out_cds "${TRANSCRIPTS_SUBDIR}/${base_name}.fna")
                    break
                fi
            done
        fi

        if [ -n "$CDNA_DIR" ]; then
            for ext in fna fa fasta; do
                if [ -f "${CDNA_DIR}/${base_name}.${ext}" ]; then
                    cdna_input=(--cdna "${CDNA_DIR}/${base_name}.${ext}")
                    cdna_out=(--out_cdna "${TRANSCRIPTS_SUBDIR}/${base_name}.cdna.fasta")
                    break
                fi
            done
        fi

        if [ ${#cds_input[@]} -gt 0 ] || [ ${#cdna_input[@]} -gt 0 ]; then
            "$SCRIPT_DIR/../bin/extract_longest_transcripts.py" \
                --gff "$input_gff" \
                --protein "$output_primary_fasta" \
                "${cds_input[@]}" \
                "${cdna_input[@]}" \
                "${cds_out[@]}" \
                "${cdna_out[@]}" \
                >> "$LOG_FILE" 2>&1
        else
            echo "[WARNING] No CDS/cDNA file found for ${base_name}, skipping transcript extraction" >> "$LOG_FILE"
        fi
    fi
}

# =========================================
# OUTPUT STRUCTURE SETUP
# =========================================

CLEANED_FASTA_SUBDIR="${OUTPUT_DIR}/1_cleaned_proteins"
QC_DETAILS_SUBDIR="${OUTPUT_DIR}/2_qc_details"
ISOFORM_DETAILS_SUBDIR="${OUTPUT_DIR}/3_isoform_details"
PRIMARY_FASTA_SUBDIR="${OUTPUT_DIR}/4_longest_proteins"
FINAL_FASTA_SUBDIR="${OUTPUT_DIR}/5_final_proteins"

# =========================================
# SCRIPT START
# =========================================

if [ -n "$SPECIES_MAP_FILE" ]; then
    echo "[INFO] Starting analysis using species map: $SPECIES_MAP_FILE"
else
    echo "[INFO] Starting analysis without species map"
fi

# Create output directories
mkdir -p "$CLEANED_FASTA_SUBDIR"
mkdir -p "$QC_DETAILS_SUBDIR"
mkdir -p "$ISOFORM_DETAILS_SUBDIR"
mkdir -p "$PRIMARY_FASTA_SUBDIR"
mkdir -p "$FINAL_FASTA_SUBDIR"

echo "$(date "+%Y-%m-%d %H-%M-%S") [INFO] Output will be saved in the '${OUTPUT_DIR}' directory."

# =========================================
# MAIN EXECUTION
# =========================================

if [ -n "$SPECIES_MAP_FILE" ]; then

    echo "$(date "+%Y-%m-%d %H-%M-%S") [INFO] Running with species map"

    DELIM=$(head -n 1 "$SPECIES_MAP_FILE" | grep -q $'\t' && echo $'\t' || echo ',')

    tail -n +2 "$SPECIES_MAP_FILE" | tr -d '\r' | while IFS="$DELIM" read -r species_name base_name || [[ -n "$species_name" ]]; do

        species_name=$(echo "$species_name" | xargs)
        base_name=$(echo "$base_name" | xargs)

        process_one_species "$species_name" "$base_name"

    done

else

    echo "$(date "+%Y-%m-%d %H-%M-%S") [INFO] Running without species map..."
    echo "[INFO] Running with ${THREADS} parallel jobs"

    set +e  # allow failures

    job_count=0

    for fasta_file in "$FASTA_DIR"/*"$FASTA_EXT"; do
        base_name=$(basename "$fasta_file" "$FASTA_EXT")
        species_name="$base_name"

        process_one_species "$species_name" "$base_name" \
            > "$LOG_DIR/${base_name}.log" 2>&1 &

        ((job_count++))

        # simple throttle (robust)
        if (( job_count % THREADS == 0 )); then
            wait
        fi
    done

    wait  # final wait

    set -e

fi

# =========================================
# BUILD MASTER REPORTS
# =========================================

echo "[INFO] Building master summary..."

echo -e "Species\tN_Input_Proteins\tN_Cleaned_Proteins\tFiltered_Total\tS\tnoM\t<30\tS/noM\tS/<30\tS/noM/<30\tnoM/<30\tN_input_proteins\tN_genes_processed\tN_output_proteins" > "$SUMMARY_FILE"

for file in "$OUTPUT_DIR"/tmp_*_summary.tsv; do
    [ -s "$file" ] || continue
    cat "$file" >> "$SUMMARY_FILE"
done

# optional cleanup
rm -f "$OUTPUT_DIR"/tmp_*_summary.tsv

# -------- MASTER DETAILED --------
echo -e "Species\tGene_ID\tTranscript_ID\tProtein_ID\tProtein_Length\tmRNA_length\tCDS_Length\tDelta_nt\tStopFeature\tFlags\tGene_Strand\tGene_Isoforms\tmRNA_exons" > "$MASTER_DETAILED_TSV"

for file in "$ISOFORM_DETAILS_SUBDIR"/*_isoform_details.tsv; do
    species=$(basename "$file" _isoform_details.tsv)
    awk -v sp="$species" 'NR>1 {print sp "\t" $0}' "$file" >> "$MASTER_DETAILED_TSV"
done

echo "----------------------------------------------------"
echo "$(date "+%Y-%m-%d %H-%M-%S") [INFO] ✅ All species processed. Master summary report: ${SUMMARY_FILE}"
echo "$(date "+%Y-%m-%d %H-%M-%S") [INFO] ✅ Master detailed report: ${MASTER_DETAILED_TSV}"
