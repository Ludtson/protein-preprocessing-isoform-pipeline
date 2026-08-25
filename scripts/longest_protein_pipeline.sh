#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ==========================================
# LINGUA Protein Pipeline
# Version: 1.0.0
# Author: Adekola Owoyemi (Casola Lab, Texas A&M University)
# ==========================================

# --- SCRIPT DESCRIPTION ---
#
# longest_protein_pipeline.sh
#
# Runs the protein QC and isoform selection pipeline for one or more species.
#
# The pipeline for each species consists of:
# 1. Protein QC: Cleans input FASTA files (using `qc_proteins.py`)
# 2. Isoform Selection: Selects the longest protein isoform per gene (using `select_longest_proteins.py`)
# 3. Optional FASTA header formatting (using `fasta_header_tool.py`)
# 4. Optional CDS/cDNA extraction for the selected isoforms (using `extract_longest_transcripts.py`)
#
# Species can be supplied via a tab- or comma-separated map file (species name,
# file basename), or the script can run directly over every FASTA file in
# --fasta-dir, using each filename (without extension) as the species name.
#
# Main outputs:
# 1. `master_summary_report.tsv`: A high-level summary of processed files and gene/protein counts.
# 2. `master_longest_isoforms.tsv`: A combined detailed report from all species.
# 3. `run_summary.txt`: Aggregate run-level metrics.
# 4. Individual files for each species saved in dedicated subdirectories.
#
# Usage:
#   bash longest_protein_pipeline.sh -g <gff_dir> -f <fasta_dir> [-m <map_file>] [options]
#
# Example:
#   bash longest_protein_pipeline.sh -m species.map -g gff_files -f fasta_files -o pipeline_output
#
# --- Script Directory (makes the script portable) ---
# This command finds the directory where the current script is located.
# This allows the script to find the Python files even if they are not in the PATH.
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

SECONDS=0

# --- DEFAULT VALUES ---
OUTPUT_DIR=""
GFF_DIR=""
FASTA_DIR=""
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

# --- Logging helper ---
# Prints a timestamped, leveled message to stdout, and also appends it to
# $LOG_FILE once that path has been established.
log() {
    local level="$1"; shift
    local ts
    ts=$(date "+%Y-%m-%d %H:%M:%S")
    echo "${ts} [${level}] $*"
    if [ -n "${LOG_FILE:-}" ]; then
        echo "${ts} [${level}] $*" >> "$LOG_FILE"
    fi
}

# --- Function for Help Message ---
show_help() {
    echo "Usage: $0 -g <gff_dir> -f <fasta_dir> [-m <map_file>] [options]"
    echo ""
    echo "Runs the protein isoform selection pipeline for one or more species."
    echo ""
    echo "Required arguments:"
    echo "  -g, --gff-dir <dir>       Directory containing GFF files."
    echo "  -f, --fasta-dir <dir>     Directory containing original protein FASTA files."
    echo ""
    echo "Optional arguments:"
    echo "  -m, --map-file <file>     Tab- or comma-separated file with columns: species_name, file_basename."
    echo "                            If omitted, every FASTA file in --fasta-dir is processed, using its"
    echo "                            filename (without extension) as the species name."
    echo "  -o, --out-dir <dir>       Main output directory for all results. (Default: ./analysis_results)"
    echo "  -k, --keep_noM_only       Keep sequences with only 'noM' defect"
    echo "  -l, --min-len <int>       Minimum protein length cutoff for QC (default: 30)"
    echo "  --format-headers          Enable FASTA header formatting"
    echo "  --add-species             Tag FASTA headers with the species name (the map-file species name,"
    echo "                            or the filename without extension when no map file is used)"
    echo "  --remove-species          Remove species suffix from headers"
    echo "  --threads <int>           Number of species to process in parallel (default: 1)"
    echo ""
    echo "Optional transcript extraction:"
    echo "  --with-transcripts        Enable CDS/cDNA extraction"
    echo "  --cds-dir <dir>           Directory containing CDS FASTA files"
    echo "  --cdna-dir <dir>          Directory containing cDNA FASTA files"
    echo ""
    echo "  -h, --help                Show this help message and exit."
}

# --- Parse Command-Line Options ---
MIN_LEN=30  # Default value
FORMAT_HEADERS=false
ADD_SPECIES=false
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
    -k|--keep_noM_only) KEEP_NOM_FLAG="-k" ;;
    -h|--help) show_help; exit 0 ;;
    --format-headers) FORMAT_HEADERS=true ;;
    --add-species) ADD_SPECIES=true ;;
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
    missing=()
    [ -z "$GFF_DIR" ] && missing+=("--gff-dir")
    [ -z "$FASTA_DIR" ] && missing+=("--fasta-dir")
    echo "Error: Missing required argument(s): ${missing[*]}"
    show_help
    exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [ "$THREADS" -lt 1 ]; then
    echo "Error: --threads must be a positive integer (got: '$THREADS')"
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
RUN_SUMMARY_FILE="${OUTPUT_DIR}/run_summary.txt"
LOG_FILE="${OUTPUT_DIR}/pipeline.log"
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "$LOG_DIR"

if [ "$WITH_TRANSCRIPTS" = true ]; then
    CDS_SUBDIR="${OUTPUT_DIR}/6_longest_cds"
    CDNA_SUBDIR="${OUTPUT_DIR}/7_longest_cdna"
    [ -n "$CDS_DIR" ] && mkdir -p "$CDS_SUBDIR"
    [ -n "$CDNA_DIR" ] && mkdir -p "$CDNA_SUBDIR"
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
        log WARNING "Missing input files for ${base_name}, skipping"
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

    # Preserve the original tab delimiter (no whitespace collapsing).
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

    if [ -z "$qc_line" ] || [ -z "$iso_line" ]; then
        log WARNING "Missing summary for ${base_name}, skipping"
        return
    fi

    # Field-split strictly on tabs: a species name containing a space (e.g.
    # "Arabidopsis thaliana") would otherwise shift every field under awk's
    # default whitespace splitting.
    iso_trimmed=$(echo "$iso_line" | awk -F'\t' -v OFS='\t' '{print $2, $3, $4}')

    printf "%s\t%s\n" "$qc_line" "$iso_trimmed" > "$OUTPUT_DIR/tmp_${base_name}_summary.tsv"

    if [ ! -s "$output_primary_fasta" ]; then
        log ERROR "Isoform output empty for ${base_name}, skipping"
        return
    fi

    # -------- Step 3: Header (optional) --------
    output_final_fasta="${FINAL_FASTA_SUBDIR}/${base_name}_final${FASTA_EXT}"

    if [ "$FORMAT_HEADERS" = true ]; then
        header_args=()
        if [ "$ADD_SPECIES" = true ]; then
            header_args+=(--add-species "$species_name")
        fi
        if [ "$REMOVE_SPECIES" = true ]; then
            header_args+=(--remove-species)
        fi
        "$HEADER_TOOL_CMD" -i "$output_primary_fasta" -o "$output_final_fasta" "${header_args[@]}" >> "$LOG_FILE" 2>&1
    else
        cp "$output_primary_fasta" "$output_final_fasta"
    fi

    # -------- Step 4: Transcript extraction --------
    if [ "$WITH_TRANSCRIPTS" = true ]; then

        cds_input=()
        cdna_input=()
        cds_out=()
        cdna_out=()

        if [ -n "$CDS_DIR" ]; then
            for ext in fna fa fasta; do
                if [ -f "${CDS_DIR}/${base_name}.${ext}" ]; then
                    cds_input=(--cds "${CDS_DIR}/${base_name}.${ext}")
                    cds_out=(--out_cds "${CDS_SUBDIR}/${base_name}.fna")
                    break
                fi
            done
        fi

        if [ -n "$CDNA_DIR" ]; then
            for ext in fna fa fasta; do
                if [ -f "${CDNA_DIR}/${base_name}.${ext}" ]; then
                    cdna_input=(--cdna "${CDNA_DIR}/${base_name}.${ext}")
                    cdna_out=(--out_cdna "${CDNA_SUBDIR}/${base_name}.cdna.fasta")
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
            log WARNING "No CDS/cDNA file found for ${base_name}, skipping transcript extraction"
        fi
    fi
}

# A species job runs in its own subshell so failures are isolated: stdout
# and stderr are captured to a per-species log, and a single-line status is
# echoed to the main script's stdout (unredirected) once it finishes, so a
# multi-species run shows live progress instead of going silent.
launch_job () {
    local species_name="$1"
    local base_name="$2"
    (
        if process_one_species "$species_name" "$base_name" > "$LOG_DIR/${base_name}.log" 2>&1; then
            status=0
        else
            status=$?
        fi

        ts=$(date "+%Y-%m-%d %H:%M:%S")
        if [ "$status" -ne 0 ]; then
            echo "${ts} [ERROR] ${species_name}: processing failed (exit ${status}) - see ${LOG_DIR}/${base_name}.log"
        elif grep -qE '\[(WARNING|ERROR)\]' "$LOG_DIR/${base_name}.log" 2>/dev/null; then
            echo "${ts} [WARNING] ${species_name}: completed with warnings - see ${LOG_DIR}/${base_name}.log"
        else
            echo "${ts} [INFO] Done: ${species_name}"
        fi
    ) &
}

# Runs $1 (a shell snippet that calls launch_job per species) while keeping
# at most $THREADS jobs in flight, using a rolling pool (wait -n) rather
# than fixed-size sequential batches.
run_pool () {
    local running=0
    while IFS= read -r -u 3 line; do
        [ -z "$line" ] && continue
        local species_name base_name
        species_name="${line%%$'\t'*}"
        base_name="${line#*$'\t'}"

        launch_job "$species_name" "$base_name"

        running=$((running + 1))
        if [ "$running" -ge "$THREADS" ]; then
            wait -n
            running=$((running - 1))
        fi
    done 3<"$1"
    wait
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
    log INFO "Starting analysis using species map: $SPECIES_MAP_FILE"
else
    log INFO "Starting analysis without species map"
fi

# Create output directories
mkdir -p "$CLEANED_FASTA_SUBDIR"
mkdir -p "$QC_DETAILS_SUBDIR"
mkdir -p "$ISOFORM_DETAILS_SUBDIR"
mkdir -p "$PRIMARY_FASTA_SUBDIR"
mkdir -p "$FINAL_FASTA_SUBDIR"

log INFO "Output will be saved in the '${OUTPUT_DIR}' directory."

# =========================================
# MAIN EXECUTION
# =========================================

JOB_LIST=$(mktemp)
trap 'rm -f "$JOB_LIST"' EXIT

if [ -n "$SPECIES_MAP_FILE" ]; then

    log INFO "Running with species map (threads: ${THREADS})"

    DELIM=$(head -n 1 "$SPECIES_MAP_FILE" | grep -q $'\t' && echo $'\t' || echo ',')

    tail -n +2 "$SPECIES_MAP_FILE" | tr -d '\r' | while IFS="$DELIM" read -r species_name base_name || [[ -n "$species_name" ]]; do
        species_name=$(echo "$species_name" | xargs)
        base_name=$(echo "$base_name" | xargs)
        [ -z "$species_name" ] && continue
        printf '%s\t%s\n' "$species_name" "$base_name"
    done > "$JOB_LIST"

else

    log INFO "Running without species map (threads: ${THREADS})"

    for fasta_file in "$FASTA_DIR"/*"$FASTA_EXT"; do
        [ -e "$fasta_file" ] || continue
        base_name=$(basename "$fasta_file" "$FASTA_EXT")
        printf '%s\t%s\n' "$base_name" "$base_name"
    done > "$JOB_LIST"

fi

run_pool "$JOB_LIST"

# =========================================
# BUILD MASTER REPORTS
# =========================================

log INFO "Building master summary..."

# N_Cleaned_Proteins (from qc_proteins.py) and N_Isoform_Step_Input (from
# select_longest_proteins.py) are computed independently by two separate
# scripts counting the same cleaned FASTA file, and should therefore always
# be equal. This is a deliberate cross-check, not a redundant column: if a
# row ever shows them differing, that means the isoform-selection step read
# a different cleaned FASTA than the QC step actually wrote (stale file,
# truncated write, wrong path), and is worth investigating.
echo -e "Species\tN_Input_Proteins\tN_Cleaned_Proteins\tFiltered_Total\tS\tnoM\t<30\tS/noM\tS/<30\tS/noM/<30\tnoM/<30\tN_Isoform_Step_Input\tN_genes_processed\tN_output_proteins" > "$SUMMARY_FILE"

for file in "$OUTPUT_DIR"/tmp_*_summary.tsv; do
    [ -s "$file" ] || continue
    cat "$file" >> "$SUMMARY_FILE"
done

# optional cleanup
rm -f "$OUTPUT_DIR"/tmp_*_summary.tsv

n_processed=$(( $(wc -l < "$SUMMARY_FILE") - 1 ))
if [ "$n_processed" -le 0 ]; then
    log ERROR "No species were successfully processed. Check ${LOG_FILE} and per-species logs under ${LOG_DIR}."
    exit 1
fi

# -------- MASTER DETAILED --------
echo -e "Species\tGene_ID\tTranscript_ID\tProtein_ID\tProtein_Length\tmRNA_length\tCDS_Length\tDelta_nt\tStopFeature\tFlags\tGene_Strand\tGene_Isoforms\tmRNA_exons" > "$MASTER_DETAILED_TSV"

for file in "$ISOFORM_DETAILS_SUBDIR"/*_isoform_details.tsv; do
    species=$(basename "$file" _isoform_details.tsv)
    awk -v sp="$species" 'NR>1 {print sp "\t" $0}' "$file" >> "$MASTER_DETAILED_TSV"
done

# =========================================
# RUN SUMMARY
# =========================================

# Only count species skipped entirely (missing input / missing summary), not
# per-step skips like "skipping transcript extraction" for one optional step.
n_skipped=$(grep -cE "Missing (input files|summary) for .*, skipping" "$LOG_FILE" 2>/dev/null || true)
n_skipped=${n_skipped:-0}

read -r total_input total_cleaned total_output < <(
    awk -F'\t' 'NR>1 {tin+=$2; tcln+=$3; tout+=$14} END {printf "%d %d %d\n", tin+0, tcln+0, tout+0}' "$SUMMARY_FILE"
)

qc_retention="NA"
if [ "$total_input" -gt 0 ]; then
    qc_retention=$(awk -v c="$total_cleaned" -v i="$total_input" 'BEGIN{printf "%.2f", (c/i)*100}')
fi

iso_reduction="NA"
if [ "$total_cleaned" -gt 0 ]; then
    iso_reduction=$(awk -v c="$total_cleaned" -v o="$total_output" 'BEGIN{printf "%.2f", ((c-o)/c)*100}')
fi

elapsed=$SECONDS
elapsed_fmt=$(printf '%02dh:%02dm:%02ds' $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60)))

{
    echo "Run summary"
    echo "  Species processed:      ${n_processed}"
    echo "  Species skipped:        ${n_skipped}"
    echo "  Total input proteins:   ${total_input}"
    echo "  Total cleaned proteins: ${total_cleaned}"
    echo "  Total final proteins:   ${total_output}"
    echo "  QC retention rate:      ${qc_retention}%"
    echo "  Isoform reduction rate: ${iso_reduction}%"
    echo "  Runtime:                ${elapsed_fmt}"
} | tee "$RUN_SUMMARY_FILE"

log INFO "All species processed. Master summary report: ${SUMMARY_FILE}"
log INFO "Master detailed report: ${MASTER_DETAILED_TSV}"
log INFO "Run summary: ${RUN_SUMMARY_FILE}"
