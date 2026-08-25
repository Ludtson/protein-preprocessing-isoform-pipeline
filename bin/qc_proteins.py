#!/usr/bin/env python3
"""
qc_proteins.py

Performs quality control on protein FASTA files.

QC checks and labels:
  - 'S': Internal stop codons ('*') or ambiguous 'X' before the last residue
  - 'noM': Sequence does not start with Methionine ('M')
  - '<30': Sequence length is below minimum cutoff (default: 30 aa)

Outputs:
  1. Filtered FASTA (-o): Sequences that pass QC (optionally retain 'noM'-only)
  2. Detailed TSV (-d): Problematic sequences with QC flags (optional)
  3. Summary TSV (stdout): One-line QC summary per species

Usage:
  python qc_proteins.py -i <input.faa> -o <output.faa> -s <species_name> [-l 30] [-d detailed.tsv] [-k]

Author: Adekola Owoyemi (Casola Lab, ECCB, Texas A&M University)
Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "Adekola Owoyemi (Casola Lab, ECCB, Texas A&M University)"

import sys
import argparse
import gzip
from pathlib import Path
from collections import defaultdict
from typing import Tuple, List, Dict
import datetime

# Define the possible filtering categories
CATS = ["S", "noM", "<30", "S/noM", "S/<30", "S/noM/<30", "noM/<30"]

def timestamp_message(message, stream=sys.stderr):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    stream.write(f"[{current_time}] {message}\n")
    stream.flush()

def qc_sequence(seq: str, min_len: int) -> Tuple[bool, str]:
    """
    Performs QC checks on a single protein sequence.

    Args:
        seq: The protein sequence string.
        min_len: The minimum allowed length.

    Returns:
        A tuple: (is_clean: bool, status_label: str)
    """
    is_stop = False
    is_no_M = False
    is_short = False
    
    # 1. Internal Stop Check
    # Check for '*' or 'X' anywhere except the last residue
    if seq[:-1] and ('*' in seq[:-1] or 'X' in seq[:-1]):
        is_stop = True
    
    # 2. Initial Methionine Check
    if not seq.startswith('M'):
        is_no_M = True
        
    # 3. Short Sequence Check
    if len(seq) < min_len:
        is_short = True

    # Determine status label
    labels = []
    if is_stop:
        labels.append('S')
    if is_no_M:
        labels.append('noM')
    if is_short:
        labels.append('<30')

    status_label = "/".join(labels)
    
    # Determine if the sequence is clean (should be kept by default logic)
    # A sequence is 'clean' if it has no defects.
    # The filtering logic related to 'noM' is handled in the main function
    is_clean = not is_stop and not is_no_M and not is_short
    
    # If there are no defects, status is empty string, but we return True for is_clean
    return is_clean, status_label

def main():
    parser = argparse.ArgumentParser(
        description="Performs quality control on a single protein FASTA file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-i', '--input_fasta', required=True, help='Path to the input protein FASTA file (can be gzipped).')
    parser.add_argument('-o', '--output_fasta', required=True, help='Path to the output filtered protein FASTA file.')
    parser.add_argument('-s', '--species', required=True, help='Species name for the summary report.')
    parser.add_argument('-l', '--min_len', type=int, default=30, help='Minimum protein length cutoff (default: 30 aa).')
    parser.add_argument('-d', '--detailed_tsv', help='Optional path to write a detailed TSV of problematic sequences.')
    parser.add_argument(
        '-k', 
        '--keep_noM_only', 
        action='store_true', 
        default=False, 
        help=(
            'If set, sequences with ONLY the "noM" defect will be KEPT in the output '
            'FASTA file. They will still be counted in the "noM" column of the summary.'
        )
    )

    args = parser.parse_args()

    # --- Initialization ---
    file_summary: Dict[str, int] = {"N_before": 0, "N_after": 0}
    combo_counts: Dict[str, int] = defaultdict(int)
    problem_rows: List[Tuple[str, str, str, str]] = [] # For detailed TSV: (report_name, prot_id, status, seq)
    report_name = args.species # Use species name for the report column

    # --- Read Input and Process ---
    timestamp_message(f"Starting QC for species: {args.species}")
    timestamp_message(f"Minimum length cutoff: {args.min_len}")
    if args.keep_noM_only:
        timestamp_message("NOTE: Sequences with ONLY 'noM' defect will be KEPT in the output FASTA.")

    n_kept = 0
    n_filtered = 0

    try:
        # Determine opener: gzip.open for .gz, built-in open otherwise
        opener = gzip.open if args.input_fasta.endswith('.gz') else open
        
        with opener(args.input_fasta, 'rt', encoding='utf-8') as ifh, \
             open(args.output_fasta, 'w', encoding='utf-8') as ofh:

            header = ""
            seq = ""
            for line in ifh:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    # Process previous entry if it exists
                    if header and seq:
                        file_summary["N_before"] += 1
                        
                        # Run QC checks
                        is_clean, status_label = qc_sequence(seq, args.min_len)
                        
                        # --- Filtering Logic ---
                        is_defective = status_label != ""
                        should_keep = is_clean

                        # Apply the optional rule: keep sequences with ONLY 'noM' defect
                        if is_defective:
                            if status_label == 'noM' and args.keep_noM_only:
                                # This sequence has ONLY 'noM' and the user asked to keep it
                                should_keep = True
                            else:
                                # The sequence has 'S', '<30', or a combination (including noM), 
                                # or it has 'noM' but the user did not enable the keep option.
                                should_keep = False
                        
                        # If sequence should be kept (either passed all QC, or passed the 'noM' exception)
                        if should_keep:
                            ofh.write(f"{header}\n{seq}\n")
                            n_kept += 1
                        else:
                            # If filtered, record the defect status and count
                            n_filtered += 1
                            combo_counts[status_label] += 1
                            if args.detailed_tsv:
                                problem_rows.append((report_name, header.lstrip('>').split()[0], status_label, seq))

                    # Start new entry
                    header = line
                    seq = ""
                else:
                    # Accumulate sequence
                    seq += line

            # Process the last entry in the file
            if header and seq:
                file_summary["N_before"] += 1
                is_clean, status_label = qc_sequence(seq, args.min_len)
                
                # --- Filtering Logic (repeated for last sequence) ---
                is_defective = status_label != ""
                should_keep = is_clean

                if is_defective:
                    if status_label == 'noM' and args.keep_noM_only:
                        should_keep = True
                    else:
                        should_keep = False

                if should_keep:
                    ofh.write(f"{header}\n{seq}\n")
                    n_kept += 1
                else:
                    n_filtered += 1
                    combo_counts[status_label] += 1
                    if args.detailed_tsv:
                        problem_rows.append((report_name, header.lstrip('>').split()[0], status_label, seq))
            
    except Exception as e:
        timestamp_message(f"An error occurred during file processing: {e}")
        sys.exit(1)
        
    file_summary["N_after"] = n_kept

    # --- Print Summary TSV to stdout ---
    timestamp_message("Writing summary TSV to stdout.")

    # Print summary TSV to stdout for the orchestrator
    print("Species\tN_Input_Proteins\tN_Cleaned_Proteins\tFiltered_Total\tS\tnoM\t<30\tS/noM\tS/<30\tS/noM/<30\tnoM/<30")

    # Construct the row using the exact order defined in CATS
    row = [
        args.species,
        str(file_summary["N_before"]),
        str(file_summary["N_after"]),
        str(n_filtered) # Filtered_total is n_filtered
    ]
    
    # Ensure the combined status counts are extracted in the exact order CATS specifies
    for combo_label in CATS:
        row.append(str(combo_counts.get(combo_label, 0)))

    # Print the final summary line
    print("\t".join(row))

    # --- Write detailed TSV if requested ---
    if args.detailed_tsv:
        timestamp_message(f"Writing detailed TSV to {args.detailed_tsv}")
        try:
            with open(Path(args.detailed_tsv), "w", encoding="utf-8") as oh:
                headers = [("Species"), "Protein_ID", "Status", "Sequence"]
                oh.write("\t".join(headers) + "\n")
                for _, prot_id, status, seq in problem_rows:
                    oh.write(f"{report_name}\t{prot_id}\t{status}\t{seq}\n")
        except Exception as e:
            timestamp_message(f"Error writing detailed TSV: {e}")
            # Do not exit, as the main task (summary and filtered FASTA) is complete.

    timestamp_message("QC script finished successfully.")

if __name__ == "__main__":
    main()