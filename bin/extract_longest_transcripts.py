#!/usr/bin/env python3

import re
import sys
import os
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict

"""
extract_longest_transcripts.py

Extracts nucleotide sequences (CDS or cDNA) corresponding to already-selected
longest protein isoforms.

NOTE:
This script does NOT perform isoform selection.
It uses the protein set as the source of truth and retrieves matching
CDS/cDNA sequences.

Input:
  - Protein FASTA (cleaned / longest proteins)
  - GFF annotation
  - CDS and/or cDNA FASTA

Output:
  - FASTA of matched CDS and/or cDNA sequences

Optional:
  - Deduplication logs
  - Missing ID reports
  - Master table (validation/debugging)

Usage (single species):
  python extract_longest_transcripts.py \
    --gff <input.gff> \
    --protein <input.faa> \
    --cds <cds.fna> \
    --cdna <cdna.fna>

Usage (batch mode):
  python extract_longest_transcripts.py \
    --batch \
    --prt-dir <protein_dir> \
    --cds-dir <cds_dir> \
    --cdna-dir <cdna_dir> \
    --gff-dir <gff_dir> \
    --out-dir <output_dir>

Author: Adekola Owoyemi (Casola Lab, ECCB, Texas A&M University)
Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "Adekola Owoyemi (Casola Lab, ECCB, Texas A&M University)"

# --- Core Processing Functions ---

def get_protein_info(protein_fasta_path, load_lengths=False):
    """
    Reads the cleaned Protein FASTA.
    If load_lengths=False: Returns {ProteinID: FullHeader} (Fast mode)
    If load_lengths=True: Returns {ProteinID: {'header': FullHeader, 'len': int}} (Validation mode)
    """
    data = {}
    try:
        with open(protein_fasta_path, 'r', encoding='utf-8', errors='ignore') as f:
            current_id = None
            current_header = None
            current_seq_len = 0
            
            for line in f:
                line = line.strip()
                if not line: continue
                
                if line.startswith('>'):
                    # Save previous if in length mode
                    if load_lengths and current_id:
                        data[current_id] = {'header': current_header, 'len': current_seq_len}
                    
                    current_header = line
                    # Take the first word after '>' as the ID
                    current_id = line.split()[0][1:]
                    current_seq_len = 0
                    
                    if not load_lengths:
                        data[current_id] = current_header
                else:
                    if load_lengths:
                        current_seq_len += len(line)
            
            # Save last entry
            if load_lengths and current_id:
                data[current_id] = {'header': current_header, 'len': current_seq_len}
                
    except FileNotFoundError:
        print(f"Error: Protein file '{protein_fasta_path}' not found.")
        return {}
        
    return data

def parse_gff_map(gff_path):
    """
    Parses GFF to create a comprehensive mapping.
    1. Standard ID mapping (Locus, GeneID -> ProteinID)
    2. Structural Mapping (Exon Coordinates -> ProteinID) for disambiguating isoforms.
    """
    gene_id_to_keys = defaultdict(list)
    mrna_id_to_gene_id = {}
    mrna_id_to_keys = defaultdict(list)
    cds_parent_to_prot = {}
    
    # New: Map mRNA ID to list of Exon Intervals [(start, end), ...]
    mrna_id_to_exons = defaultdict(list)
    
    # Global ambiguity tracker: Key -> Set of ProteinIDs
    key_to_prots = defaultdict(set)
    
    # Location Fingerprint Map: "start-end|start-end..." -> ProteinID
    loc_fingerprint_to_prot = {}

    re_id = re.compile(r'ID=([^;]+)')
    re_parent = re.compile(r'Parent=([^;]+)')
    re_locus = re.compile(r'locus_tag=([^;]+)')
    re_name = re.compile(r'Name=([^;]+)')
    re_gene_attr = re.compile(r'gene=([^;]+)')
    re_protein = re.compile(r'protein_id=([^;]+)')
    re_dbxref_geneid = re.compile(r'GeneID:(\d+)')
    
    try:
        with open(gff_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue
                
                feat_type = parts[2]
                try:
                    start = int(parts[3])
                    end = int(parts[4])
                except ValueError:
                    continue
                    
                attributes = parts[8]
                
                m_id = re_id.search(attributes)
                m_parent = re_parent.search(attributes)
                obj_id = m_id.group(1) if m_id else None
                parent_id = m_parent.group(1) if m_parent else None

                # 1. Harvest Keys from GENE features
                if feat_type == 'gene' and obj_id:
                    if re_name.search(attributes): gene_id_to_keys[obj_id].append(re_name.search(attributes).group(1))
                    if re_locus.search(attributes): gene_id_to_keys[obj_id].append(re_locus.search(attributes).group(1))
                    if re_gene_attr.search(attributes): gene_id_to_keys[obj_id].append(re_gene_attr.search(attributes).group(1))
                    
                    m_geneid = re_dbxref_geneid.search(attributes)
                    if m_geneid: gene_id_to_keys[obj_id].append(f"GeneID:{m_geneid.group(1)}")

                # 2. Harvest Keys from mRNA features
                elif feat_type in ['mRNA', 'transcript', 'mrna'] and obj_id:
                    if parent_id: mrna_id_to_gene_id[obj_id] = parent_id
                    if re_locus.search(attributes): mrna_id_to_keys[obj_id].append(re_locus.search(attributes).group(1))

                # 3. Map CDS -> Protein ID
                elif feat_type == 'CDS' and parent_id:
                    m_prot = re_protein.search(attributes)
                    if m_prot:
                        # Handle multiple parents if strictly necessary, but usually comma separated
                        for pid in parent_id.split(','):
                            cds_parent_to_prot[pid] = m_prot.group(1)

                # 4. Harvest Exons for Structural Fingerprinting
                elif feat_type == 'exon' and parent_id:
                    for pid in parent_id.split(','):
                        mrna_id_to_exons[pid].append((start, end))

        # --- Build Structural Fingerprints ---
        # Map: Fingerprint -> Set of Proteins (to check ambiguity)
        fingerprint_to_prots = defaultdict(set)
        
        for mrna_id, prot_id in cds_parent_to_prot.items():
            if mrna_id in mrna_id_to_exons:
                # Sort exons by start position to ensure consistent string
                exons = sorted(mrna_id_to_exons[mrna_id])
                # Create string representation: "start-end|start-end"
                fingerprint = "|".join([f"{s}-{e}" for s, e in exons])
                fingerprint_to_prots[fingerprint].add(prot_id)

        # Resolve Fingerprint Ambiguity
        for fp, prots in fingerprint_to_prots.items():
            if len(prots) == 1:
                loc_fingerprint_to_prot[fp] = list(prots)[0]

        # --- Aggregate ID Mappings ---
        final_map = {}
        for mrna_id, prot_id in cds_parent_to_prot.items():
            final_map[mrna_id] = prot_id
            if mrna_id.startswith('rna-'):
                final_map[mrna_id.replace('rna-', '')] = prot_id
            
            if mrna_id in mrna_id_to_keys:
                for key in mrna_id_to_keys[mrna_id]:
                    key_to_prots[key].add(prot_id)

            if mrna_id in mrna_id_to_gene_id:
                gene_id = mrna_id_to_gene_id[mrna_id]
                if gene_id in gene_id_to_keys:
                    for key in gene_id_to_keys[gene_id]:
                        key_to_prots[key].add(prot_id)

        # --- Apply Strict Ambiguity Filter to IDs ---
        for key, prots in key_to_prots.items():
            if len(prots) == 1:
                final_map[key] = list(prots)[0]
            
        return final_map, loc_fingerprint_to_prot

    except Exception as e:
        print(f"Error parsing GFF '{gff_path}': {e}")
        return {}, {}

def get_authenticity_score(header):
    """
    Scoring to prioritize Authentic mRNA over misc_RNA.
    """
    header_upper = header.upper()
    if "GBKEY=MRNA" in header_upper or "GBKEY=CDS" in header_upper:
        return 2
    if "XM_" in header or "NM_" in header:
        return 2
    if "GBKEY=MISC_RNA" in header_upper or "XR_" in header:
        return 0
    return 1

def generate_header_fingerprint(header):
    """
    Parses [location=...] from header and generates "start-end|start-end" string.
    Fix: Now handles '<' and '>' characters in NCBI partial coordinates (e.g. 123..>456).
    """
    # Look for [location=...]
    match = re.search(r'\[location=([^\]]+)\]', header)
    if not match:
        return None
    
    loc_str = match.group(1)
    # Remove < and > characters to normalize coordinates
    loc_clean = loc_str.replace('<', '').replace('>', '')
    
    # Extract all number pairs "123..456"
    intervals = re.findall(r'(\d+)\.\.(\d+)', loc_clean)
    
    if not intervals:
        return None
    
    # Convert to ints
    int_intervals = []
    for s, e in intervals:
        int_intervals.append((int(s), int(e)))
    
    # Sort to match GFF order logic
    int_intervals.sort()
    
    # Build string
    return "|".join([f"{s}-{e}" for s, e in int_intervals])

def process_fasta(input_path, output_path, mapping_dict, loc_map, valid_ids_dict, mode='cds', dedup=False, dedup_out_path=None, write_missing=False, missing_out_path=None, track_logic=False):
    """
    Processes headers, buffers matches, and selects the BEST authentic sequence.
    Now supports logic tracking for master table.
    """
    re_protein_tag = re.compile(r'\[protein_id=([^\]]+)\]')
    re_polypeptide = re.compile(r'polypeptide=([^ ]+)')
    re_transcript_id = re.compile(r'\[transcript_id=([^\]]+)\]')
    re_locus_tag = re.compile(r'\[locus_tag=([^\]]+)\]')
    re_gene_tag = re.compile(r'\[gene=([^\]]+)\]')
    re_dbxref_geneid = re.compile(r'GeneID:(\d+)')
    
    # Fallback: Extract ID from NCBI lcl| header (e.g. lcl|NC_123_cds_XP_456.1_789)
    re_lcl_fallback = re.compile(r'_cds_([A-Za-z0-9]+\.\d+)_')

    best_matches = {} # ProteinID -> (score, seq_len, header, full_seq, logic)
    all_candidates_log = defaultdict(list)
    count_total = 0
    
    with open(input_path, 'r', encoding='utf-8', errors='ignore') as fin:
        current_header = None
        current_seq = []
        
        def process_record(header, seq):
            # List of tuples: (candidate_id, logic_string)
            candidates = []
            
            # Priority 0: Direct Explicit IDs
            m_pid = re_protein_tag.search(header)
            if m_pid: candidates.append((m_pid.group(1), "ExplicitID"))
            
            m_poly = re_polypeptide.search(header)
            if m_poly: candidates.append((m_poly.group(1), "ExplicitID"))

            # Priority 1: Structural Fingerprint
            fingerprint = generate_header_fingerprint(header)
            if fingerprint and fingerprint in loc_map:
                candidates.append((loc_map[fingerprint], "Fingerprint"))
            
            # Priority 2: Mappings
            m_trans = re_transcript_id.search(header)
            if m_trans and m_trans.group(1) in mapping_dict:
                candidates.append((mapping_dict[m_trans.group(1)], "TransMap"))

            m_locus = re_locus_tag.search(header)
            if m_locus and m_locus.group(1) in mapping_dict:
                candidates.append((mapping_dict[m_locus.group(1)], "LocusMap"))
            
            m_gene = re_gene_tag.search(header)
            if m_gene and m_gene.group(1) in mapping_dict:
                candidates.append((mapping_dict[m_gene.group(1)], "GeneMap"))
            
            m_geneid = re_dbxref_geneid.search(header)
            if m_geneid:
                key = f"GeneID:{m_geneid.group(1)}"
                if key in mapping_dict:
                    candidates.append((mapping_dict[key], "GeneIDMap"))
                
            # Priority 3: Fallbacks
            simple_id = header.strip().split()[0][1:]
            candidates.append((simple_id, "SimpleID"))
            candidates.append((simple_id + ".p", "SimpleID_p"))
            
            m_lcl = re_lcl_fallback.search(header)
            if m_lcl:
                candidates.append((m_lcl.group(1), "LCLFallback"))
            
            final_id = None
            final_logic = "NA"
            
            # Iterate candidates and pick first valid one
            for cand_id, logic in candidates:
                if cand_id in valid_ids_dict:
                    final_id = cand_id
                    final_logic = logic
                    break
            
            if final_id:
                full_seq = "".join(seq)
                score = get_authenticity_score(header)
                seq_len = len(full_seq)
                
                all_candidates_log[final_id].append(header.strip())
                
                # Logic to keep best score/length
                update = False
                if final_id not in best_matches:
                    update = True
                else:
                    curr_score, curr_len, _, _, _ = best_matches[final_id]
                    if score > curr_score:
                        update = True
                    elif score == curr_score and seq_len > curr_len:
                        update = True
                
                if update:
                    best_matches[final_id] = (score, seq_len, header, full_seq, final_logic)

        for line in fin:
            line = line.strip()
            if not line: continue
            
            if line.startswith('>'):
                if current_header:
                    process_record(current_header, current_seq)
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_header:
            process_record(current_header, current_seq)

    # --- Write Output ---
    metadata = {} # ID -> {'len': int, 'logic': str}

    with open(output_path, 'w') as fout:
        for pid, (_, seq_len, head, seq, logic) in best_matches.items():
            clean_head = head.lstrip('>')
            # fout.write(f"{head}\n{seq}\n")
            if clean_head != pid:
                fout.write(f">{pid} {clean_head}\n{seq}\n")
            else:
                fout.write(f">{pid}\n{seq}\n")
            count_total += 1
            if track_logic:
                metadata[pid] = {'len': seq_len, 'logic': logic}

    # --- Write Dedup Map ---
    if dedup_out_path:
        duplicates_found = {k: v for k, v in all_candidates_log.items() if len(v) > 1}
        if duplicates_found:
            try:
                with open(dedup_out_path, 'w') as flog:
                    flog.write(f"# ProteinID\tKept_Header\tSkipped_Headers...\n")
                    for pid, headers in duplicates_found.items():
                        kept_header = best_matches[pid][2].strip()
                        sorted_headers = [kept_header] + [h for h in headers if h != kept_header]
                        line = f"{pid}\t" + "\t".join(sorted_headers)
                        flog.write(f"{line}\n")
            except Exception as e:
                print(f"Warning: Could not write dedup log to {dedup_out_path}: {e}")

    # --- Write Missing IDs ---
    if write_missing and missing_out_path:
        matched_set = set(best_matches.keys())
        all_prots = set(valid_ids_dict.keys())
        missing_ids = all_prots - matched_set
        
        if missing_ids:
            try:
                with open(missing_out_path, 'w') as fmiss:
                    fmiss.write("# ProteinID\tGeneID\tOriginalProteinHeader\n")
                    for mid in sorted(missing_ids):
                        # Extract header from valid_ids_dict (structure depends on mode)
                        val = valid_ids_dict[mid]
                        if isinstance(val, dict):
                             header = val['header']
                        else:
                             header = val
                             
                        gene_match = re_gene_tag.search(header)
                        gene_id = gene_match.group(1) if gene_match else "NA"
                        fmiss.write(f"{mid}\t{gene_id}\t{header}\n")
            except Exception as e:
                 print(f"Warning: Could not write missing log to {missing_out_path}: {e}")

    return count_total, len(best_matches), metadata

# --- Batch Processing Logic ---

def process_species_task(task_data):
    species, paths, out_dirs, dedup, write_dedup, write_missing, do_master_table = task_data
    
    try:
        # Load proteins (with lengths if master table requested)
        prot_data = get_protein_info(paths['protein'], load_lengths=do_master_table)
        
        # Helper to normalize dict access
        if do_master_table:
            valid_ids = {k: v['header'] for k, v in prot_data.items()}
            input_count = len(prot_data)
        else:
            valid_ids = prot_data
            input_count = len(prot_data)

        if input_count == 0:
            return (f"[{species}] SKIPPED: No valid protein IDs found.", [])

        # Parse GFF
        id_map, loc_map = parse_gff_map(paths['gff'])
        
        results_str = []
        
        # Validation Data Containers
        cds_meta = {}
        cdna_meta = {}

        # 1. Process CDS
        if paths['cds']:
            out_cds = os.path.join(out_dirs['cds'], f"{species}.fna")
            dedup_log_cds = os.path.join(out_dirs['dedup'], f"{species}_cds_dedup.txt") if write_dedup else None
            missing_log_cds = os.path.join(out_dirs['missing'], f"{species}_cds_missing.txt") if write_missing else None
            
            total, unique, cds_meta = process_fasta(
                paths['cds'], out_cds, id_map, loc_map, valid_ids, mode='cds', 
                dedup=True, dedup_out_path=dedup_log_cds,
                write_missing=write_missing, missing_out_path=missing_log_cds,
                track_logic=do_master_table
            )
            results_str.append(f"CDS: {total}")
        else:
            results_str.append("CDS: None")
        
        # 2. Process cDNA
        if paths['cdna']:
            out_cdna = os.path.join(out_dirs['cdna'], f"{species}.fas") 
            dedup_log_cdna = os.path.join(out_dirs['dedup'], f"{species}_cdna_dedup.txt") if write_dedup else None
            missing_log_cdna = os.path.join(out_dirs['missing'], f"{species}_cdna_missing.txt") if write_missing else None

            total, unique, cdna_meta = process_fasta(
                paths['cdna'], out_cdna, id_map, loc_map, valid_ids, mode='cdna', 
                dedup=True, dedup_out_path=dedup_log_cdna,
                write_missing=write_missing, missing_out_path=missing_log_cdna,
                track_logic=do_master_table
            )
            results_str.append(f"cDNA: {total}")
        else:
            results_str.append("cDNA: None")
        
        # 3. Generate Master Table Rows (if requested)
        table_rows = []
        if do_master_table:
            for pid in sorted(prot_data.keys()):
                p_len = prot_data[pid]['len']
                
                # CDS Info
                if pid in cds_meta:
                    c_len = cds_meta[pid]['len']
                    c_log = cds_meta[pid]['logic']
                else:
                    c_len = "NA"
                    c_log = "Missing"
                
                # cDNA Info
                if pid in cdna_meta:
                    d_len = cdna_meta[pid]['len']
                    d_log = cdna_meta[pid]['logic']
                else:
                    d_len = "NA"
                    d_log = "Missing"
                
                # Row: Species, ProteinID, PrtLen, CDSLen, CDNALen, CDS_Logic, cDNA_Logic
                row = f"{species}\t{pid}\t{p_len}\t{c_len}\t{d_len}\t{c_log}\t{d_log}"
                table_rows.append(row)

        return (f"[{species}] Input_Prots: {input_count} | {', '.join(results_str)}", table_rows)
        
    except Exception as e:
        import traceback
        return (f"[{species}] ERROR: {str(e)}\n{traceback.format_exc()}", [])

def match_files(prt_dir, cds_dir, cdna_dir, gff_dir):
    species_tasks = {}
    
    def find_file(directory, species_name, extensions):
        if not directory: return None
        for ext in extensions:
            path = os.path.join(directory, species_name + ext)
            if os.path.exists(path):
                return path
        return None

    print(f"Scanning Protein Directory: {prt_dir}")
    prt_exts = ['.faa', '.fa', '.fasta']
    files = [f for f in os.listdir(prt_dir) if any(f.endswith(e) for e in prt_exts)]
    
    for f in files:
        species_name = os.path.splitext(f)[0]
        p_path = os.path.join(prt_dir, f)
        c_path = find_file(cds_dir, species_name, ['.fna', '.fa', '.fasta', '.cds'])
        d_path = find_file(cdna_dir, species_name, ['.fas', '.fasta', '.fa', '.cdna'])
        g_path = find_file(gff_dir, species_name, ['.gff', '.gff3'])
        
        if g_path and (c_path or d_path):
            species_tasks[species_name] = {
                'protein': p_path,
                'cds': c_path,
                'cdna': d_path,
                'gff': g_path
            }
        else:
            missing = []
            if not c_path and not d_path: missing.append('CDS/cDNA')
            if not g_path: missing.append('GFF')
            print(f"Warning: Skipping {species_name}. Critical files missing: {', '.join(missing)}")
            
    return species_tasks

def resolve_output_path(filename, out_dir):
    if os.path.isabs(filename):
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        return filename

    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        return os.path.join(out_dir, os.path.basename(filename))

    return os.path.abspath(filename)

# --- Main Execution ---

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean sequences using GFF Structural Fingerprinting and ID mapping.")
    
    mode_group = parser.add_argument_group('Mode Selection')
    mode_group.add_argument("--batch", action="store_true", help="Enable batch directory processing mode")
    mode_group.add_argument("--dedup", action="store_true", help="If multiple CDS/cDNA map to one Protein ID, keep only the best one.")
    mode_group.add_argument("--write-dedup", action="store_true", help="Write a log file mapping Protein IDs to the skipped duplicate headers.")
    mode_group.add_argument("--write-missing", action="store_true", help="Write a file listing Protein IDs that were not matched.")
    mode_group.add_argument("--master-table", action="store_true", help="Generate a master_table.tsv with lengths and matching logic.")
    
    single_group = parser.add_argument_group('Single File Mode')
    single_group.add_argument("--gff", help="Path to GFF file")
    single_group.add_argument("--protein", help="Path to CLEANED Protein FASTA")
    single_group.add_argument("--cds", help="Input CDS FASTA")
    single_group.add_argument("--cdna", help="Input cDNA FASTA")
    single_group.add_argument("--out_cds", default="clean_cds.fasta")
    single_group.add_argument("--out_cdna", default="clean_cdna.fasta")
    
    batch_group = parser.add_argument_group('Batch Directory Mode')
    batch_group.add_argument("--prt-dir", help="Directory containing Protein FASTAs")
    batch_group.add_argument("--cds-dir", help="Directory containing CDS FASTAs")
    batch_group.add_argument("--cdna-dir", help="Directory containing cDNA FASTAs")
    batch_group.add_argument("--gff-dir", help="Directory containing GFF files")
    batch_group.add_argument("--out-dir", default=None, help="Root directory for output files")
    batch_group.add_argument("--threads", type=int, default=4, help="Number of parallel threads/processes")
    
    args = parser.parse_args()

    if args.batch:
        if not all([args.prt_dir, args.gff_dir]):
            parser.error("Batch mode requires at least --prt-dir and --gff-dir")
        if not args.cds_dir and not args.cdna_dir:
             parser.error("Batch mode requires at least one of --cds-dir or --cdna-dir")

        cds_out = os.path.join(args.out_dir, "longest_cds")
        cdna_out = os.path.join(args.out_dir, "longest_cdna")
        dedup_out = os.path.join(args.out_dir, "dedup_logs")
        missing_out = os.path.join(args.out_dir, "missing_ids")
        
        os.makedirs(cds_out, exist_ok=True)
        os.makedirs(cdna_out, exist_ok=True)
        if args.write_dedup:
            os.makedirs(dedup_out, exist_ok=True)
        if args.write_missing:
            os.makedirs(missing_out, exist_ok=True)
            
        out_dirs = {'cds': cds_out, 'cdna': cdna_out, 'dedup': dedup_out, 'missing': missing_out}
        
        tasks = match_files(args.prt_dir, args.cds_dir, args.cdna_dir, args.gff_dir)
        print(f"Found {len(tasks)} valid species sets.")
        
        task_list = [(species, paths, out_dirs, args.dedup, args.write_dedup, args.write_missing, args.master_table) for species, paths in tasks.items()]
        
        if not task_list:
            sys.exit(0)

        print(f"Starting processing with {args.threads} threads (Dedup: {args.dedup}, WriteLog: {args.write_dedup}, MasterTable: {args.master_table})...")
        
        all_table_rows = []
        
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            futures = [executor.submit(process_species_task, t) for t in task_list]
            for future in as_completed(futures):
                log_msg, rows = future.result()
                print(log_msg)
                if rows:
                    all_table_rows.extend(rows)
        
        # Write Master Table
        if args.master_table and all_table_rows:
            master_path = os.path.join(args.out_dir, "master_table.tsv")
            print(f"Writing Master Table to: {master_path}")
            try:
                with open(master_path, 'w') as fmaster:
                    fmaster.write("Species\tProteinID\tPrtLen\tCDSLen\tCDNALen\tCDS_Logic\tcDNA_Logic\n")
                    fmaster.write("\n".join(all_table_rows))
                    fmaster.write("\n")
            except Exception as e:
                print(f"Error writing master table: {e}")
                
        print(f"\nBatch processing complete. Outputs in '{args.out_dir}'")
        
    else:
        if not (args.gff and args.protein and (args.cds or args.cdna)):
            parser.error("Single file mode requires --gff, --protein, and at least one of --cds or --cdna")

        print("Running in Single File Mode...")
        # Note: Single file mode master table logic not implemented as per user request focus on Batch
        prot_data = get_protein_info(args.protein, load_lengths=args.master_table)
        
        if args.master_table:
             valid_ids = {k: v['header'] for k, v in prot_data.items()}
        else:
             valid_ids = prot_data
             
        id_map, loc_map = parse_gff_map(args.gff)
        print(f"Input Proteins: {len(valid_ids)}")
        
        # ✅ SAFE handling: out_dir may not exist in single mode
        out_dir = getattr(args, "out_dir", None)

        if args.cds:
            dedup_log = "cds_dedup.txt" if args.write_dedup else None
            missing_log = "cds_missing.txt" if args.write_missing else None

            out_cds_path = resolve_output_path(args.out_cds, out_dir)

            t, u, meta = process_fasta(
                args.cds,
                out_cds_path,
                id_map,
                loc_map,
                valid_ids,
                mode='cds',
                dedup=args.dedup,
                dedup_out_path=dedup_log,
                write_missing=args.write_missing,
                missing_out_path=missing_log,
                track_logic=args.master_table
            )

            print(f"Processed CDS: {t} sequences ({u} unique)")


        if args.cdna:
            dedup_log = "cdna_dedup.txt" if args.write_dedup else None
            missing_log = "cdna_missing.txt" if args.write_missing else None

            out_cdna_path = resolve_output_path(args.out_cdna, out_dir)

            t, u, meta = process_fasta(
                args.cdna,
                out_cdna_path,
                id_map,
                loc_map,
                valid_ids,
                mode='cdna',
                dedup=args.dedup,
                dedup_out_path=dedup_log,
                write_missing=args.write_missing,
                missing_out_path=missing_log,
                track_logic=args.master_table
            )

            print(f"Processed cDNA: {t} sequences ({u} unique)")

   