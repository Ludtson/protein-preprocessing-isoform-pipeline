#!/usr/bin/env python3
"""
select_longest_proteins.py

Selects the longest protein isoform per gene from a GFF/FASTA pair.

- Input is a single GFF file and a single protein FASTA file.
- Outputs a new FASTA file with only the longest protein per gene.
- Outputs a single TSV with one row per longest protein.

TSV Columns:
  Species            # Name of the species (from -s or inferred)
  Gene_ID            # ID for the gene
  Transcript_ID      # The ID of the transcript for the selected protein
  Protein_ID         # ID for the selected longest protein
  Protein_Length     # Length of the selected protein
  mRNA_length        # Length of the mRNA feature
  CDS_Length         # Length of the CDS features
  Delta_nt           # Difference between CDS and 3x protein length
  StopFeature        # "1" if a stop codon feature is present, else "0"
  Flags              # Special flags like 'transl_except'
  Gene_Strand        # Strand of the gene (+ or -)
  Gene_Isoforms      # Total count of isoforms for this gene
  mRNA_exons         # Number of exons for the selected isoform
  Other_Isoforms     # Details on all other isoforms of this gene

Other_Isoforms Format:
"T_ID|mRNA=<length>|P_ID|aa=<length>|Stop=<0/1>" (semicolon-separated for multiple)

Usage:
  python select_longest_proteins.py -g <input.gff> -f <input.faa> -v <output.tsv> -o <output.faa>
  
  # With an optional species name
  python select_longest_proteins.py -g <input.gff> -f <input.faa> -v <output.tsv> -o <output.faa> -s "Homo sapiens"

Options:
  -g, --input_gff     Path to the input GFF/GFF3 file (can be .gz).
  -f, --input_fasta   Path to the input protein FASTA file (can be .gz).
  -v, --output_tsv    Path for the output TSV report.
  -o, --output_fasta  Path for the output FASTA file with longest proteins.
  -s, --species       Optional species name to include in reports.
"""

import argparse
import gzip
import os
import re
import sys
from collections import defaultdict
import datetime
from typing import Dict, List, Tuple

__version__ = "1.0.0"
__author__ = "Adekola Owoyemi (Casola Lab, ECCB, Texas A&M University)"

# --------------------------- Helpers -----------------------------------------

def timestamp_message(message: str, file=sys.stderr):
    """Prints a timestamped message to the specified file handle."""
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{current_time}] [INFO] {message}", file=file)

def open_maybe_gzip(path: str, mode: str):
    """Open a file, handling gzipped files transparently."""
    if path.lower().endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", errors="ignore")
    return open(path, mode, encoding="utf-8", errors="ignore")

def parse_attributes(attr_str):
    d = {}
    for item in attr_str.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            d[key] = val
    return d


def format_attributes(attr_dict):
    return ";".join(f"{k}={v}" for k, v in attr_dict.items())


def has_gene_features(gff_lines):
    return any(line.split("\t")[2] == "gene" for line in gff_lines if not line.startswith("#"))


def rebuild_gene_features(gff_lines):
    new_lines = []
    seen_genes = set()

    for line in gff_lines:
        if line.startswith("#"):
            new_lines.append(line)
            continue

        parts = line.strip().split("\t")
        feature_type = parts[2]
        attributes = parse_attributes(parts[8])

        if feature_type == "mRNA":
            gene_id = attributes.get("geneID")

            if gene_id:
                # ✅ create gene feature ONCE
                if gene_id not in seen_genes:
                    gene_line = parts.copy()
                    gene_line[2] = "gene"
                    gene_attrs = {"ID": gene_id}
                    gene_line[8] = format_attributes(gene_attrs)

                    new_lines.append("\t".join(gene_line))
                    seen_genes.add(gene_id)

                # ✅ fix parent relationship
                attributes["Parent"] = gene_id
                parts[8] = format_attributes(attributes)
                line = "\t".join(parts)

        new_lines.append(line)

    return new_lines

def parse_gff_attributes(attr_str):
    """Parse GFF3 column 9 into a dict (handles common escapes)."""
    attrs = {}
    if not attr_str or attr_str.strip() == ".":
        return attrs
    for field in attr_str.strip().split(";"):
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
        elif " " in field:
            k, v = field.split(" ", 1)
        else:
            attrs[field] = True
            continue
        k = k.strip()
        v = v.strip()
        v = v.replace("%2C", ",").replace("%3B", ";").replace("%3D", "=")
        attrs[k] = v
    return attrs

def first_token(s):
    return s.split()[0] if s else s

def normalize_protein_id(raw):
    """Make protein IDs comparable between GFF and FASTA."""
    if raw is None:
        return None
    tok = first_token(raw)
    tok = tok.replace("protein_id:", "")
    tok = tok.replace("UniProtKB:", "")
    return tok

VERSION_SUFFIX_RE = re.compile(r"\.v\d+(?:\.\d+)*$")

def strip_version_suffix(s: str) -> str:
    return VERSION_SUFFIX_RE.sub("", s)

def maybe_from_dbxref(attrs):
    """Try to extract a protein ID from Dbxref-like attributes."""
    dbx = attrs.get("Dbxref") or attrs.get("DbxRef") or ""
    for part in dbx.split(","):
        part = part.strip()
        if not part:
            continue
        if "protein_id:" in part:
            return part.split("protein_id:")[-1]
        if part.startswith(("RefSeq:", "Genbank:", "GenBank:", "UniProtKB:")):
            return part.split(":", 1)[-1]
    return None

def read_gff_maps(gff_path, stderr=sys.stderr):

    transcript_to_gene = {}
    protein_to_transcript = {}
    gene_biotype = {}
    transcript_biotype = {}
    transcripts_with_cds = set()

    tx_span_start = {}
    tx_span_end = {}
    exon_len_sum = defaultdict(int)
    cds_len_sum  = defaultdict(int)
    exon_counts = defaultdict(int)

    mrna_lengths = defaultdict(int)
    cds_lengths = defaultdict(int)

    stopcodon_present = defaultdict(int)
    special_flags = defaultdict(set)
    gene_strand = {}

    transcript_alias_to_id = {}
    derives_protein_to_tx = {}

    transcript_types = {"mRNA", "transcript", "lincRNA", "ncRNA"}
    protein_feature_types = {"polypeptide", "protein", "peptide"}
    biotype_keys_gene = ("gene_biotype", "gene_type", "biotype")
    biotype_keys_tx   = ("transcript_biotype", "transcript_type", "biotype")

    # ✅ Read all lines first
    with open_maybe_gzip(gff_path, "rt") as fh:
        gff_lines = fh.readlines()

    # ✅ Fix AUGUSTUS-style GFF (missing gene features)
    if not has_gene_features(gff_lines):
        timestamp_message(f"No gene features detected in {gff_path} — rebuilding from geneID", stderr)
        gff_lines = rebuild_gene_features(gff_lines)

    # ✅ MAIN PARSING LOOP
    for line in gff_lines:
        if not line or line.startswith("#"):
            continue

        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue

        _seqid, _src, ftype, start, end, _score, strand, _phase, attr = cols
        attrs = parse_gff_attributes(attr)

        # ✅ Capture special flags
        for key in ("transl_except", "selenocysteine", "ribosomal_slippage", "exception"):
            if key in attrs and "Parent" in attrs:
                for tid0 in attrs["Parent"].split(","):
                    special_flags[tid0].add(key)

        # ✅ GENE
        if ftype == "gene":
            gid = attrs.get("ID") or attrs.get("Name")
            if gid:
                gene_strand[gid] = strand
                for k in biotype_keys_gene:
                    if k in attrs:
                        gene_biotype[gid] = attrs[k]
                        break

        # ✅ TRANSCRIPT
        elif ftype in transcript_types:
            tid = attrs.get("ID") or attrs.get("transcript_id") or attrs.get("Name")
            parent = attrs.get("Parent")
            if not tid or not parent:
                continue

            gid = parent.split(",")[0]
            transcript_to_gene[tid] = gid

            for k in biotype_keys_tx:
                if k in attrs:
                    transcript_biotype[tid] = attrs[k]
                    break

            try:
                s = int(start); e = int(end)
                if tid not in tx_span_start or s < tx_span_start[tid]:
                    tx_span_start[tid] = s
                if tid not in tx_span_end or e > tx_span_end[tid]:
                    tx_span_end[tid] = e
            except ValueError:
                pass

            aliases = set()
            name = attrs.get("Name")
            if name:
                aliases.add(name)
            aliases.add(strip_version_suffix(tid))
            aliases.add(tid)

            for a in aliases:
                transcript_alias_to_id.setdefault(a, tid)

        # ✅ EXONS
        elif ftype == "exon":
            parent = attrs.get("Parent")
            if not parent:
                continue

            try:
                s = int(start); e = int(end)
                length = e - s + 1
            except ValueError:
                continue

            for tid in parent.split(","):
                exon_len_sum[tid] += max(0, length)
                exon_counts[tid] += 1

        # ✅ CDS
        elif ftype == "CDS":
            parent = attrs.get("Parent")
            if not parent:
                continue

            try:
                s = int(start); e = int(end)
                length = e - s + 1
            except ValueError:
                length = None

            for tid in parent.split(","):
                transcripts_with_cds.add(tid)
                if length is not None:
                    cds_len_sum[tid] += max(0, length)

                pid = (
                    attrs.get("protein_id")
                    or attrs.get("proteinId")
                    or maybe_from_dbxref(attrs)
                )
                pid = normalize_protein_id(pid) if pid else None

                if pid:
                    protein_to_transcript.setdefault(pid, tid)

        # ✅ STOP CODON
        elif ftype == "stop_codon":
            parent = attrs.get("Parent")
            if parent:
                for tid in parent.split(","):
                    stopcodon_present[tid] = 1

        # ✅ PROTEIN FEATURES
        elif ftype in protein_feature_types:
            pid = attrs.get("ID") or attrs.get("Name")
            tid = attrs.get("Derives_from") or attrs.get("Parent")
            pid = normalize_protein_id(pid) if pid else None
            if pid and tid:
                derives_protein_to_tx[pid] = tid

    # ✅ Merge derived protein mappings
    for pid, tid in derives_protein_to_tx.items():
        protein_to_transcript.setdefault(pid, tid)

    # ✅ mRNA lengths
    mrna_lengths = {}
    all_tids = set(transcript_to_gene) | set(exon_len_sum) | set(tx_span_start)
    for tid in all_tids:
        if exon_len_sum.get(tid, 0) > 0:
            mrna_lengths[tid] = exon_len_sum[tid]
        elif tid in tx_span_start and tid in tx_span_end:
            mrna_lengths[tid] = tx_span_end[tid] - tx_span_start[tid] + 1
        else:
            mrna_lengths[tid] = 0

    # ✅ CDS lengths
    cds_lengths = {}
    for tid in set(transcript_to_gene) | set(cds_len_sum):
        cds_lengths[tid] = cds_len_sum.get(tid, 0)

    return (
        transcript_to_gene,
        protein_to_transcript,
        gene_biotype,
        transcript_biotype,
        transcripts_with_cds,
        mrna_lengths,
        cds_lengths,
        stopcodon_present,
        special_flags,
        transcript_alias_to_id,
        exon_counts,
        gene_strand
    )

def is_protein_coding(gene_id, transcript_id, gene_biotype, transcript_biotype, transcripts_with_cds):
    tx_bt = (transcript_biotype.get(transcript_id) or "").lower()
    if tx_bt == "protein_coding": return True
    if transcript_id in transcripts_with_cds: return True
    gn_bt = (gene_biotype.get(gene_id) or "").lower()
    return gn_bt == "protein_coding"

def read_fasta_sequences_and_tags(fasta_path):
    """
    Read a protein FASTA and return:
      - seqs: dict primary_id -> sequence  (primary_id = first token, normalized)
      - transcript_tag: dict primary_id -> transcript tag value if header contains 'transcript=...'
    """
    seqs = {}
    transcript_tag = {}
    cur_id = None
    chunks = []
    with open_maybe_gzip(fasta_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(chunks)
                header = line[1:].strip()
                primary = first_token(header)
                pid = normalize_protein_id(primary)
                cur_id = pid
                chunks = []
                m = re.search(r"(?:^|\s)transcript=([^\s;]+)", header)
                if m:
                    transcript_tag[cur_id] = m.group(1)
            else:
                s = line.strip()
                if cur_id is not None and s:
                    chunks.append(s)
    if cur_id is not None:
        seqs[cur_id] = "".join(chunks)
    return seqs, transcript_tag

def main():
    parser = argparse.ArgumentParser(
        description="Select the longest protein isoform per gene from GFF and protein FASTA input."
    )
    parser.add_argument("-g", "--input_gff", required=True, help="Input GFF file.")
    parser.add_argument("-f", "--input_fasta", required=True, help="Input FASTA file.")
    parser.add_argument("-v", "--output_tsv", required=True, help="Output TSV file.")
    parser.add_argument("-o", "--output_fasta", required=True, help="Output FASTA file.")
    parser.add_argument("-s", "--species", help="Optional species name.")
    args = parser.parse_args()

    if args.species:
        report_name = args.species
    else:
        report_name = os.path.basename(args.input_gff).split('.')[0]
    
    timestamp_message(f"Starting isoform selection for '{report_name}'.")

    try:
        (transcript_to_gene, protein_to_transcript, gene_biotype, transcript_biotype, transcripts_with_cds,
         mrna_lengths, cds_lengths, stopcodon_present, special_flags, transcript_alias_to_id, exon_counts, gene_strand) = read_gff_maps(args.input_gff)
        protein_seqs, fasta_transcript_tag = read_fasta_sequences_and_tags(args.input_fasta)
    except Exception as e:
        timestamp_message(f"[ERROR] Failed to parse input files: {e}")
        sys.exit(1)
        
    protein_lengths = {pid: len(seq) for pid, seq in protein_seqs.items()}
    n_input_seqs = len(protein_seqs)

    all_isoforms_details = defaultdict(list)
    used_pids = set()

    # 1) Primary: protein_id on CDS / polypeptide-protein features
    for pid, tid in protein_to_transcript.items():
        gid = transcript_to_gene.get(tid)
        if not gid: continue
        if not is_protein_coding(gid, tid, gene_biotype, transcript_biotype, transcripts_with_cds): continue
        plen = protein_lengths.get(pid)
        if plen is None or plen == 0: continue
        all_isoforms_details[gid].append({
            "gene_id": gid, "transcript_id": tid, "protein_id": pid, "protein_length": plen,
            "mRNA_length": mrna_lengths.get(tid, 0), "CDS_length": cds_lengths.get(tid, 0),
            "stop_feature": stopcodon_present.get(tid, 0), "flags": ";".join(special_flags.get(tid, [])),
            "exon_count": exon_counts.get(tid, 0), "gene_strand": gene_strand.get(gid, ".")
        })
        used_pids.add(pid)

    # 2) Fallback A: transcript-as-protein via alias map and/or transcript= tag
    for primary_pid, seq in protein_seqs.items():
        if primary_pid in used_pids: continue
        t_tag = fasta_transcript_tag.get(primary_pid)
        tid = None
        if t_tag and t_tag in transcript_alias_to_id: tid = transcript_alias_to_id[t_tag]
        elif primary_pid in transcript_alias_to_id: tid = transcript_alias_to_id[primary_pid]
        
        if tid:
            gid = transcript_to_gene.get(tid)
            if not gid: continue
            if not is_protein_coding(gid, tid, gene_biotype, transcript_biotype, transcripts_with_cds): continue
            plen = len(seq)
            if plen == 0: continue
            all_isoforms_details[gid].append({
                "gene_id": gid, "transcript_id": tid, "protein_id": primary_pid, "protein_length": plen,
                "mRNA_length": mrna_lengths.get(tid, 0), "CDS_length": cds_lengths.get(tid, 0),
                "stop_feature": stopcodon_present.get(tid, 0), "flags": ";".join(special_flags.get(tid, [])),
                "exon_count": exon_counts.get(tid, 0), "gene_strand": gene_strand.get(gid, ".")
            })
            used_pids.add(primary_pid)

    # 3) Fallback B: gene-as-protein (Maker-style)
    gene_ids = set(transcript_to_gene.values())
    tx_by_gene = defaultdict(list)
    for tid, gid in transcript_to_gene.items():
        tx_by_gene[gid].append(tid)

    for primary_pid, seq in protein_seqs.items():
        if primary_pid in used_pids: continue
        if primary_pid in gene_ids:
            gid = primary_pid
            if not tx_by_gene.get(gid): continue
            tid = tx_by_gene[gid][0]
            if not is_protein_coding(gid, tid, gene_biotype, transcript_biotype, transcripts_with_cds): continue
            plen = len(seq)
            if plen == 0: continue
            all_isoforms_details[gid].append({
                "gene_id": gid, "transcript_id": tid, "protein_id": primary_pid, "protein_length": plen,
                "mRNA_length": mrna_lengths.get(tid, 0), "CDS_length": cds_lengths.get(tid, 0),
                "stop_feature": stopcodon_present.get(tid, 0), "flags": ";".join(special_flags.get(tid, [])),
                "exon_count": exon_counts.get(tid, 0), "gene_strand": gene_strand.get(gid, ".")
            })
            used_pids.add(primary_pid)

    n_genes_processed = len(all_isoforms_details)
    n_output_seqs = 0
    longest_proteins_to_write = {}
    tsv_rows = []

    for gene_id, isoforms in all_isoforms_details.items():
        if not isoforms: continue
        
        longest_isoform = max(isoforms, key=lambda x: x["protein_length"])

        if longest_isoform["protein_length"] > 0:
            longest_proteins_to_write[gene_id] = {
                "protein_id": longest_isoform["protein_id"],
                "sequence": protein_seqs.get(longest_isoform["protein_id"], ""),
                "length": longest_isoform["protein_length"],
            }
            n_output_seqs += 1
            
            other_isoforms = [
                iso for iso in isoforms if iso["protein_id"] != longest_isoform["protein_id"]
            ]
            
            other_iso_str = ""
            if other_isoforms:
                other_iso_list = [
                    f"{iso['transcript_id']}|mRNA={iso['mRNA_length']}|{iso['protein_id']}|aa={iso['protein_length']}|Stop={iso['stop_feature']}"
                    for iso in other_isoforms
                    ]

                other_iso_str = ";".join(other_iso_list)
            
            delta_nt = longest_isoform["CDS_length"] - (3 * longest_isoform["protein_length"]) if longest_isoform["protein_length"] > 0 else "."
            
            tsv_rows.append([
                report_name,
                gene_id,
                longest_isoform["transcript_id"],
                longest_isoform["protein_id"],
                str(longest_isoform["protein_length"]),
                str(longest_isoform["mRNA_length"]),
                str(longest_isoform["CDS_length"]),
                str(delta_nt),
                str(longest_isoform["stop_feature"]),
                longest_isoform["flags"] if longest_isoform["flags"] else ".",
                longest_isoform["gene_strand"],
                str(len(isoforms)),
                str(longest_isoform["exon_count"]),
                other_iso_str
            ])

    try:
        with open_maybe_gzip(args.output_fasta, "wt") as oh:
            for gene_id, data in longest_proteins_to_write.items():
                header = f">{data['protein_id']}"
                oh.write(f"{header}\n{data['sequence']}\n")

        with open_maybe_gzip(args.output_tsv, "wt") as oh:
            headers = ["Species", "Gene_ID", "Transcript_ID", "Protein_ID", "Protein_Length", "mRNA_length", "CDS_Length", "Delta_nt", "StopFeature", "Flags", "Gene_Strand", "Gene_Isoforms", "mRNA_exons", "Other_Isoforms_Details"]
            oh.write("\t".join(headers) + "\n")
            for row in tsv_rows:
                oh.write("\t".join(row) + "\n")
    except Exception as e:
        timestamp_message(f"[ERROR] Failed to write output files: {e}")
        sys.exit(1)

    summary_header = ["File/Species", "N_input_proteins", "N_genes_processed", "N_output_proteins"]
    summary_data = [report_name, str(n_input_seqs), str(n_genes_processed), str(n_output_seqs)]
    print("\t".join(summary_header))
    print("\t".join(summary_data))

    timestamp_message(f"Processed {n_input_seqs} input protein sequences.")
    timestamp_message(f"Selected {n_output_seqs} longest isoforms from {n_genes_processed} genes.")
    timestamp_message(f"Output FASTA written to '{args.output_fasta}'.")
    timestamp_message(f"Output TSV written to '{args.output_tsv}'.\n")

if __name__ == "__main__":
    main()
