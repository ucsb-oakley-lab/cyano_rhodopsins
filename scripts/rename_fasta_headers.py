#!/usr/bin/env python3
"""
Rename FASTA headers to: >{PREFIX}__{Sanitized_Taxon}__{SequenceID}
…and clean accidental metadata lines inside sequences (e.g., 'RepID=...').

- Taxon is parsed from 'Tax=...' if present; else 'UnknownTaxon'.
- SequenceID = first whitespace-delimited token after '>'.
- Non-header lines that contain metadata (e.g., 'RepID=') are dropped,
  but if they accidentally include sequence, we salvage from the first 'M' onward.
"""

import argparse
import os
import re
import sys
from textwrap import wrap

AA_RE = re.compile(r'^[A-Za-z*.\-]+$')  # allowed sequence characters
FIND_SEQ_FROM_M = re.compile(r'M[A-Za-z*.\-]+')  # typical start codon onward

def sanitize_taxon(tax: str) -> str:
    tax = re.sub(r"[^\w]", "_", tax)
    tax = re.sub(r"_+", "_", tax).strip("_")
    return tax or "UnknownTaxon"

def parse_taxon(header_line: str) -> str:
    m = re.search(r"Tax=([^\n]+?)\s+TaxID=", header_line)
    if m:
        return sanitize_taxon(m.group(1))
    m2 = re.search(r"Tax=([^\n]+)$", header_line)
    if m2:
        return sanitize_taxon(m2.group(1))
    return "UnknownTaxon"

def is_pure_seq_line(s: str) -> bool:
    s = s.strip()
    return bool(s) and bool(AA_RE.fullmatch(s))

def salvage_seq_if_glued(s: str) -> str:
    """
    If a metadata line (e.g., 'RepID=...') accidentally contains sequence
    glued to its end, try to salvage from the first 'M'.
    """
    s = s.strip()
    m = FIND_SEQ_FROM_M.search(s)
    return m.group(0) if m else ""

def reformat_fasta_headers(input_file: str, output_file: str, prefix: str,
                           force: bool = False, wrap_len: int | None = 60) -> int:
    if os.path.exists(output_file) and not force:
        print(f"ERROR: Output file exists: {output_file}. Use --force to overwrite.", file=sys.stderr)
        sys.exit(1)

    n_headers = 0
    dropped_meta = 0
    salvaged = 0
    filtered_chars = 0

    out = open(output_file, "w")
    with open(input_file) as infile, out:
        for line in infile:
            if line.startswith(">"):
                n_headers += 1
                seq_id = line[1:].split()[0].strip()
                tax = parse_taxon(line)
                out.write(f">{prefix}__{tax}__{seq_id}\n")
                continue

            # Non-header lines: try to keep only real sequence
            s = line.strip()
            if not s:
                continue

            if is_pure_seq_line(s):
                seq = s.upper()
            else:
                # Likely a metadata line (contains '=' or spaces, etc.)
                maybe = salvage_seq_if_glued(s)
                if maybe:
                    salvaged += 1
                    seq = maybe.upper()
                else:
                    dropped_meta += 1
                    continue

            # As an extra guard, strip any accidental non-AA chars
            cleaned = re.sub(r'[^A-Za-z*.\-]', '', seq)
            filtered_chars += len(seq) - len(cleaned)
            if not cleaned:
                continue

            if wrap_len and wrap_len > 0:
                for chunk in wrap(cleaned, wrap_len):
                    out.write(chunk + "\n")
            else:
                out.write(cleaned + "\n")

    print(f"Reformatted {n_headers} headers → {output_file}")
    if dropped_meta or salvaged or filtered_chars:
        print(f"[cleaning] dropped_meta_lines={dropped_meta}, salvaged_seq_fragments={salvaged}, filtered_chars={filtered_chars}")
    return n_headers

def main():
    p = argparse.ArgumentParser(
        description="Rename FASTA headers to >{PREFIX}__{Sanitized_Taxon}__{SequenceID} and clean metadata lines."
    )
    p.add_argument("input_fasta", help="Input FASTA file")
    p.add_argument("output_fasta", help="Output FASTA file")
    p.add_argument("--prefix", default="CYANO150", help="Header prefix (default: CYANO150)")
    p.add_argument("--force", action="store_true", help="Overwrite output if it exists")
    p.add_argument("--wrap", type=int, default=60, help="Wrap sequence to this many chars/line (0 to disable)")
    args = p.parse_args()

    wrap_len = None if args.wrap == 0 else args.wrap
    reformat_fasta_headers(args.input_fasta, args.output_fasta, args.prefix, args.force, wrap_len)

if __name__ == "__main__":
    main()
