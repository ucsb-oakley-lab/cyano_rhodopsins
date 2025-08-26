#!/usr/bin/env python3
"""
Rename FASTA headers to a clean, analysis-friendly format:

> {PREFIX}__{Sanitized_Taxon}__{SequenceID}

- Taxon is parsed from "Tax=..." if present in the header; otherwise "UnknownTaxon".
- SequenceID defaults to the first whitespace-delimited token (e.g., UniRef90_...).
"""

import argparse
import os
import re
import sys

def sanitize_taxon(tax: str) -> str:
    # Replace any non-alphanumeric/underscore with underscore, collapse repeats, trim
    tax = re.sub(r"[^\w]", "_", tax)
    tax = re.sub(r"_+", "_", tax).strip("_")
    return tax or "UnknownTaxon"

def parse_taxon(header_line: str) -> str:
    # Try to extract Tax=... TaxID=... (UniRef style)
    m = re.search(r"Tax=([^\n]+?)\s+TaxID=", header_line)
    if m:
        return sanitize_taxon(m.group(1))
    # Fallbacks: sometimes there’s just Tax=... at end of line
    m2 = re.search(r"Tax=([^\n]+)$", header_line)
    if m2:
        return sanitize_taxon(m2.group(1))
    return "UnknownTaxon"

def reformat_fasta_headers(input_file: str, output_file: str, prefix: str, force: bool = False) -> int:
    if os.path.exists(output_file) and not force:
        print(f"ERROR: Output file exists: {output_file}. Use --force to overwrite.", file=sys.stderr)
        sys.exit(1)

    n = 0
    with open(input_file) as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                n += 1
                # Sequence ID = first word after '>'
                seq_id = line[1:].split()[0].strip()
                tax = parse_taxon(line)
                new_header = f">{prefix}__{tax}__{seq_id}\n"
                outfile.write(new_header)
            else:
                outfile.write(line)
    return n

def main():
    p = argparse.ArgumentParser(
        description="Rename FASTA headers to: >{PREFIX}__{Sanitized_Taxon}__{SequenceID}"
    )
    p.add_argument("input_fasta", help="Input FASTA file")
    p.add_argument("output_fasta", help="Output FASTA file")
    p.add_argument("--prefix", default="CYANO150", help="Header prefix (default: CYANO150)")
    p.add_argument("--force", action="store_true", help="Overwrite output if it exists")
    args = p.parse_args()

    count = reformat_fasta_headers(args.input_fasta, args.output_fasta, args.prefix, args.force)
    print(f"Reformatted {count} headers → {args.output_fasta}")

if __name__ == "__main__":
    main()
