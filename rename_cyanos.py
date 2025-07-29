#!/usr/bin/env python3

import argparse
import re

def reformat_fasta_headers(input_file, output_file):
    with open(input_file) as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Extract UniRef90 ID (first word in header)
                uniref_id = line[1:].split()[0]

                # Extract Tax= value using regex
                tax_match = re.search(r'Tax=([^\n]+?)\s+TaxID=', line)
                if tax_match:
                    tax_name = tax_match.group(1)
                    tax_name = tax_name.replace(' ', '_')
                    tax_name = tax_name.replace('.', '').replace('-', '')
                else:
                    tax_name = "UnknownTaxon"

                new_header = f">CYANO150__{tax_name}__{uniref_id}\n"
                outfile.write(new_header)
            else:
                outfile.write(line)

def main():
    parser = argparse.ArgumentParser(description="Rename FASTA headers with CYANO150 format and cleaned taxon names.")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("output_fasta", help="Output FASTA file with modified headers")
    args = parser.parse_args()

    reformat_fasta_headers(args.input_fasta, args.output_fasta)
    print(f"Reformatted FASTA written to {args.output_fasta}")

if __name__ == "__main__":
    main()

