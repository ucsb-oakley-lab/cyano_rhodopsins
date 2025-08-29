#!/usr/bin/env python3

import argparse

def load_gene_names(filename):
    with open(filename) as f:
        return set(line.strip() for line in f if line.strip())

def filter_fasta_by_names(fasta_file, gene_names, output_file):
    with open(fasta_file) as f, open(output_file, 'w') as out:
        write = False
        for line in f:
            if line.startswith('>'):
                # Extract the first word after '>' to match gene name
                gene_id = line[1:].strip().split()[0]
                write = gene_id in gene_names
            if write:
                out.write(line)

def main():
    parser = argparse.ArgumentParser(description="Extract sequences from a FASTA file by gene names.")
    parser.add_argument("fasta_file", help="Input FASTA file")
    parser.add_argument("gene_list", help="File with gene names (one per line)")
    parser.add_argument("-o", "--output", default="filtered_output.fasta", help="Output FASTA file name")
    args = parser.parse_args()

    gene_names = load_gene_names(args.gene_list)
    filter_fasta_by_names(args.fasta_file, gene_names, args.output)
    print(f"Filtered sequences written to {args.output}")

if __name__ == "__main__":
    main()

