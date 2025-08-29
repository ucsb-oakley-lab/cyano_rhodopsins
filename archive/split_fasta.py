#!/usr/bin/env python3

import argparse
import os

def parse_fasta(filepath):
    with open(filepath, 'r') as f:
        sequences = []
        current_seq = []
        for line in f:
            if line.startswith(">"):
                if current_seq:
                    sequences.append("".join(current_seq))
                    current_seq = []
                current_seq.append(line)
            else:
                current_seq.append(line)
        if current_seq:
            sequences.append("".join(current_seq))
    return sequences

def split_list(lst, n):
    """Split list into n parts as evenly as possible"""
    k, m = divmod(len(lst), n)
    return [lst[i*k + min(i, m):(i+1)*k + min(i+1, m)] for i in range(n)]

def write_splits(sequences, base, ext):
    parts = split_list(sequences, 4)
    for i, chunk in enumerate(parts, 1):
        out_file = f"{base}{i}{ext}"
        with open(out_file, 'w') as f:
            f.writelines(chunk)
        print(f"Wrote {len(chunk)} sequences to {out_file}")

def split_filename(filename):
    base, ext = os.path.splitext(filename)
    if ext.lower() in ['.gz', '.bz2', '.xz']:
        base, ext2 = os.path.splitext(base)
        ext = ext2 + ext
    return base, ext

def main():
    parser = argparse.ArgumentParser(description="Split a FASTA file into 4 roughly equal parts.")
    parser.add_argument("input_file", help="Input FASTA file")
    args = parser.parse_args()

    base, ext = split_filename(args.input_file)
    sequences = parse_fasta(args.input_file)
    write_splits(sequences, base, ext)

if __name__ == "__main__":
    main()

