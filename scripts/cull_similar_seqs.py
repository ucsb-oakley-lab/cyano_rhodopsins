#!/usr/bin/env python3
import argparse
import csv
import os
import sys
from collections import OrderedDict
from typing import List, Tuple, Set

from Bio import SeqIO
from Bio.Align import PairwiseAligner

# ---------------------- similarity functions ----------------------

def kmer_set(seq: str, k: int) -> Set[str]:
    if k <= 0 or k > len(seq):
        return set()
    return {seq[i:i+k] for i in range(len(seq) - k + 1)}

def kmer_similarity(a: str, b: str, k: int) -> float:
    """Jaccard-like similarity on k-mer sets, normalized to [0,1]."""
    A = kmer_set(a, k)
    B = kmer_set(b, k)
    if not A or not B:
        return 0.0
    inter = len(A & B)
    denom = min(len(A), len(B))  # conservative normalization
    return inter / denom

def make_aligner(mode: str, match: float, mismatch: float, gap_open: float, gap_extend: float) -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = mode  # 'global' or 'local'
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    return aligner

def align_similarity(a: str, b: str, aligner: PairwiseAligner, match: float) -> float:
    """
    Normalized alignment similarity in [0,1].
    With default params (match=1, mismatch=0, gaps<0), the raw score ~ #matches,
    so normalize by min(len(a), len(b)) * match.
    """
    if not a or not b:
        return 0.0
    score = aligner.align(a, b).score
    denom = max(1.0, min(len(a), len(b)) * (match if match > 0 else 1.0))
    sim = score / denom
    # Clamp to [0,1] for numerical robustness
    if sim < 0: sim = 0.0
    if sim > 1: sim = 1.0
    return sim

# ---------------------- culling logic ----------------------

def compute_max_sims(records: List[SeqIO.SeqRecord],
                     method: str,
                     aligner: PairwiseAligner = None,
                     match: float = 1.0,
                     k: int = 3) -> List[float]:
    """
    For each sequence i, compute max_j!=i similarity(i,j).
    Computes each pair once and updates both i and j.
    """
    n = len(records)
    max_sims = [0.0] * n

    if method == "kmer":
        # Precompute k-mer sets
        ksets = [kmer_set(str(r.seq), k) for r in records]

    for i in range(n):
        if n >= 20 and (i % max(1, n // 10) == 0):
            print(f"Progress: {i}/{n} sequences processed...", file=sys.stderr)
        ai = str(records[i].seq)

        for j in range(i + 1, n):
            bj = str(records[j].seq)
            if method == "align":
                sim = align_similarity(ai, bj, aligner, match)
            else:
                # kmer
                A = ksets[i]; B = ksets[j]
                if not A or not B:
                    sim = 0.0
                else:
                    inter = len(A & B)
                    denom = min(len(A), len(B))
                    sim = inter / denom

            if sim > max_sims[i]:
                max_sims[i] = sim
            if sim > max_sims[j]:
                max_sims[j] = sim

    return max_sims

def cull(records: List[SeqIO.SeqRecord],
         max_sims: List[float],
         discard_k: int = None,
         max_sim_threshold: float = None) -> Tuple[List[SeqIO.SeqRecord], List[SeqIO.SeqRecord]]:
    """
    Either discard the top-K most similar sequences, OR discard all with max_sim >= threshold.
    If both provided, threshold takes precedence.
    """
    indexed = list(enumerate(records))

    if max_sim_threshold is not None:
        to_discard_idx = {i for i, _ in indexed if max_sims[i] >= max_sim_threshold}
    else:
        k = discard_k or 0
        # Sort by max similarity (desc), take first k indices
        order = sorted(range(len(records)), key=lambda i: max_sims[i], reverse=True)
        to_discard_idx = set(order[:k])

    retained = [r for i, r in indexed if i not in to_discard_idx]
    discarded = [r for i, r in indexed if i in to_discard_idx]
    return retained, discarded

# ---------------------- main CLI ----------------------

def main():
    ap = argparse.ArgumentParser(
        description="Cull similar sequences by all-vs-all similarity (alignment or k-mer)."
    )
    ap.add_argument("input_fasta", help="Input FASTA")
    ap.add_argument("retained_fasta", help="Output FASTA for retained sequences")
    ap.add_argument("--discarded-fasta", help="Optional FASTA to write discarded sequences")
    group = ap.add_mutually_exclusive_group()
    group.add_argument("--discard", type=int, default=20,
                       help="Number of most similar sequences to discard (default: 20)")
    group.add_argument("--max-sim", type=float,
                       help="Discard all sequences whose max similarity >= this threshold (0..1)")

    ap.add_argument("--method", choices=["align", "kmer"], default="align",
                    help="Similarity method: 'align' (Biopython PairwiseAligner) or 'kmer' (Jaccard-like) [default: align]")

    # Alignment params
    ap.add_argument("--align-mode", choices=["global", "local"], default="global", help="Alignment mode [default: global]")
    ap.add_argument("--match", type=float, default=1.0, help="Match score [default: 1.0]")
    ap.add_argument("--mismatch", type=float, default=0.0, help="Mismatch score [default: 0.0]")
    ap.add_argument("--gap-open", type=float, default=-1.0, help="Gap open score [default: -1.0]")
    ap.add_argument("--gap-extend", type=float, default=-0.5, help="Gap extend score [default: -0.5]")

    # k-mer params
    ap.add_argument("--k", type=int, default=3, help="k-mer size for kmer method [default: 3]")

    # Outputs / provenance
    ap.add_argument("--summary-csv", help="Write a CSV with per-sequence max similarity")
    ap.add_argument("--force", action="store_true", help="Overwrite output files if they exist")

    args = ap.parse_args()

    if (os.path.exists(args.retained_fasta) or (args.discarded_fasta and os.path.exists(args.discarded_fasta))) and not args.force:
        print("ERROR: Output file exists. Use --force to overwrite.", file=sys.stderr)
        sys.exit(1)

    records = list(SeqIO.parse(args.input_fasta, "fasta"))
    if not records:
        print("No sequences found in input.", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(records)} sequences from {args.input_fasta}", file=sys.stderr)

    # Prepare similarity machinery
    aligner = None
    if args.method == "align":
        aligner = make_aligner(args.align_mode, args.match, args.mismatch, args.gap_open, args.gap_extend)

    # Compute per-sequence max similarity
    max_sims = compute_max_sims(records, args.method, aligner=aligner, match=args.match, k=args.k)

    # Optional summary CSV
    if args.summary_csv:
        with open(args.summary_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["seq_id", "max_similarity"])
            for rec, ms in zip(records, max_sims):
                w.writerow([rec.id, f"{ms:.6f}"])

    # Decide what to discard
    retained, discarded = cull(records, max_sims, discard_k=args.discard, max_sim_threshold=args.max_sim)

    # Write outputs
    SeqIO.write(retained, args.retained_fasta, "fasta")
    if args.discarded_fasta:
        SeqIO.write(discarded, args.discarded_fasta, "fasta")

    print(f"Retained {len(retained)} sequences → {args.retained_fasta}", file=sys.stderr)
    if args.discarded_fasta:
        print(f"Discarded {len(discarded)} sequences → {args.discarded_fasta}", file=sys.stderr)

if __name__ == "__main__":
    main()
