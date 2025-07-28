from Bio import SeqIO
from Bio.Align import PairwiseAligner
import argparse

def max_similarity(seq, others, aligner):
    max_score = 0
    for other in others:
        if seq.id == other.id:
            continue
        score = aligner.align(seq.seq, other.seq).score
        if score > max_score:
            max_score = score
    return max_score

def discard_most_similar(input_fasta, output_fasta, num_to_discard=20):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5

    print("Computing max similarity scores...")
    scores = []
    for i, seq in enumerate(records):
        others = records[:i] + records[i+1:]
        max_score = max_similarity(seq, others, aligner)
        scores.append((seq, max_score))

    print("Selecting sequences to discard...")
    scores.sort(key=lambda x: x[1], reverse=True)
    to_discard = set(r[0].id for r in scores[:num_to_discard])

    retained = [r for r in records if r.id not in to_discard]
    SeqIO.write(retained, output_fasta, "fasta")
    print(f"Done. Retained {len(retained)} sequences and wrote to {output_fasta}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta", help="Input FASTA file with opsin sequences")
    parser.add_argument("output_fasta", help="Output FASTA file (after discarding similar ones)")
    parser.add_argument("--discard", type=int, default=20, help="Number of most similar sequences to discard")
    args = parser.parse_args()

    discard_most_similar(args.input_fasta, args.output_fasta, args.discard)

