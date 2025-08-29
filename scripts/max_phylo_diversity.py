import argparse
from Bio import Phylo
import itertools

def parse_args():
    parser = argparse.ArgumentParser(description="Select y taxa from x to maximize phylogenetic diversity.")
    parser.add_argument("--tree", required=True, help="Input tree file (Newick format)")
    parser.add_argument("--taxa", required=True, help="File with taxon names (one per line)")
    parser.add_argument("-y", type=int, required=True, help="Number of taxa to select")
    return parser.parse_args()

def pairwise_distance(tree, taxa):
    # Compute all pairwise distances between taxa
    dists = {}
    for t1, t2 in itertools.combinations(taxa, 2):
        try:
            d = tree.distance(t1, t2)
            dists[(t1, t2)] = d
        except Exception:
            pass  # skip if taxa not found
    return dists

def greedy_max_diversity(tree, taxa, y):
    # Start with the pair with the largest distance
    dists = pairwise_distance(tree, taxa)
    if not dists:
        raise ValueError("No valid pairwise distances found.")
    (t1, t2), _ = max(dists.items(), key=lambda x: x[1])
    selected = {t1, t2}
    while len(selected) < y:
        best_taxon = None
        best_score = -1
        for t in set(taxa) - selected:
            score = sum(tree.distance(t, s) for s in selected)
            if score > best_score:
                best_score = score
                best_taxon = t
        selected.add(best_taxon)
    return list(selected)

def main():
    args = parse_args()
    tree = Phylo.read(args.tree, "newick")
    with open(args.taxa) as f:
        taxa = [line.strip() for line in f if line.strip()]
    if args.y > len(taxa):
        raise ValueError("y cannot be greater than the number of taxa provided.")
    selected = greedy_max_diversity(tree, taxa, args.y)
    print("Selected taxa maximizing phylogenetic diversity:")
    for t in selected:
        print(t)

if __name__ == "__main__":
    main()