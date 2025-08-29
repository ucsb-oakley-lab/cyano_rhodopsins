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
import re
def _sanitize_name(name):
    # Remove brackets and their contents, and strip whitespace
    return re.sub(r'[\[\]]', '', name).strip()

def _find_matching_tip(name, tree_tips):
    # Try exact match first
    if name in tree_tips:
        return name
    # Try sanitized match
    sanitized = _sanitize_name(name)
    if sanitized in tree_tips:
        return sanitized
    # Try fuzzy match: ignore case and strip all non-alphanum
    def normalize(s):
        return re.sub(r'[^A-Za-z0-9]', '', s).lower()
    norm_name = normalize(name)
    for tip in tree_tips:
        if normalize(tip) == norm_name:
            return tip
    return None

def greedy_max_diversity(tree, taxa, y):
    """
    Greedy algorithm to select y taxa maximizing phylogenetic diversity.
    tree: Bio.Phylo tree object
    taxa: list of taxon names to select from
    y: number of taxa to select
    """
    if y > len(taxa):
        raise ValueError("y cannot be greater than the number of taxa")
    tree_tips = {term.name for term in tree.get_terminals()}
    # Map input taxa to actual tree tip names (with fuzzy/sanitized matching)
    mapped_taxa = []
    for t in taxa:
        match = _find_matching_tip(t, tree_tips)
        if match:
            mapped_taxa.append(match)
        else:
            print(f"Warning: taxon '{t}' not found in tree tips. Skipping.")
    if not mapped_taxa:
        raise ValueError("No valid taxa found in tree.")
    selected = set()
    # Start with a random taxon
    selected.add(mapped_taxa[0])
    while len(selected) < min(y, len(mapped_taxa)):
        best_score = -1
        best_taxon = None
        for t in set(mapped_taxa) - selected:
            try:
                score = sum(tree.distance(t, s) for s in selected)
            except Exception as e:
                print(f"Warning: could not compute distance for '{t}': {e}")
                continue
            if score > best_score:
                best_score = score
                best_taxon = t
        if best_taxon is None:
            break
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