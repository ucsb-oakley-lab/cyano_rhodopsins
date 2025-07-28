import argparse
from Bio import Phylo
import matplotlib.pyplot as plt

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Midpoint root and plot a Newick tree")
parser.add_argument("tree_file", help="Path to Newick tree file")
args = parser.parse_args()

# Read and midpoint-root the tree
tree = Phylo.read(args.tree_file, "newick")
tree.root_at_midpoint()

plt.rcParams.update({'font.size': 6}) 

# Plot the tree
fig = plt.figure(figsize=(12, 20))
Phylo.draw(tree, do_show=False)
plt.tight_layout()
plt.show()

