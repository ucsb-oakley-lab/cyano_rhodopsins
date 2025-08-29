from Bio import Phylo
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

def plot_tree(tree_file, color_tsv=None, save_path=None, figsize=(12, 30)):
    """
    Plot a phylogenetic tree and color tip labels according to a 2-column TSV.
    TSV columns: taxon_name<TAB>color
    """
    # Read tree
    tree = Phylo.read(tree_file, "newick")
    tree.root_at_midpoint()
    plt.rcParams.update({'font.size': 9})
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)
    ax.axis('off')
    plt.tight_layout()

    # Read color mapping from TSV if provided
    color_map = {}
    if color_tsv:
        df = pd.read_csv(color_tsv, sep="\t", header=None, names=["taxon", "color"])
        for _, row in df.iterrows():
            color_map[str(row["taxon"]).strip()] = str(row["color"]).strip()

    # Color tip labels
    for text in ax.texts:
        label = text.get_text().strip()
        if label in color_map:
            text.set_color(color_map[label])

    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Saved tree figure to: {save_path}")
    plt.show()