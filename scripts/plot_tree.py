from Bio import Phylo
from pathlib import Path
import matplotlib.pyplot as plt

def plot_tree(tree_file, save_path=None, figsize=(12, 30)):
    tree = Phylo.read(tree_file, "newick")
    tree.root_at_midpoint()
    plt.rcParams.update({'font.size': 9}) 
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1)
    # Draw tree without showing
    Phylo.draw(tree, do_show=False, axes=ax)
    ax.axis('off')
    plt.tight_layout()

    # Color tip labels containing "BAIT" in red
    for text in ax.texts:
        if "BAIT" in text.get_text():
            text.set_color('red')

    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Saved tree figure to: {save_path}")
    plt.show()