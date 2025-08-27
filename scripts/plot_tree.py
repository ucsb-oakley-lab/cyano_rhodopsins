from Bio import Phylo
from pathlib import Path
import matplotlib.pyplot as plt  # <-- Make sure this is imported

def plot_tree(tree_file, save_path=None, figsize=(12, 30)):
    tree = Phylo.read(tree_file, "newick")
    tree.root_at_midpoint()
    plt.rcParams.update({'font.size': 9}) 
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)
    ax.axis('off')  # <-- This hides the bounding box and axes
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Saved tree figure to: {save_path}")
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Midpoint root and plot a Newick tree")
    parser.add_argument("tree_file", help="Path to Newick tree file")
    parser.add_argument("--save", help="Path to save the figure (PDF/PNG)")
    parser.add_argument("--height", type=int, default=30, help="Figure height (inches)")
    parser.add_argument("--width", type=int, default=12, help="Figure width (inches)")
    args = parser.parse_args()
    plot_tree(
        args.tree_file,
        save_path=args.save,
        figsize=(args.width, args.height)
    )