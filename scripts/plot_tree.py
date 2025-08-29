from Bio import Phylo
from pathlib import Path
import matplotlib.pyplot as plt

def plot_tree(tree_file, save_path=None, figsize=(12, 30), optics_tsv=None, n=20, lmax_table_path=None):
    # Read tree
    tree = Phylo.read(tree_file, "newick")
    tree.root_at_midpoint()
    plt.rcParams.update({'font.size': 9}) 
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)
    ax.axis('off')
    plt.tight_layout()

    # Prepare coloring info
    top_n = set()
    bottom_n = set()
    if optics_tsv:
        import pandas as pd
        df = pd.read_csv(optics_tsv, sep="\t")
        # Drop rows with missing predictions
        df = df.dropna(subset=["Single_Prediction"])
        # Sort by lmax
        df_sorted = df.sort_values("Single_Prediction")
        bottom_n_df = df_sorted.head(n)
        top_n_df = df_sorted.tail(n)
        # Normalize names by stripping whitespace
        bottom_n = set(name.strip() for name in bottom_n_df["Names"])
        top_n = set(name.strip() for name in top_n_df["Names"])

        # Print tables of n lowest and n highest lmax genes
        print(f"\nLowest {n} λmax predictions:")
        print(bottom_n_df[["Names", "Single_Prediction"]].to_string(index=False))
        print(f"\nHighest {n} λmax predictions:")
        print(top_n_df[["Names", "Single_Prediction"]].to_string(index=False))

        # Write table to file if requested
        if lmax_table_path:
            with open(lmax_table_path, "w") as f:
                f.write(f"Lowest {n} λmax predictions:\n")
                f.write(bottom_n_df[["Names", "Single_Prediction"]].to_string(index=False))
                f.write("\n\nHighest {n} λmax predictions:\n")
                f.write(top_n_df[["Names", "Single_Prediction"]].to_string(index=False))
                f.write("\n")
            print(f"Saved λmax table to: {lmax_table_path}")

        # Color tip labels
        for text in ax.texts:
            label = text.get_text().strip()
            if label in bottom_n:
                text.set_color('blue')
            elif label in top_n:
                text.set_color('red')
            elif "BAIT" in label:
                text.set_color('green')
    else:
        # Only color BAIT if no optics_tsv
        for text in ax.texts:
            if "BAIT" in text.get_text():
                text.set_color('green')

    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Saved tree figure to: {save_path}")
    plt.show()