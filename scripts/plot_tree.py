#!/usr/bin/env python3
import argparse, re
from pathlib import Path
from Bio import Phylo
import matplotlib.pyplot as plt

def extract_label(name: str, mode: str) -> str:
    if not name:
        return ""
    if mode == "full":
        return name
    if mode == "taxon":
        parts = name.split("__")  # e.g. CYANO150__Taxon__UniRef90_...
        if len(parts) >= 3:
            return parts[1]
        m = re.search(r"Tax=([^ ]+)", name)  # fallback for old headers
        return m.group(1).replace(" ", "_") if m else name
    if mode == "id":
        m = re.search(r"(UniRef90_[^\s\)]+)", name)
        return m.group(1) if m else name
    return name

def plot_tree(tree_path, fmt="newick", midpoint=True, outgroup=None, ladderize=True,
              label_mode="taxon", label_internal=False, show_support=False, support_min=70.0,
              width=None, height=None, per_tip=0.20, font_size=6, output=None, dpi=300, no_show=False):
    tree = Phylo.read(tree_path, fmt)

    if outgroup:
        # try substring match first, then exact
        candidates = [c for c in tree.find_clades() if c.name and outgroup in c.name]
        if candidates:
            tree.root_with_outgroup(candidates[0])
        else:
            term = tree.find_any(name=outgroup)
            if term:
                tree.root_with_outgroup(term)
            elif midpoint:
                print(f"[warn] Outgroup '{outgroup}' not found; using midpoint.")
                tree.root_at_midpoint()
    elif midpoint:
        tree.root_at_midpoint()

    if ladderize:
        tree.ladderize()

    def label_func(clade):
        if clade.is_terminal() or label_internal:
            return extract_label(clade.name, label_mode)
        return None

    def support_func(clade):
        # IQ-TREE stores supports in clade.confidence
        val = getattr(clade, "confidence", None)
        try:
            if val is not None and float(val) >= support_min:
                return f"{int(round(float(val)))}"
        except Exception:
            pass
        return None

    n = tree.count_terminals()
    fig_w = width or 12
    fig_h = height or max(6, min(30, n * per_tip))

    plt.rcParams.update({'font.size': font_size})
    fig = plt.figure(figsize=(fig_w, fig_h))
    Phylo.draw(
        tree,
        do_show=False,
        label_func=label_func,
        branch_labels=(support_func if show_support else None)
    )
    plt.tight_layout()

    if output:
        out = Path(output)
        out.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"Saved tree to {out.resolve()}")
    if not no_show:
        plt.show()

def main():
    p = argparse.ArgumentParser(description="Root and plot a tree with nicer labels.")
    p.add_argument("tree_file", help="Path to tree (.treefile/.contree/.nwk)")
    p.add_argument("--format", default="newick", choices=["newick","nexus","phyloxml"])
    p.add_argument("--midpoint", action="store_true", help="Midpoint root (default if no --outgroup)")
    p.add_argument("--outgroup", help="Outgroup name (substring or exact tip)")
    p.add_argument("--ladderize", action="store_true", help="Ladderize the tree")
    p.add_argument("--label-mode", default="taxon", choices=["taxon","id","full"])
    p.add_argument("--label-internal", action="store_true", help="Label internal nodes too")
    p.add_argument("--show-support", action="store_true", help="Draw node supports")
    p.add_argument("--support-min", type=float, default=70.0)
    p.add_argument("--width", type=float)
    p.add_argument("--height", type=float)
    p.add_argument("--per-tip", type=float, default=0.20, help="inches per tip for auto height")
    p.add_argument("--font-size", type=int, default=6)
    p.add_argument("--output", "-o", help="Output image (.png/.pdf/.svg)")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--no-show", action="store_true")
    args = p.parse_args()

    # default to midpoint if neither specified
    if not args.midpoint and not args.outgroup:
        args.midpoint = True

    plot_tree(
        args.tree_file, args.format, args.midpoint, args.outgroup, args.ladderize,
        args.label_mode, args.label_internal, args.show_support, args.support_min,
        args.width, args.height, args.per_tip, args.font_size, args.output, args.dpi, args.no_show
    )

if __name__ == "__main__":
    main()
