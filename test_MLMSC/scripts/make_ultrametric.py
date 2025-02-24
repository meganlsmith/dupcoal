"""This script takes a non-ultrametric tree (due to rounding issues) and corrects it."""

import dendropy
import sys
import numpy as np
from scipy.stats import pearsonr

def main(input_tree, output_tree):
    # Load the tree
    tree = dendropy.Tree.get(path=input_tree, schema="newick")

    # Store original branch lengths
    original_lengths = [edge.length for edge in tree.preorder_edge_iter() if edge.length is not None]

    try:
        tree.calc_node_ages()
        print("No need for adjustments.")
    except:
        print("Failed to calculate node ages before adjustment.")

        # Find the maximum root-to-tip distance
        max_root_to_tip = max(n.distance_from_root() for n in tree.leaf_node_iter())

        # Adjust all terminal branch lengths
        for leaf in tree.leaf_node_iter():
            current_distance = leaf.distance_from_root()
            correction = max_root_to_tip - current_distance
            leaf.edge.length += correction  # Adjust branch length to make it ultrametric

        try:
            tree.calc_node_ages()
        except:
            print("Failed to calculate node ages after adjustment.")
            sys.exit(1)

        # Store new branch lengths
        adjusted_lengths = [edge.length for edge in tree.preorder_edge_iter() if edge.length is not None]

        # Compute correlation if branch lengths exist
        correlation, _ = pearsonr(original_lengths, adjusted_lengths)
        print(f"Pearson correlation between original and adjusted branch lengths: {correlation:.6f}")

        # Save the corrected tree
        tree.write(path=f"{output_tree}.newick", schema="newick")
        tree.write(path=f"{output_tree}.nexus", schema="nexus")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fix_ultrametric.py <input_tree> <output_tree>")
        sys.exit(1)

    input_tree_file = sys.argv[1]
    output_tree_file = sys.argv[2]
    main(input_tree_file, output_tree_file)
