import dendropy

# Original trees in Newick format
tree1_newick = "((B_3:1.7781073387695634,A_1:1.7781073387695634):0.06909242386938486);"
tree2_newick = "((C_1:0.6307895249054567,B_1:0.6307895249054567):1.2164102377334916,A_1:1.7781073387695634):0.07364269584997318;"

# A1_B1 subtree in Newick format
subtree_newick = "(A_1:0.5,B_1:0.5);"

# Load the trees
tree1 = dendropy.Tree.get(data=tree1_newick, schema="newick")
tree2 = dendropy.Tree.get(data=tree2_newick, schema="newick")
subtree = dendropy.Tree.get(data=subtree_newick, schema="newick")

# Function to replace node
def replace_node_with_subtree(tree, node_label, subtree):
    # Find the node to replace
    node_to_replace = tree.find_node_with_taxon_label(node_label)
    if not node_to_replace:
        print(f"Node {node_label} not found in the tree.")
        return

    # Re-root the subtree to ensure it's properly attached
    new_subtree = subtree.clone(depth=2)
    
    # Attach the new subtree to the parent node
    parent_node = node_to_replace.parent_node
    edge_length = node_to_replace.edge_length
    parent_node.remove_child(node_to_replace)

    for child in new_subtree.seed_node.child_nodes():
        child.edge_length += edge_length  # Add the original edge length
        parent_node.add_child(child)

# Replace A_1 with the A1_B1 subtree in both trees
replace_node_with_subtree(tree1, "A_1", subtree)
replace_node_with_subtree(tree2, "A_1", subtree)

# Print the updated trees
print("Updated Tree 1:", tree1.as_string(schema="newick"))
print("Updated Tree 2:", tree2.as_string(schema="newick"))
