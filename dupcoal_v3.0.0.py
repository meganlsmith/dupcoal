import dendropy
from dendropy.simulate import treesim
import numpy as np
import random
import sys
import copy
import io
import argparse
import os
from collections import deque


class SubtreeNode:
    def __init__(self, tree, name=None):
        """
        A node to represent a subtree with metadata.
        
        Args:
        - tree: The actual tree (e.g., a DendroPy tree or a string placeholder).
        - name: Optional string for naming the subtree (for identification).
        """
        self.tree = tree
        self.name = name
        self.children = []    # List of child SubtreeNodes
        self.parent = None    # Reference to the parent SubtreeNode

    def add_child(self, child_node):
        """Add a child node and set its parent."""
        child_node.parent = self
        self.children.append(child_node)

    def find_node_by_name(self, node, name):

        if node.name == name:
            return node


        for child in node.children:
            found_node = self.find_node_by_name(child, name)  # Recursive call
            if found_node:
                return found_node

    def find_node_by_tree(self, node, tree):

        if node.tree == tree:
            return node


        for child in node.children:
            found_node = self.find_node_by_tree(child, tree)  # Recursive call
            if found_node:
                return found_node

    def level_order_traversal_bottom_up(self, root, sp_tree):

        # set a counter for naming copies
        copyid = 1

        # set other ocunters
        ils_joining_dlcpar_count = 0
        nni_joining_count = 0

        # Create a queue for level-order traversal
        queue = deque([(root, 0)])  # (node, depth)
        levels = {}  # To store nodes by level (depth)

        # Standard level-order traversal, recording levels
        while queue:
            node, depth = queue.popleft()

            if depth not in levels:
                levels[depth] = []

            levels[depth].append(node)

            # Add children to the queue
            for child in node.children:
                queue.append((child, depth + 1))

        # set the root edge in the species tree to infinity
        for edge in sp_tree.preorder_edge_iter():
            edge.length = np.inf
            break

        # Process the levels in reverse order (tips to root)
        for depth in sorted(levels.keys(), reverse=True):
            
            levelheights = []
            for theitem in levels[depth]:
                levelheights.append(get_tree_height(theitem.tree))

            sorted_pairs = sorted(zip(levelheights, levels[depth]), key=lambda x: x[0], reverse=True)
            sorted_ages, sorted_levels = zip(*sorted_pairs)
            levels[depth] = list(sorted_levels)

            for node in levels[depth]:
                #print(f"Depth {depth}: Node: {node.name}, Tree: {node.tree}")

                # Do what we need to do (i.e., coalesce to parent)
                if not depth==0:
                    #print(f"Coalesce this daughter tree {node.tree}.")
                    #print(f"Coalesce it to its parent: {node.parent.tree}")

                    # get the age of the present copy
                    age = get_tree_height(node.tree)
                    #original_age = age.copy()

                    # get the age of its parent
                    parent_height = get_tree_height(node.parent.tree)
                    if node.parent.name == "Root":
                        parent_height = np.inf

                    # get the list of leaves in the subtree.
                    subtree_leaves = [str(x.taxon).split()[0].strip("'") for x in node.tree.leaf_nodes()]

                    # Boolean to track whether we have managed to coalesce.
                    coalesced = False

                    # Use a while loop to find a place to coalesce the subtree.
                    while coalesced == False and age < parent_height:

                        # iterate over the edges, and check the timing of events against the age of the duplication and the leaves that should be present.
                        for edge in sp_tree.preorder_edge_iter():
                            min_age = edge.head_node.distance_from_tip()
                            try:
                                max_age = edge.tail_node.distance_from_tip()
                            except:
                                max_age = np.inf
                            if age >= min_age and age < max_age:
                                leaves = get_leaves(edge)
                                if len(intersection(leaves, subtree_leaves))> 0:
                                    # record the branch that is present at the correct time and has the correct leaves as the branch_to_coalesce.
                                    branch_to_coalesce = [edge, leaves, min_age, max_age, max_age - age]
                                    original_branch_leaves = leaves
                                    #print(f"Can coalsece in species tree branch {leaves} with ages ({min_age},{max_age})")

                        # Boolean to note whether we still need to look for the coalescence edge.
                        find_edge_to_coal = True
                        # while loop to find the edge of the parent gene tree that we should coalesce to.
                        while find_edge_to_coal == True:
                            # list to store edge heights (so we can know when first coalescent event happens)
                            edge_heights = []
                            possible_edges_to_coalesce = []
                            # iterate over the edges, and check the timing of events against the age of the duplication and the leaves that should be present.
                            for edge in node.parent.tree.levelorder_edge_iter():
                                min_age = edge.head_node.distance_from_tip()
                                try:
                                    max_age = edge.tail_node.distance_from_tip()
                                except:
                                    max_age = parent_height


                                if age >= min_age and age < max_age:
                                    leaves = get_leaves(edge)
                                    intersect_leaves = [x.split()[0].strip("'") for x in leaves]
                                    if len(intersection(intersect_leaves, branch_to_coalesce[1]))> 0:
                                        possible_edges_to_coalesce.append([edge, leaves, min_age, max_age])
                                        edge_heights.append(max_age)
                                        #print(f"Can coalesce to parent tree branch {leaves} with ages ({min_age, max_age})")


                            # draw a coalescent time
                            combinations = len(possible_edges_to_coalesce)
                            tcoal = np.random.exponential(scale = 1/combinations)
                            my_edge_min= min(edge_heights)
                            #print(f"Drew a coalescent time of {tcoal}, for a height of {tcoal+age}, which we will now compare to {branch_to_coalesce[3]} and {parent_height} and {my_edge_min}")

                            # do we coalesce before the shortest edge ends?
                            if tcoal+age > my_edge_min and find_edge_to_coal == True:
                                if tcoal+age > branch_to_coalesce[3] or tcoal+age>parent_height:
                                    # if we have failed and are now in a different edge of the species tree or outside the parent tree, we can break the parent gene tree loop and go get the new branch of the species tree.
                                    #print(f"We could not coalesce here, because {tcoal+age} exceeded the age of either the species tree branch {branch_to_coalesce[3]} or the parent height {parent_height}.")
                                    find_edge_to_coal = False
                                else:
                                    # otherwise, we failed to coalesce but still have time set the new age to the minimum edge height.
                                    #print(f"We didn't coalesce before the first event, but we still have time.")
                                    age = my_edge_min
                            else:
                                # if we haven't failed, then we know which edge to coalesce to.
                                #print(f"Successfully chose an edge for coalescence at height {tcoal+age}.")
                                find_edge_to_coal = False



                        # did coalescence happen in this branch?
                        if tcoal + age < branch_to_coalesce[3] and tcoal+age < parent_height:

                            #print(f"We were successful, and will coalesce the trees.")
                            copyid+=1

                            # select an edge from available at random
                            the_coalesced_edge_info = random.sample(possible_edges_to_coalesce, k=1)[0]
                            the_coalesced_edge = the_coalesced_edge_info[0]

                            # get relevant nodes
                            original_child = the_coalesced_edge.head_node
                            original_edge_length = the_coalesced_edge.length
                            original_parent_node = the_coalesced_edge.tail_node

                            # check for any mismatch
                            parent_leaves_check = [x.split()[0].strip("'") for x in get_leaves_node(original_child)]
                            daughter_leaves_check = [x.split()[0].strip("'") for x in get_leaves_node(node.tree.seed_node)]
                            if set(parent_leaves_check) != set(daughter_leaves_check):
                                ils_joining_dlcpar_count += 1
                                #print(f'The parent tree leaves are {get_leaves_node(original_child)}')
                                #print(f'The daughter tree leaves are {get_leaves_node(node.tree.seed_node)}')
                            
                                nni_joining_count_current = self.check_nni_joining(the_coalesced_edge, get_leaves_node(node.tree.seed_node), sp_tree)
                                nni_joining_count += nni_joining_count_current

                            # update taxon namespace of new tree
                            for leaf in node.tree.leaf_node_iter():
                                if str(leaf.taxon.label).split()[1].strip("'") == "1":
                                    new_taxon = str(leaf.taxon.label).split()[0].strip("'") + ' '+ str(copyid)
                                    leaf.taxon.label = new_taxon
                                    leaf.taxon = dendropy.Taxon(new_taxon)

                            # adjust branch length to account for any failures to coalesce to the initial parent branches
                            for edge in node.tree.preorder_edge_iter():
                                if edge.length is None:
                                    edge.length = 0
                                new_length = age - (edge.head_node.distance_from_tip()+edge.length) + edge.length
                                edge.length = new_length
                                break

                            # adjust branch length in subtree to include coalescence time
                            for edge in node.tree.preorder_edge_iter():
                                edge.length = edge.length + tcoal
                                total_edge = edge.length
                                total_height = edge.head_node.distance_from_tip() + edge.length
                                break

                            #print(f"The tree to coalesce is {node.tree} with height {total_height}")

                            # adjust branch length in parent tree
                            for edge in node.parent.tree.preorder_edge_iter():
                                 if edge == the_coalesced_edge:
                                    if edge.length == None:
                                        edge.length = total_height - edge.head_node.distance_from_tip()
                                        
                                    elif edge.length < total_height - edge.head_node.distance_from_tip():
                                        edge.length = total_height - edge.head_node.distance_from_tip()
                                        previous_edge = edge.length   
                                    else:
                                        previous_edge = edge.length
                                        edge.length = total_height - edge.head_node.distance_from_tip()

                                    # branch lengths for daughter nodes

                                    this_new_length = edge.length

                            #print(f"The tree to coalesce to is {node.parent.tree} at branch {the_coalesced_edge_info}")

                            # create new subtree
                            new_subtree = dendropy.Tree()
                            for item in node.tree.taxon_namespace:
                                if item not in new_subtree.taxon_namespace:
                                    new_subtree.taxon_namespace.add_taxon(item)
                            for item in node.parent.tree.taxon_namespace:
                                if item not in new_subtree.taxon_namespace:
                                    new_subtree.taxon_namespace.add_taxon(item)
                            cnode = new_subtree.seed_node.new_child(edge_length = previous_edge - this_new_length)
                            cnode.add_child(node.tree.seed_node)
                            cnode.add_child(original_child)
                            new_subtree.reroot_at_node(new_subtree.seed_node)

                            #print(f"We have created the new subtree: {new_subtree}.")

                            for taxon in new_subtree.taxon_namespace:
                                if not taxon in node.parent.tree.taxon_namespace:
                                    node.parent.tree.taxon_namespace.add_taxon(taxon)
                            # add new subtree to old subtree
                            if the_coalesced_edge.head_node == node.parent.tree.seed_node:
                                node.parent.tree = new_subtree
                            else:
                                original_parent_node.remove_child(original_child)
                                newch = original_parent_node.new_child(edge_length = previous_edge - this_new_length)
                                newch.add_child(node.tree.seed_node)
                                newch.add_child(original_child)

                            #node.parent.tree.calc_node_ages()
                            coalesced = True
                            node.parent.tree.calc_node_ages()

                            # reformat tree
                            output_stream = io.StringIO()
                            sys.stdout = output_stream

                            # convert the tree to a string and print it
                            print(str(node.parent.tree))

                            # redirect stdout back to its original destination
                            sys.stdout = sys.__stdout__

                            # read the captured output from the IOStream object
                            output = output_stream.getvalue()
                            output = output + ';'

                            # print the captured output
                            node.parent.tree = dendropy.Tree.get(data=output, schema="newick")
                            node.parent.tree.calc_node_ages()

                        elif tcoal+age > branch_to_coalesce[3] and tcoal+age <= parent_height:
                            # we failed to coalesce in this branch of the species tree
                            #print(f"We will move to a new branch of the species tree.")
                            age = tcoal+age
                            del branch_to_coalesce

                        else:
                            # we failed to coalesce before the origin of the parent tree
                            age = tcoal+age

                            # find parent's parent
                            parent_info = self.find_node_by_tree(root, node.parent.tree)
                            node.parent = parent_info.parent
                            
                            # this is where we should update thge height inthe subtree and add it to the dictionary at the right depth
                            for edge in node.tree.preorder_edge_iter():
                                edge.length = parent_height
                                break

                            levels[depth-1].append(node)


        return(ils_joining_dlcpar_count, nni_joining_count)

    def check_nni_joining(self, the_coalesced_edge, subtreeleaves, sp_tree):

        ref_tree = copy.deepcopy(sp_tree)
        nni = 0

        # get leaf info
        joining_leaves = [x.split()[0] for x in get_leaves(the_coalesced_edge)]
        subtreeleaves = [x.split()[0] for x in subtreeleaves]
        all_leaves_to_join = set(joining_leaves + subtreeleaves)

        # get mrca
        count=0
        for thenode in ref_tree.levelorder_node_iter():

            node_leaves = get_leaves_node(thenode)
            #node_leaves = [x.strip("'").split()[0] for x in node_leaves]

            if count == 0:
                mrca_leaves = node_leaves.copy()
                count+=1
            elif set(all_leaves_to_join).issubset(set(node_leaves)):
                mrca_leaves = node_leaves.copy()

        for thenode in ref_tree.levelorder_node_iter():
            node_leaves = get_leaves_node(thenode)

            # is it the mrca
            if check_lists_equality(node_leaves, mrca_leaves):
                nni+=1
            elif set(mrca_leaves).issubset(set(node_leaves)):
                continue
            elif set(subtreeleaves).issubset(set(node_leaves)):
                if not check_lists_equality(subtreeleaves, node_leaves):
                    nni+=1

        for thenode in ref_tree.levelorder_node_iter():
            node_leaves = get_leaves_node(thenode)

            # is it the mrca
            if check_lists_equality(node_leaves, mrca_leaves):
                continue
            elif set(mrca_leaves).issubset(set(node_leaves)):
                continue
            elif set(joining_leaves).issubset(set(node_leaves)):
                if not check_lists_equality(joining_leaves, node_leaves):
                    nni+=1

        return nni

class mlmsc:

    def __init__(self, sp_tree, lambda_par, mu_par, og_spheight):
        self.sp_tree = sp_tree
        self.lambda_par = lambda_par
        self.mu_par = mu_par
        self.og_spheight = og_spheight

    def generate(self):

        # counter to track duplications
        duplication_count = 0
        cnh_count = 0
        rkh_count = 0
        ils_count = 0
        ils_dlcpar_count = 0

        # sample the parent tree from the MSC
        parent_tree = self._get_gene_tree()
        original_parent = copy.deepcopy(parent_tree)
        current_ils_count, current_ils_dlcpar_count = self._tree_differences(parent_tree)
        ils_count += current_ils_count
        ils_dlcpar_count += current_ils_dlcpar_count

        # store in tree
        root = SubtreeNode(tree=parent_tree, name="Root")

        # return this if there is no duplication
        if self.lambda_par == 0:
            return(root, duplication_count, cnh_count, rkh_count, ils_count, ils_dlcpar_count, [], [], original_parent, [])

        # place duplications on tree
        original_subtrees, mutated_subtrees, current_duplication_count, current_cnh_count, current_rkh_count, current_ils_count, current_ils_dlcpar_count, ages = self._birth_mlmsc(parent_tree)
        duplication_count+= current_duplication_count
        cnh_count+=current_cnh_count
        rkh_count += current_rkh_count
        ils_count += current_ils_count
        ils_dlcpar_count += current_ils_dlcpar_count


        # process mutated subtrees iteratively
        to_process = mutated_subtrees.copy()  # List of subtrees to process
        all_original_subtrees = original_subtrees.copy()  # Keep track of all original subtrees
        all_mutated_subtrees = mutated_subtrees.copy()
        all_ages = ages.copy()

        # add to tree
        for subtree in mutated_subtrees:
            subtree_str = SubtreeNode(subtree, name=subtree.as_string(schema="newick").strip())
            root.add_child(subtree_str)

        while to_process:
            # get the next subtree to process
            current_subtree = to_process.pop(0)

            # generate duplications for the current subtree
            new_original, new_mutated, current_duplication_count, current_cnh_count, current_rkh_count, current_ils_count, current_ils_dlcpar_count, current_ages = self._birth_mlmsc(current_subtree)
            duplication_count+=current_duplication_count
            cnh_count+=current_cnh_count
            rkh_count += current_rkh_count
            ils_count += current_ils_count
            ils_dlcpar_count += current_ils_dlcpar_count

            # add new original subtrees to the overall list
            all_original_subtrees.extend(new_original)
            all_mutated_subtrees.extend(new_mutated)
            all_ages.extend(current_ages)

            # add to tree
            for subtree in new_mutated:
                subtree_str = SubtreeNode(subtree, name=subtree.as_string(schema="newick").strip())
                get_parent = root.find_node_by_name(root, current_subtree.as_string(schema="newick").strip())
                get_parent.add_child(subtree_str)
                to_process.append(subtree)

        return(root, duplication_count, cnh_count, rkh_count, ils_count, ils_dlcpar_count, all_original_subtrees, all_mutated_subtrees, original_parent, all_ages)

    def coalesce(self, root):
        """Coalesce subtrees."""

        # set the root edge in the species tree to infinity
        for edge in self.sp_tree.preorder_edge_iter():
            edge.length = np.inf
            break

        # iterate over the loci. Begin at the tips (with the last round of duplication), and process all tips before moving to the next level. End at the root.
        ils_joining_dlcpar_count, nni_joining_count = root.level_order_traversal_bottom_up(root, self.sp_tree)

        #print(f"We have coalesced everything into: {root.tree}.")
        return (root.tree, ils_joining_dlcpar_count, nni_joining_count)

    def losses(self, combined_tree):
        """Place losses."""
        loss_count = 0
        cnh_count = 0
        rkh_count = 0

        # return tree if loss is zero
        if self.mu_par == 0:
            return(combined_tree, 0, 0, 0)

        edges_to_lose = []
        leaves_to_remove = []

        for edge in combined_tree.preorder_edge_iter():

            if edge.length == 0 or edge.length == None: # skip the root edge
                continue
            
            #if (edge.length + edge.head_node.age) > self.og_spheight:
            #    print('do something')

            # draw the time of the next event from an exponential distribution.
            current_scale_parameter =  1/(self.mu_par)
            ti = np.random.exponential(scale = current_scale_parameter)

            # check whether event occurs on edge
            if ti < edge.length:
                    
                event_age = (edge.length - (ti)) + edge.head_node.age

                if set(get_leaves(edge)).isdisjoint(leaves_to_remove):
                    edges_to_lose.append(edge)
                    for item in edge.head_node.leaf_iter():
                        leaves_to_remove.append(item.taxon.label)
                    loss_count += 1
                
                    # check for CNH
                    match=False
                    edge_leaves = [x.split()[0].strip("'") for x in get_leaves(edge)]
                    for spnode in self.sp_tree:
                        leaves_of_node = get_leaves_node(spnode)
                        if set(leaves_of_node) == set(edge_leaves):
                            match=True
                    cnh_count+= (1-match)

                    # check for RK CNH
                    # find sp branch
                    sp_leaves = self._get_leaves_from_species(edge_leaves, event_age)
                    if set(sp_leaves) != set(edge_leaves):
                        rkmatch=True
                    else:
                        rkmatch=False
                    rkh_count += rkmatch



        # remove any daughters of the selected edge
        taxa_to_remove = []
        for edge in edges_to_lose:
            for item in edge.head_node.leaf_iter():
                taxa_to_remove.append(item.taxon)
        taxa_to_remove = set(taxa_to_remove)

        if len(set(taxa_to_remove)) == len(combined_tree.taxon_namespace):
            return (None, loss_count, cnh_count, rkh_count)

        else:
            updated_tree = combined_tree.extract_tree_without_taxa(taxa_to_remove)
            return(updated_tree, loss_count, cnh_count, rkh_count)

    def _get_gene_tree(self):
        gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
            containing_taxon_namespace=self.sp_tree.taxon_namespace,
            num_contained=1)
        gene_tree = treesim.contained_coalescent_tree(containing_tree=self.sp_tree, gene_to_containing_taxon_map=gene_to_species_map)
        gene_tree.calc_node_root_distances()
        gene_tree.calc_node_ages()

        return(gene_tree)

    def _birth_mlmsc(self, parent_tree):

        current_duplication_count = 0
        current_cnh_count = 0
        current_rkh_count = 0
        current_ils_count = 0
        current_ils_dlcpar_count = 0

        original_locus_trees = []
        mutated_locus_trees = []
        ages_locus_trees = []


        for edge in parent_tree.preorder_edge_iter():

            locus_leaves = get_leaves(edge)

            if edge.length == None: # skip root edge
                continue

            # set tc to 0
            tc = 0.0

            while tc < edge.length:

                # draw the time of the next event from an exponential distribution.
                current_scale_parameter =  1/(self.lambda_par)
                ti = np.random.exponential(scale = current_scale_parameter)

                tc+=ti

                # check whether event occurs on edge
                if tc < edge.length:



                    # calculate the event age (time from the present)
                    event_age = (edge.length - (tc)) + edge.head_node.age

                    # count the event
                    current_duplication_count+=1

                    ages_locus_trees.append(event_age)

                    # get the possible leaves from the species tree
                    sp_leaves = self._get_leaves_from_species(locus_leaves, event_age)

                    # draw a daugther gene tree
                    new_locus_tree = self._get_gene_tree()
                    original_locus_trees.append(new_locus_tree)

                    # add the duplication
                    mutated_subtree, cnh_indicator, rkh_indicator = self._add_duplication(event_age, sp_leaves, new_locus_tree)
                    current_cnh_count+= cnh_indicator
                    current_rkh_count += rkh_indicator

                    # count ils
                    temp_ils_count, temp_ils_dlcpar_count = self._count_ils(mutated_subtree)
                    current_ils_count += temp_ils_count
                    current_ils_dlcpar_count += temp_ils_dlcpar_count

                    mutated_locus_trees.append(mutated_subtree)


        if len(mutated_locus_trees) > 0:
            sorted_pairs = sorted(zip(ages_locus_trees, mutated_locus_trees), key=lambda x: x[0], reverse=True)
            sorted_ages, sorted_mutated_locus_trees = zip(*sorted_pairs)
            mutated_locus_trees = list(sorted_mutated_locus_trees)

        return(original_locus_trees, mutated_locus_trees, current_duplication_count, current_cnh_count, current_rkh_count, current_ils_count, current_ils_dlcpar_count, ages_locus_trees)

    def _get_leaves_from_species(self, locus_leaves, event_age):

        locus_leaves_comp = [x.split()[0] for x in locus_leaves]
        relevant_leaves = []

        # iterate over the edges of the gene tree. If an edge is a) alive at the right time and b) contains a relevant descedent, add it to the options
        for edge in self.sp_tree.preorder_edge_iter():
            if edge.tail_node != None:
                relevant_time_frame = [edge.head_node.age, edge.head_node.age + edge.length]
                if event_age > relevant_time_frame[0] and event_age < relevant_time_frame[1] and set(locus_leaves_comp).issubset(set(get_leaves(edge))):
                    relevant_leaves = get_leaves(edge)
            else:
                all_leaves = get_leaves(edge)
        
        if len(relevant_leaves) == 0:
            relevant_leaves = all_leaves
        
        return(relevant_leaves)

    def _add_duplication(self, event_age, sp_leaves, new_subtree):

        # iterate over the edges of the gene tree. If an edge is a) alive at the right time and b) contains a relevant descedent, add it to the options
        potential_edges = []
        for gtedge in new_subtree.preorder_edge_iter():
            if gtedge.tail_node != None:
                relevant_time_frame = [gtedge.head_node.age, gtedge.head_node.age + gtedge.length]
                if event_age > relevant_time_frame[0] and event_age < relevant_time_frame[1]:
                    gt_leaves = [x.split()[0] for x in get_leaves(gtedge)]
                    if set(gt_leaves).issubset(sp_leaves):
                        potential_edges.append(gtedge)

        if len(potential_edges) > 0:

            # sample an edge
            the_edge = random.sample(potential_edges, 1)

            # extract the subtree
            leaves_to_sample = get_leaves(the_edge[0])
            subtree = new_subtree.extract_tree_with_taxa_labels(leaves_to_sample)
            subtree.calc_node_ages()

            for subtree_node in subtree.preorder_node_iter():
                subtending_length = event_age - subtree_node.age
                break

            # truncate the branch to the duplication time'
            for subtree_edge in subtree.preorder_edge_iter():
                subtree_edge.length = subtending_length
                break

        else:
            subtree = new_subtree
            for thenode in subtree.preorder_node_iter():
                treeheight = thenode.age
                break
            for edge in subtree.preorder_edge_iter():
                edge.length = event_age - treeheight
                break

        subtree.annotations['age'] = float(event_age)

        # check whether the root branch of the subtree is in the species tree
        leaves_of_root = get_leaves_node(subtree.seed_node)
        leaves_of_root = [x.split()[0].strip("'") for x in leaves_of_root]

        # check for hemiplasy
        match=False
        for spnode in self.sp_tree:
            leaves_of_node = get_leaves_node(spnode)
            if set(leaves_of_node) == set(leaves_of_root):
                match=True
        
        # check for RK hemiplasy
        if set(leaves_of_root) != set(sp_leaves):
            rkmatch =True
        else:
            rkmatch = False
        return(subtree, 1-match, rkmatch)

    def _count_ils(self, gene_tree):
        # prune species tree
        leaves_to_keep = []
        for edge in gene_tree.postorder_edge_iter():
            leaves_to_keep.append(get_leaves(edge))
        leaves_to_keep = [item for row in leaves_to_keep for item in row]
        leaves_to_keep = set(leaves_to_keep)
        leaves_to_keep = list(leaves_to_keep)
        leaves_to_keep = [x.split()[0] for x in leaves_to_keep]

        species_subtree = self.sp_tree.extract_tree_with_taxa_labels(leaves_to_keep)
        species_subtree.calc_node_ages()

        ils_count, ils_dlcpar_count = self._tree_differences(gene_tree=gene_tree)

        return(ils_count, ils_dlcpar_count)           

    def _tree_differences(self, gene_tree):
        ils_count = 0
        ils_dlcpar = 0
        for gtedge in gene_tree.postorder_edge_iter():
            match = False
            gt_relevant_leaves = get_leaves(gtedge)
            gt_relevant_leaves = [x.split()[0] for x in gt_relevant_leaves]
            for edge in self.sp_tree.postorder_edge_iter():
                relevant_leaves = get_leaves(edge)
                if check_lists_equality(relevant_leaves, gt_relevant_leaves) == True:
                    match = True
            if match == False:
                ils_count += 1
        if ils_count > 0:
            ils_dlcpar = 1

        return(ils_count, ils_dlcpar)


# utility functions using

def get_leaves_node(node):
    leaves = []
    for leaf in node.leaf_iter():
        leaves.append(str(leaf.taxon).strip("'"))
    return(leaves)

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def get_leaves(edge):
    leaves = []
    for leaf in edge.head_node.leaf_iter():
        leaves.append(str(leaf.taxon).strip("'"))
    return(leaves)

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate gene family trees.")
    parser.add_argument('--stree', type=str, help="Path to file with newick formatted species tree to use for simulation. Branch lengths should be in coalescent units.")
    parser.add_argument('--mu_par', type=float, help="Loss rate.")
    parser.add_argument('--lambda_par', type=float, help="Duplication rate.")
    parser.add_argument('--reps', type=int, help="Number of gene families to simulate.")
    parser.add_argument('--output', type=str, help="Folder for storing results.")
    parser.add_argument('--verbose', type=int, help="Whether to output additional information. 1=output original and mutated trees.", default=0)
    args= parser.parse_args()
    return(args)

def write_trees(tree, treefile):
    if tree is None:
        treefile.write("None\n")
    else:
        treefile.write(tree.as_string(schema="newick", suppress_rooting=True))

def write_log(rep,total_dups, total_losses,  hemiplasy, rk_hemiplasy, all_nni, all_ils,logfile):
    logfile.write(f"{rep},{total_dups},{total_losses},{hemiplasy},{rk_hemiplasy},{all_nni},{all_ils}\n")


def check_lists_equality(list1, list2):
    return set(list1) == set(list2)

def get_tree_height(tree):
    for thenode in tree.preorder_node_iter():
        height = thenode.age
        break
    for edge in tree.preorder_edge_iter():
        if edge.length is None:
            break
        height+=edge.length
        break
    return(height)

def write_subtrees(original_subtrees, mutated_subtrees, output, rep, parent_tree, ages):
    otrees = os.path.join(output, f"original_trees_{rep}.trees")
    mtrees = os.path.join(output, f"mutated_trees_{rep}.trees")
    agesfile = os.path.join(output, f"ages_{rep}.txt")
    with open(otrees, 'w') as f:
        f.write(parent_tree.as_string(schema="newick"))
        for item in original_subtrees:
            f.write(item.as_string(schema="newick"))
    with open(mtrees, 'w') as f:
        f.write(parent_tree.as_string(schema="newick"))
        for item in mutated_subtrees:
            f.write(item.as_string(schema="newick"))
    with open(agesfile, 'w') as f:
        for item in ages:
            f.write(f"{item}\n")


def main():
             
    # get arguments
    args = parse_args()

    # create output folder
    if os.path.exists(args.output):
        sys.exit('The output folder already exists. Please remove, or change the output folder name and try again.')
    else:
        os.system('mkdir %s' % args.output)
        
    # specify species tree
    sp_tree_main = dendropy.Tree.get(path=args.stree, schema="newick")
    sp_tree_main.calc_node_root_distances()
    sp_tree_main.calc_node_ages()

    # set the arugments
    lambda_par = args.lambda_par
    mu_par = args.mu_par

    # create tree file and tables
    treefile = open(os.path.join(args.output, 'trees.tre'), 'w')
    logfile = open(os.path.join(args.output, 'summary.csv'), 'w')
    logfile.write('replicate,duplications,losses,CNH,RKH,ILS,ILS_DLCPAR\n')

    for i in range(args.reps):

        # copy of species tree
        sp_tree = sp_tree_main.clone()

        # perform top-down birth death under the modified model
        mlmsc_simulator = mlmsc(sp_tree, lambda_par, mu_par, get_tree_height(sp_tree))
        simulated_trees, duplication_count, cnh_dupcount, rkh_dupcount, ils_count, ils_dlcpar_count, original_subtrees, mutated_subtrees, parent_tree, ages = mlmsc_simulator.generate()
        combined_tree, ils_joining_dlcpar_count, nni_joining_count = mlmsc_simulator.coalesce(copy.deepcopy(simulated_trees))
        final_trees, loss_count, cnh_losscount, rkh_losscount = mlmsc_simulator.losses(combined_tree)

        cnh_count = cnh_dupcount + cnh_losscount
        rkh_count = rkh_dupcount + rkh_losscount
        all_ils = nni_joining_count + ils_count
        all_ils_dlcpar = ils_joining_dlcpar_count + ils_dlcpar_count

        print('replicate: %s' % i)

        write_trees(final_trees, treefile)

        write_log(str(i),duplication_count, loss_count,  cnh_count, rkh_count, all_ils, all_ils_dlcpar, logfile)

        if args.verbose == 1:
            write_subtrees(original_subtrees, mutated_subtrees, args.output, str(i), parent_tree, ages)

        del final_trees 
        del simulated_trees
        del duplication_count
        del cnh_dupcount
        del rkh_dupcount
        del ils_count
        del ils_dlcpar_count
        del combined_tree
        del ils_joining_dlcpar_count
        del nni_joining_count
        del loss_count
        del cnh_losscount
        del rkh_losscount
        del original_subtrees
        del mutated_subtrees

    treefile.close()
    logfile.close()


if __name__ == "__main__":
    main()
