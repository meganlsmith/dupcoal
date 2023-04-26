# Copy Number Hemiplasy Simulator
Simulate gene trees under a model that allows copy number hemiplasy.

# Arguments
* stree
    + type=str
    + help=Path to file with newick formatted species tree to use for simulation. Branch lengths should be in coalescent units.

* mu_par
    + type=float
    + help="Duplication rate."

* lambda_par
    + type=float
    + help="Loss rate."

* reps
    + type=int
    + help="Number of gene families to simulate."

* output
    + type=str
    + help="Folder for storing results."

* root_length
    + type=float
    + help="Arbitrarily long root for the species tree (should be large enough that all coalescences should occur.)"


# Example Usage
python3 simulator_v1e.py --stree sp_tree.tre --mu_par 0.3 --lambda_par 0.5 --reps 100 --output example_trees --root_length 100    
