# dupcoal
Simulate gene trees under a model that allows copy number hemiplasy.

For an extended description of the simulator, please see 'Extended_Description_v1a.pdf'.

# Arguments
* stree
    + type = str
    + help = Path to file with newick formatted species tree to use for simulation. Branch lengths should be in coalescent units.

* mu_par
    + type = float
    + help = loss rate.

* lambda_par
    + type = float
    + help = duplication rate.

* reps
    + type = int
    + help = Number of gene families to simulate.

* output
    + type = str
    + help = Folder for storing results.

# Outputs

* summary.csv = Table containing the following: replicate,duplications,losses,CNH,RKH,ILS,ILS_DLCPAR (see definitions below)
* trees.tre = All generated gene family trees.

# Table file explanation

* duplications: number of duplications
* observed_duplications: number of observable duplications (i.e., those not entirely obscured by subsequent losses)
* losses: number of losses
* Copy Number Hemiplasy (CNH): number of cases in which a mutation is placed on a branch that does not exist in the species tree.
* Rasmussen and Kellis CNH: number of cases in which a mutation is placed on a branch that does exist in the species tree, but does not match the current branch of the species tree. This is the phenomenon described as hemiplasy in Rasmussen and Kellis (2012, Figure 2B)
* All ILS (ILS): This includes the number of discordant branches in the parent and daughter trees and the number of NNIs that should be needed to reconcile discordance that occurs when coalescing subtrees.
* All ILS DLCPar (ILS_DLCPAR): This includes the number of parent or daughter trees that are discordant and the number of times that a subtree joined a branch that did not match it in terms of taxon composition.

# Example Usage
python3 dupcoal_v3.0.0.py --stree sp_tree.tre --mu_par 0 --lambda_par 0.3 --reps 10 --output example


