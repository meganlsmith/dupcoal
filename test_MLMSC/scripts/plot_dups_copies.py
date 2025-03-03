"""This script will plot the number of duplications and copies across programs, Ne, ld, and lb values."""

import argparse
import os
import glob
import sqlite3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from ete3 import Tree
import sys
import numpy as np

def create_map(tree, temp_output, sptree):
    thenames = tree.get_leaf_names()
    spdict = {}
    for item in thenames:
        species = item.split('_')[0]
        if species in spdict:
            if not item in spdict[species]:
                spdict[species].append(item)
        else:
            spdict[species] = [item]
    
    with open(temp_output, 'w') as f:
        for key in spdict.keys():
            string = f"{key}:"
            for item in spdict[key]:
                string += f"{item},"
            string.strip(",")
            f.write(string + "\n")

    
    sptree = Tree(sptree)
    to_prune = [x for x in sptree.get_leaf_names() if x in spdict.keys()]
    sptree.prune(to_prune)
    with open('sp.temp.tre', 'w') as f:
        f.write(sptree.write())

def astral(line, args, program):

    totalquartets=np.nan
    quartetscore=np.nan
    normquartetscore=np.nan

    # rename tree tips if mlmsc
    tree = Tree(line.strip())
    if program == "mlmsc":
        tree = rename_tree(tree)
        tree.write(outfile="temp.tre")

    else:
        with open("temp.tre", 'w') as f:
            f.write(line.strip())

    # create map
    create_map(tree, 'temp_map.txt', args.species)


    # astral score
    command = f"java -jar {args.astral} -i temp.tre -o temp_out.tre -q sp.temp.tre -a temp_map.txt 2> temp.log"
    os.system(command)

    # get score
    with open('temp.log', 'r') as logf:
        for line in logf.readlines():
            if line.startswith("Number of quartet trees in the gene trees: "):
                totalquartets=line.split(": ")[1].strip()
            elif line.startswith("Final quartet score is: "):
                quartetscore=line.split(": ")[1].strip()
            elif line.startswith("Final normalized quartet score is: "):
                normquartetscore=line.split(": ")[1].strip()
    os.system("rm temp_map.txt temp.tre temp.log sp.temp.tre temp_out.tre")
    
    return(totalquartets, quartetscore, normquartetscore)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--simphy", required=False, help="Path to SimPhy results folder")
    parser.add_argument("--dupcoal", required=False, help="Path to DupCoal results folder")
    parser.add_argument("--mlmsc", required=False, help="Path to MLMSCII results folder")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output graph files")
    parser.add_argument("--astral", required=True, help="Path to ASTRAL.")
    parser.add_argument("--species", required=True, help="Path to species tree.")
    parser.add_argument("--discordance", action="store_true", help="Flag to calculate discordance.")
    return parser.parse_args()

def extract_simphy_data(folder, args):
    data = []
    for subfolder in os.listdir(folder):
        if not subfolder.startswith("simphy_"):
            continue
        params = subfolder.split("_")[1:]
        db_path = os.path.join(folder, subfolder, f"simphy_{params[0]}_{params[1]}_{params[2]}_{params[3]}.db")
        tree_path = os.path.join(folder, subfolder, "g_trees.tre")
        
        if os.path.exists(db_path) and os.path.exists(tree_path):
            conn = sqlite3.connect(db_path)
            df = pd.read_sql_query("SELECT n_dup FROM Locus_Trees", conn)

            conn.close()
            
            trees = [Tree(line.strip()) for line in open(tree_path) if line.strip()]
            tip_counts = [len(tree.get_leaves()) for tree in trees]
            unique_copies = []
            for tree in trees:
                unique_copies.append(len(set([x.split('_')[1] for x in tree.get_leaf_names()]))-1)

            if args.discordance: 
                totalquartets = []
                quartetscore = []
                normalizedscore = []
                for line in open(tree_path):
                    curr_totalquartets, curr_quartetscore, curr_normalizedscore = astral(line, args, program="simphy")
                    totalquartets.append(curr_totalquartets)
                    quartetscore.append(curr_quartetscore)
                    normalizedscore.append(curr_normalizedscore)

            
                for n_dup, n_tips, obs_dup, totalq, qscore, normalqscore in zip(df["n_dup"], tip_counts,unique_copies, totalquartets, quartetscore, normalizedscore):
                    data.append(["simphy"] + params + [n_dup, n_tips, obs_dup, totalq, qscore, normalqscore])
            else:
                for n_dup, n_tips, obs_dup in zip(df["n_dup"], tip_counts, unique_copies):
                    data.append(["simphy"] + params + [n_dup, n_tips, obs_dup])
    return data

def extract_csv_data(folder, program, dup_col, tree_pattern, args):
    data = []
    for subfolder in os.listdir(folder):
        if not subfolder.startswith(f"{program}_"):
            continue
        params = subfolder.split("_")[1:]
        csv_path = os.listdir(os.path.join(folder, subfolder))
        try:
            csv_path = [x for x in csv_path if x.endswith('.txt')][0]
        except:
            csv_path = [x for x in csv_path if x.endswith('.csv')][0]
        tree_path = glob.glob(os.path.join(folder, subfolder, tree_pattern))[0]
        csv_path = os.path.join(folder, subfolder, csv_path)
        
        if os.path.exists(csv_path) and tree_path:
            df = pd.read_csv(csv_path)
            tip_counts = []
            unique_copies = []

            if args.discordance: 
                totalquartets = []
                quartetscore = []
                normalizedscore = []
            with open(tree_path, 'r') as f:

                for line in f:

                    if line != "None\n" and line != "None":
                        tree = Tree(line)
                        tip_counts.append(len(tree.get_leaves()))
                        if program=="dupcoal":
                            unique_copies.append(len(set([x.split('_')[1] for x in tree.get_leaf_names()]))-1)
                        else:
                            unique_copies.append(np.nan)
                        if args.discordance:
                            curr_totalquartets, curr_quartetscore, curr_normalizedscore = astral(line, args, program=program)
                            totalquartets.append(curr_totalquartets)
                            quartetscore.append(curr_quartetscore)
                            normalizedscore.append(curr_normalizedscore)

                    else:
                        tip_counts.append(0)
                        if program=="dupcoal":
                            unique_copies.append(0)
                        else:
                            unique_copies.append(np.nan)
                        if args.discordance:
                            totalquartets.append(np.nan)
                            quartetscore.append(np.nan)
                            normalizedscore.append(np.nan)

            if program == "dupcoal":
                if args.discordance:
                    for n_dup, n_tips, obs_dup, totalq, qscore, normalqscore in zip(df[dup_col], tip_counts, unique_copies, totalquartets, quartetscore, normalizedscore):
                        data.append([program] + params + [n_dup, n_tips, obs_dup, totalq, qscore, normalqscore])
                else:
                    for n_dup, n_tips, obs_dup in zip(df[dup_col], tip_counts, unique_copies):
                        data.append([program] + params + [n_dup, n_tips, obs_dup])
            else:
                if args.discordance:
                    for n_dup, n_tips, obs_dup, totalq, qscore, normalqscore in zip(unique_copies, tip_counts, df[dup_col], totalquartets, quartetscore, normalizedscore):
                        data.append([program] + params + [n_dup, n_tips, obs_dup, totalq, qscore, normalqscore])
                else:
                    for n_dup, n_tips, obs_dup in zip(unique_copies, tip_counts, df[dup_col]):
                        data.append([program] + params + [n_dup, n_tips, obs_dup])
    return data

def plot_summaries(df, args):

    df["branch_length"] = 1 / df["c"]
    fig, axes = plt.subplots(df["ld"].nunique(), 2, figsize=(15, 6), sharex=True)

    for j, ld_val in enumerate(df['ld'].unique()):
        subset = df[(df["ld"]==ld_val)].copy()
        for k, metric in enumerate(["obs_dup", "gene_copies"]):
            ax = axes[j,k]
            sns.lineplot(data=subset, x="lb", y=metric, hue="program", ax=ax, errorbar="se", err_style="bars", style="branch_length")
            ax.set_title(f"{metric.replace('_', ' ').title()} (ld={ld_val})")
   
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_summary.png")

def plot_astral(df, args):
    fig, axes = plt.subplots(df["ld"].nunique(), 2, figsize=(15, 6), sharex=True)

    for j, ld_val in enumerate(df['ld'].unique()):
        subset = df[(df["ld"]==ld_val) & df["c"]==1].copy()
        for k, metric in enumerate(["total_quartets", "normalized_quartet_score"]):
            if df["ld"].nunique()>1:
                ax = axes[j,k]
            else:
                ax = axes[k]
            sns.lineplot(data=subset, x="duplications", y=metric, hue="program", ax=ax, errorbar="se", err_style="bars")
            ax.set_title(f"{metric.replace('_', ' ').title()} (ld={ld_val})")
   
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_astral.png")

def rename_tree(tree):
    sp_dict = {}
    for node in tree.traverse():
        if  node.is_leaf():
            if node.name in sp_dict:
                last_integer = sp_dict[node.name][-1]
                sp_dict[node.name].append(last_integer+1)
                node.name = f"{node.name}_{last_integer+1}"
            else:
                sp_dict[node.name] = [0]
                node.name = f"{node.name}_0"
    return tree

def main():
    args = parse_arguments()

    # read data
    data = []
    if args.simphy:
        data.extend(extract_simphy_data(args.simphy, args))
    if args.dupcoal:
        data.extend(extract_csv_data(args.dupcoal, "dupcoal", "duplications", "trees.tre", args))
    if args.mlmsc:
        data.extend(extract_csv_data(args.mlmsc, "mlmsc", "n_d", "gene_tree_*.newick", args))

    if args.discordance:
        df = pd.DataFrame(data, columns=["program", "Ne", "lb", "ld", "c", "duplications", "gene_copies", "observed_duplications", "total_quartets", "quartet_score", "normalized_quartet_score"])
    else:
        df = pd.DataFrame(data, columns=["program", "Ne", "lb", "ld", "c", "duplications", "gene_copies", "observed_duplications"])

    df.to_csv(f"{args.output_prefix}.csv", index=False)
    df = pd.read_csv(filepath_or_buffer=f"{args.output_prefix}.csv")
    df["lb"] = df["lb"].astype(float)
    df["Ne"] = df["Ne"].astype(float)
    df["ld"] = df["ld"].astype(str)
    df["c"] = df["c"].astype(float)
    df["duplications"] = df["duplications"].astype(int)
    df["gene_copies"] = df["gene_copies"].astype(int)
    if args.discordance:
        df["total_quartets"] = df["total_quartets"].astype(int)
        df["quartet_score"] = df["quartet_score"].astype(int)
        df["normalized_quartet_score"] = df["normalized_quartet_score"].astype(float)


    if args.discordance:
        plot_astral(df, args)
    else:
        plot_summaries(df, args)


if __name__ == "__main__":
    main()
