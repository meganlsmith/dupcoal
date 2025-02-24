#!/bin/bash
# conda activate dupcoal



#!/bin/bash

# Define the parameter values
lb_values=(2)
c_values=(1)
ld_values=(1)

# Path to SimPhy executable
DUPCOAL_BIN="$HOME/Documents/GitHub/dupcoal/dupcoal_v3.0.0.py"

for lb in "${lb_values[@]}"; do

    # adjust lb
    lb_adj=$(echo "$lb * 1e-3" | bc)

    for ld in "${ld_values[@]}"; do

        # adjust ld
        ld_adj=$(echo "$ld * $lb_adj" | bc)
        # set ld_name
        if [ "$ld" -eq 1 ]; then
            ld_name="lb"
        else
            ld_name="$ld"
        fi

        for c in "${c_values[@]}"; do
            sp_tree="fungi/fungi_corrected_dupcoal_${c}.newick"


            # Output directory name
            output_dir="dupcoal_fungi_results_disco/dupcoal_1_${lb}_${ld_name}_${c}"

            # Run SimPhy
            python3 $DUPCOAL_BIN --stree ${sp_tree} --mu_par ${ld_adj} --lambda_par ${lb_adj} --reps 1000 --output ${output_dir}
        done
    done
done