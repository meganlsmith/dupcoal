#!/bin/bash

# Define the parameter values
lb_values=(0 1 2 3 4 5 6)
c_values=(1 9)
ld_values=(0.5 1)

# Path to SimPhy executable
MLMLSC_BIN="$HOME/Documents/MLMSC/MLMSC-II/MLMSC.py"

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

            # adjust sp tree
            sp_tree="fungi/fungi_corrected_dupcoal_${c}.newick"


            # Output directory name
            output_dir="mlmsc_fungi_results/mlmsc_1_${lb}_${ld_name}_${c}"

            # Run SimPhy
            python ${MLMLSC_BIN} -i ${sp_tree} -l ${ld_adj} -d ${lb_adj} -n 10000 -t 0 -c 1
            mv output ${output_dir}

        done

    done

done