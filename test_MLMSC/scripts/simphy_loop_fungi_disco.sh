#!/bin/bash
# Define the parameter values
lb_values=(2)
c_values=(1)
ld_values=(1)

# Path to SimPhy executable
SIMPHY_BIN="$HOME/Documents/programs/SimPhy_1.0.2/bin/simphy_mac64"
    
for lb in "${lb_values[@]}"; do


    # adjust lb
    lb_adj=$(echo "$lb * 0.0000000001" | bc)

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

            sp_tree=$(cat fungi/fungi_simphy_${c}.newick)


            # Output directory name
            output_dir="simphy_fungi_results_disco/simphy_1_${lb}_${ld_name}_${c}"

            # adjust sp size
            sp_adj=$(echo "1 * 10000000" | bc)

            # adjust gen time
            n_gen=$(echo "scale=10; 10 / 1" | bc -l)

            # Run SimPhy
            $SIMPHY_BIN -S "$sp_tree" -LB f:${lb_adj} -LD f:${ld_adj} -RL f:1000 -SP f:${sp_adj} -OD 1 -O "$output_dir" -SG f:${n_gen}


            # Concatenate g_trees1.trees files
            find "$output_dir" -mindepth 2 -maxdepth 2 -type f -name "g_trees*.trees" | sort | xargs cat > "$output_dir/g_trees.tre"
        done
    done
done