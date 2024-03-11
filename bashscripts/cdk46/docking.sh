#!/bin/bash

PROTEIN="cdk46"
PDB_IDS=("5L2S" "7SJ3")
PDB_BINDING_POCKETS=("A:18 A:19 A:20 A:21 A:24 A:27 A:41 A:43 A:77 A:98 A:99 A:100 A:101 A:102 A:103 A:104 A:107 A:149 A:150 A:152 A:162 A:163" \
                     "A:11 A:12 A:13 A:17 A:20 A:33 A:35 A:56 A:72 A:93 A:94 A:95 A:96 A:97 A:98 A:99 A:102 A:144 A:145 A:147 A:157 A:158")
MODEL="moad"
PDBQT_DIR="/data/user/mas/DiffSBDD/data/cancer/$PROTEIN/"

BASE_DIR="/data/user/mas/DiffSBDD/"

cd $BASE_DIR

# Loop over the PDB IDs and queue the docking jobs.
counter2=0
for pdb_id_gen in "${PDB_IDS[@]}"; do
    counter=0
    pdb_length=${#PDB_IDS[@]}
    while [ $counter -lt $pdb_length ]; do
        # Defining the dirs
        pdb_id_pdbfile=${PDB_IDS[$counter]}
        binding_pocket=${PDB_BINDING_POCKETS[$counter]}
        SDF_DIR="/data/user/mas/DiffSBDD/data/cancer/$PROTEIN/results/$pdb_id_gen/$MODEL/"
        OUT_DIR="/data/user/mas/DiffSBDD/data/cancer/$PROTEIN/results/$pdb_id_gen/$MODEL/sdf_out_$pdb_id_pdbfile/"
        # Use tsp to queue the job
        tsp -N 1 python analysis/docking.py \
                --pdbqt_dir $PDBQT_DIR \
                --sdf_dir $SDF_DIR \
                --out_dir $OUT_DIR \
                --new_receptor_target $pdb_id_pdbfile \
                --new_target_binding_pocket $binding_pocket \
                --write_csv \
                --write_dict
        ((counter++))
    done
    # Also dock for the known inhibitors
    binding_pocket=${PDB_BINDING_POCKETS[$counter2]}
    SDF_DIR="/data/user/mas/DiffSBDD/data/cancer/$PROTEIN/inhibitors/"
    OUT_DIR="/data/user/mas/DiffSBDD/data/cancer/$PROTEIN/results/$pdb_id_gen/$MODEL/sdf_out_inhibitors/"
    #Make dir if it doesn't exist
    mkdir -p $OUT_DIR
    tsp -N 1 python analysis/docking.py \
            --pdbqt_dir $PDBQT_DIR \
            --sdf_dir $SDF_DIR \
            --out_dir $OUT_DIR \
            --new_receptor_target $pdb_id_gen \
            --new_target_binding_pocket $binding_pocket \
            --write_csv \
            --write_dict
    ((counter2++))
done