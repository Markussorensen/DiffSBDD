#!/bin/bash

PROTEIN="bcl2"
PDB_IDS=("6O0K" "6O0L" "6O0M" "6O0P")
PDB_BINDING_POCKETS=("A:99 A:100 A:101 A:103 A:104 A:105 A:107 A:108 A:111 A:112 A:115 A:133 A:136 A:137 A:143 A:144 A:145 A:146 A:148 A:149 A:152 A:153 A:156 A:198 A:202" \
                     "A:99 A:100 A:101 A:103 A:104 A:107 A:108 A:111 A:112 A:115 A:133 A:136 A:137 A:143 A:144 A:145 A:146 A:148 A:149 A:152 A:153 A:156 A:198 A:202" \
                     "A:99 A:100 A:101 A:103 A:104 A:107 A:108 A:111 A:112 A:115 A:133 A:136 A:137 A:143 A:144 A:145 A:146 A:148 A:149 A:152 A:153 A:156 A:198 A:202" \
                     "A:99 A:100 A:101 A:103 A:104 A:107 A:108 A:111 A:112 A:114 A:115 A:133 A:136 A:137 A:143 A:144 A:145 A:146 A:148 A:149 A:152 A:153 A:156 A:198 A:202" \
                     )
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