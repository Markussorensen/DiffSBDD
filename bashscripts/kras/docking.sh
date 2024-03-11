#!/bin/bash

PROTEIN="kras"
PDB_IDS=("5VPI" "6OIM" "7OK3" "7RT1" "7TLG" "8AZZ" "8CX5" "8B00")
PDB_BINDING_POCKETS=("A:10 A:11 A:12 A:13 A:14 A:15 A:16 A:17 A:18 A:28 A:29 A:30 A:31 A:32 A:33 A:34 A:57 A:58 A:59 A:116 A:117 A:118 A:119 A:120 A:144 A:145 A:146 A:147" \
                     "A:9 A:10 A:11 A:12 A:13 A:16 A:34 A:58 A:59 A:60 A:61 A:62 A:63 A:68 A:69 A:72 A:92 A:95 A:96 A:99 A:100 A:103" \
                     "A:11 A:12 A:13 A:14 A:15 A:16 A:17 A:18 A:19 A:28 A:29 A:30 A:31 A:32 A:57 A:59 A:116 A:117 A:118 A:119 A:120 A:144 A:145 A:146 A:147" \
                     "A:9 A:10 A:11 A:12 A:16 A:59 A:60 A:61 A:62 A:63 A:64 A:65 A:68 A:69 A:72 A:78 A:92 A:95 A:96 A:99 A:100 A:102 A:103" \
                     "A:9 A:10 A:11 A:12 A:13 A:16 A:34 A:58 A:59 A:60 A:61 A:62 A:63 A:64 A:65 A:68 A:69 A:72 A:78 A:88 A:92 A:95 A:96 A:99 A:100 A:102 A:103" \
                     "A:7 A:9 A:10 A:11 A:12 A:58 A:60 A:61 A:62 A:63 A:64 A:65 A:68 A:69 A:72 A:78 A:88 A:92 A:95 A:96 A:99 A:100 A:102 A:103" \
                     "A:8 A:9 A:10 A:11 A:12 A:16 A:34 A:58 A:59 A:60 A:61 A:62 A:63 A:64 A:68 A:69 A:72 A:78 A:88 A:92 A:95 A:96 A:99 A:100 A:102 A:103" \
                     "A:7 A:9 A:10 A:11 A:58 A:60 A:61 A:62 A:63 A:64 A:65 A:68 A:69 A:72 A:78 A:92 A:95 A:96 A:99 A:100 A:102 A:103" \
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