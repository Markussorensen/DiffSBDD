#!/bin/bash

PROTEIN="cdk46"
DOCKING_TARGETS=("7SJ3" "5L2S")
OUTDIR="/data/user/mas/DiffSBDD/data/cancer/$PROTEIN/"
BASE_DIR="/data/user/mas/DiffSBDD/"
MOL_DIRS=()

for target in "${DOCKING_TARGETS[@]}"; do
    MOL_DIRS+=("/data/user/mas/DiffSBDD/data/cancer/$PROTEIN/results/$target/moad/")
done

# Go to the base directory
cd $BASE_DIR

# Run the evaluate script

tsp -N 1 python evaluate_new.py \
    --outdir $OUTDIR \
    --mol_dirs "${MOL_DIRS[@]}" \
    --docking_targets "${DOCKING_TARGETS[@]}"