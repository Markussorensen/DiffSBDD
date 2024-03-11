#!/bin/bash

PROTEIN="kras"
PDB_ID="8AZV"
MODEL="moad"
RESI_LIST="A:7 A:9 A:10 A:11 A:58 A:60 A:61 A:62 A:63 A:64 A:65 A:68 A:69 A:72 A:78 A:88 A:92 A:95 A:96 A:99 A:100 A:102 A:103"
JOB_ID_START=261
MODEL_CHECKPOINT="/data/user/mas/DiffSBDD/models/models/${MODEL}_fullatom_joint.ckpt"
PDB_FILE="/data/user/mas/DiffSBDD/data/cancer/$PROTEIN/$PDB_ID.pdb"
OUT_DIR="/data/user/mas/DiffSBDD/data/cancer/$PROTEIN/results/$PDB_ID/$MODEL/"
BATCH_SIZE=16
N_SAMPLES=16
N_CORES=1
GPU_ID=1
N_JOBS_IN_PARALLEL=7
NUM_NODES_START=5
NUM_NODES_END=80

BASE_DIR="/data/user/mas/DiffSBDD/"

# Divide the number og jobs by 4. If it is not divisible by 4, add 1 to the first.

queued_jobs=$((JOB_ID_START-1))
parallel_job=1

# Go to the base directory
cd $BASE_DIR

echo "Queued jobs: $queued_jobs"

# Make out directory if it does not exist
mkdir -p $OUT_DIR

while [ $parallel_job -le $N_JOBS_IN_PARALLEL ]; do
    if [ $((parallel_job % 2)) -eq 0 ]; then
        jobs_run=0
        for ((num_nodes = (($NUM_NODES_START - $parallel_job)); num_nodes <= $NUM_NODES_END; num_nodes+=$N_JOBS_IN_PARALLEL)); do
            if (($num_nodes < $NUM_NODES_START)); then
                continue
            fi
            echo "Generating ligands for num_nodes_lig = $num_nodes"
            # Use tsp to queue the job
            if (($jobs_run == 0)); then
                tsp -N $N_CORES python generate_ligands.py $MODEL_CHECKPOINT \
                    --gpu_id $GPU_ID \
                    --pdbfile $PDB_FILE \
                    --outdir $OUT_DIR \
                    --resi_list $RESI_LIST \
                    --batch_size $BATCH_SIZE \
                    --n_samples $N_SAMPLES \
                    --num_nodes_lig $num_nodes \
                    --sanitize \
                    --relax \
                    --all_frags
            else
                tsp -N $N_CORES -D $queued_jobs python generate_ligands.py $MODEL_CHECKPOINT \
                    --gpu_id $GPU_ID \
                    --pdbfile $PDB_FILE \
                    --outdir $OUT_DIR \
                    --resi_list $RESI_LIST \
                    --batch_size $BATCH_SIZE \
                    --n_samples $N_SAMPLES \
                    --num_nodes_lig $num_nodes \
                    --sanitize \
                    --relax \
                    --all_frags
            fi
            # Increment the number of queued jobs
            ((jobs_run++))
            ((queued_jobs++))
        done
        
    else
        #Getting right starting index when descending
        INC_START_DES=$((NUM_NODES_START - parallel_job))
        while [ $INC_START_DES -le $NUM_NODES_END ]; do
            ((INC_START_DES+=$N_JOBS_IN_PARALLEL))
        done
        jobs_run=0
        for ((num_nodes = $INC_START_DES; num_nodes >= $NUM_NODES_START; num_nodes-=$N_JOBS_IN_PARALLEL)); do
            if (($num_nodes > $NUM_NODES_END)); then
                continue
            fi
            echo "Generating ligands for num_nodes_lig = $num_nodes"
            # Use tsp to queue the job
            if (( $jobs_run == 0 )); then
                tsp -N $N_CORES python generate_ligands.py $MODEL_CHECKPOINT \
                    --gpu_id $GPU_ID \
                    --pdbfile $PDB_FILE \
                    --outdir $OUT_DIR \
                    --resi_list $RESI_LIST \
                    --batch_size $BATCH_SIZE \
                    --n_samples $N_SAMPLES \
                    --num_nodes_lig $num_nodes \
                    --sanitize \
                    --relax \
                    --all_frags
            else
                tsp -N $N_CORES -D $queued_jobs python generate_ligands.py $MODEL_CHECKPOINT \
                    --gpu_id $GPU_ID \
                    --pdbfile $PDB_FILE \
                    --outdir $OUT_DIR \
                    --resi_list $RESI_LIST \
                    --batch_size $BATCH_SIZE \
                    --n_samples $N_SAMPLES \
                    --num_nodes_lig $num_nodes \
                    --sanitize \
                    --relax \
                    --all_frags
            fi
            # Increment the number of queued jobs
            ((jobs_run++))
            ((queued_jobs++))
        done
    fi
    ((parallel_job++))
done