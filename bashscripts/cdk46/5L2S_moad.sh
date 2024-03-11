#!/bin/bash

MODEL_CHECKPOINT="/data/user/mas/DiffSBDD/models/models/moad_fullatom_joint.ckpt"
PDB_FILE="/data/user/mas/DiffSBDD/data/cancer/cdk46/5L2S_cdk6_only.pdb"
OUT_DIR="/data/user/mas/DiffSBDD/data/cancer/cdk46/results/5L2S/moad/"
RESI_LIST="A:18 A:19 A:20 A:21 A:24 A:27 A:41 A:43 A:77 A:98 A:99 A:100 A:101 A:102 A:103 A:104 A:107 A:149 A:150 A:152 A:162 A:163"
JOB_ID_START=413
BATCH_SIZE=16
N_SAMPLES=16
N_CORES=2
GPU_ID=1
N_JOBS_IN_PARALLEL=4
NUM_NODES_START=5
NUM_NODES_END=80

BASE_DIR="/data/user/mas/DiffSBDD/"

# Divide the number og jobs by 4. If it is not divisible by 4, add 1 to the first.

queued_jobs=$((JOB_ID_START-1))
parallel_job=1

# Go to the base directory
cd $BASE_DIR

echo "Queued jobs: $queued_jobs"

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