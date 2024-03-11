PROTEIN="bcl2"
PDB_IDS=("6O0K")
PROTEIN_DIR="/data/user/mas/DiffSBDD/data/cancer/$PROTEIN/"
OUT_DIR="/data/user/mas/DiffSBDD/data/cancer/$PROTEIN/conformations/"
BASE_DIR="/data/user/mas/DiffSBDD/"

# Go to the base directory
cd $BASE_DIR

# Make out directory if it does not exist
mkdir -p $OUT_DIR

#Make sure Molkit can be imported adding the directory to the enironment path for the bash script
export PYTHONPATH=$PYTHONPATH:/data/user/mas/DiffSBDD/autodock/autodocktools-prepare-py3k

# Run the prepare_pdb_split_alt_confs for all proteins
for pdb_id in "${PDB_IDS[@]}"; do
    tsp -N 1 python "/data/user/mas/DiffSBDD/autodock/autodocktools-prepare-py3k/AutoDockTools/Utilities24/prepare_pdb_split_alt_confs.py" \
        -r "${PROTEIN_DIR}${pdb_id}.pdb" \
        -o "${OUT_DIR}${pdb_id}"
done
