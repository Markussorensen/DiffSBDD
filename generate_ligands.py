import argparse
from pathlib import Path
from tqdm import tqdm

import torch

import utils
from lightning_modules import LigandPocketDDPM


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('checkpoint', type=Path)
    parser.add_argument('--pdbfile', type=str)
    parser.add_argument('--resi_list', type=str, nargs='+', default=None)
    parser.add_argument('--ref_ligand', type=str, default=None)
    parser.add_argument('--outdir', type=Path)
    parser.add_argument('--n_samples', type=int, default=20)
    parser.add_argument('--batch_size', type=int, default=None)
    parser.add_argument('--num_nodes_lig', type=int, default=None)
    parser.add_argument('--all_frags', action='store_true')
    parser.add_argument('--sanitize', action='store_true')
    parser.add_argument('--relax', action='store_true')
    parser.add_argument('--resamplings', type=int, default=10)
    parser.add_argument('--jump_length', type=int, default=1)
    parser.add_argument('--timesteps', type=int, default=None)
    parser.add_argument('--gpu_id', type=int, default=1)
    args = parser.parse_args()

    pdb_id = Path(args.pdbfile).stem

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    torch.cuda.set_device(args.gpu_id)

    print(f"Using device {torch.cuda.get_device_name(args.gpu_id)}")

    if args.batch_size is None:
        args.batch_size = args.n_samples
    assert args.n_samples % args.batch_size == 0

    # Load model
    model = LigandPocketDDPM.load_from_checkpoint(
        args.checkpoint, map_location=device)
    model = model.to(device)

    if args.num_nodes_lig is not None:
        num_nodes_lig = torch.ones(args.n_samples, dtype=int) * \
                        args.num_nodes_lig
    else:
        num_nodes_lig = None

    molecules = []
    for i in tqdm(range(args.n_samples // args.batch_size)):
        try:
            molecules_batch = model.generate_ligands(
                args.pdbfile, args.batch_size, args.resi_list, args.ref_ligand,
                num_nodes_lig, args.sanitize, largest_frag=not args.all_frags,
                relax_iter=(200 if args.relax else 0),
                resamplings=args.resamplings, jump_length=args.jump_length,
                timesteps=args.timesteps) 
            molecules.extend(molecules_batch)
        except RuntimeError as e:
            print(f'Error in batch {i}: {e}')
            continue

    # Make SDF files
    utils.write_sdf_file(Path(args.outdir, f'{pdb_id}_mol_size{args.num_nodes_lig}.sdf'), molecules)
