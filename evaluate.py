import argparse
from pathlib import Path

import torch
import ast
from rdkit import Chem

import pandas as pd

from analysis.metrics import MoleculeProperties
from utils import read_sdf_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('mol_file', type=str)
    parser.add_argument('--mol_dir', type=str, default=None)
    parser.add_argument('--docking_results', type=str, default=None)
    parser.add_argument('--outdir', type=Path, default=None)
    args = parser.parse_args()

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    
    # Load SDF file
    if args.mol_dir is not None:
        sdf_files = list(Path(args.mol_dir).glob('[!.]*.sdf'))
    else:
        molecules = read_sdf_file(args.mol_file)

    # Load docking results
    if args.docking_results is not None:
        docking_results = pd.read_csv(args.docking_results)
        if args.mol_dir is not None:
            filenames = docking_results.iloc[:,2]
            generated = [True if not 'knowninhib' in filename else False for filename in filenames]
            docking_results["generated"] = generated
            filenames = []
            scores = []
            generated = []
            for i, row in docking_results.iterrows():
                score = row["scores"]
                filename = row["ligand"]
                gen = row["generated"]
                score = score.replace("nan", "None")
                score = ast.literal_eval(score)
                scores.append(score)
                filenames.append(filename)
                generated.append(gen)
        else:
            docking_results = docking_results["scores"][0]
            docking_results = docking_results.replace("nan", "None")
            docking_results = ast.literal_eval(docking_results)

    new_docking_results = []
    sdf_indexes = []
    new_generated = []
    new_filenames = []

    # Find valid molecules
    valid_molecules = []
    for k, filename in enumerate(filenames):
        molecules = read_sdf_file(Path(filename), sanitize=True)
        score = scores[k]
        gen = generated[k]
        if gen == False and len(score) == 1:
            score = [score[0] for i in range(len(molecules))]
        #Do not include directory in filename
        name = filename.split("/")[-1]
        for i, mol in enumerate(molecules):
            if mol is not None:
                valid_molecules.append(mol)
                sdf_indexes.append(i)
                new_generated.append(gen)
                new_filenames.append(name)
                new_docking_results.append(score[i])

    smiles = [Chem.MolToSmiles(mol) for mol in valid_molecules]
    
    print(f"Found {len(valid_molecules)}/{len(molecules)} valid molecules")

    mol_metrics = MoleculeProperties()

    valid_molecules = [valid_molecules]

    all_qed, all_sa, all_logp, all_lipinski, all_lipinski_violations, all_lipinski_values = mol_metrics.evaluate(valid_molecules, return_lipinski_violations=True, return_diversity=False)
    # all_qed, all_sa, all_logp, all_lipinski, per_pocket_diversity = mol_metrics.evaluate(valid_molecules, return_lipinski_violations=True, return_diversity=False)

    #Getting values out of the lipinski vlaues.
    all_weight = [x[0] for x in all_lipinski_values[0]]
    all_hdonor = [x[1] for x in all_lipinski_values[0]]
    all_hacceptor = [x[2] for x in all_lipinski_values[0]]
    all_rotatable = [x[4] for x in all_lipinski_values[0]]

    #Save metrics as csv
    metrics = {'Smile': smiles, 
               'QVinaScore': new_docking_results, 
               'sdfFileIdx':sdf_indexes,
               'Filename': new_filenames,
               'Generated': new_generated,
               'QED': all_qed[0], 
               'SA': all_sa[0], 
               'LogP': all_logp[0],
               'Weight': all_weight,
               'HDonor': all_hdonor,
               'HAcceptor': all_hacceptor,
               'RotatableBonds': all_rotatable,
               'Lipinski': all_lipinski[0], 
               'Lipinski_violation_rule1': [x[0] for x in all_lipinski_violations[0]],
               'Lipinski_violation_rule2': [x[1] for x in all_lipinski_violations[0]],
               'Lipinski_violation_rule3': [x[2] for x in all_lipinski_violations[0]],
               'Lipinski_violation_rule4': [x[3] for x in all_lipinski_violations[0]],
               'Lipinski_violation_rule5': [x[4] for x in all_lipinski_violations[0]],
            #    'Diversity': per_pocket_diversity[0],
               }
    metrics_df = pd.DataFrame(metrics)
    metrics_df = metrics_df.sort_values(by=['QVinaScore'], ascending=True)

    if args.outdir is not None:
        metrics_df.to_csv(Path(args.outdir, 'metrics.csv'))
    else:
        if args.mol_dir is not None:
            out_file = Path(args.mol_dir).stem + '_metrics.csv'
        else:
            out_file = Path(args.mol_file).stem + '_metrics.csv'
            out_file = Path(args.mol_file).parent / out_file
        metrics_df.to_csv(out_file)

