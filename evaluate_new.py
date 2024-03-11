import argparse
from pathlib import Path

import torch
import ast
from rdkit import Chem
import numpy as np

import pandas as pd

from analysis.metrics import MoleculeProperties
from utils import read_sdf_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--mol_dirs', type=str, nargs='+', default=None)
    parser.add_argument('--docking_targets', type=str, nargs='+', default=None)
    parser.add_argument('--outdir', type=Path, default=None)
    args = parser.parse_args()

    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    # Receptor name that the molecule has been generated for
    all_receptor_names_gen = []
    # Receptor name that the molecule has been docked to
    all_receptor_names_dock = []
    all_metric_dfs = []
    for j, mol_dir in enumerate(args.mol_dirs):
        # Load SDF file
        sdf_files = list(Path(mol_dir).glob('[!.]*.sdf'))
        # Add inhibitor sdf files
        receptor_name_gen = mol_dir.split("/")
        results_idx = [i for i, s in enumerate(receptor_name_gen) if 'results' in s][0]
        receptor_name_gen = receptor_name_gen[results_idx+1]
        for i, docking_target in enumerate(args.docking_targets):
            docking_result_path = Path(mol_dir, f'sdf_out_{docking_target}/qvina2_scores.csv')
            docking_result = pd.read_csv(docking_result_path)
            filenames = docking_result.iloc[:,2]
            generated = [True if not 'knowninhib' in filename else False for filename in filenames]
            docking_result["generated"] = generated
            filenames = []
            scores = []
            generated = []
            for i, row in docking_result.iterrows():
                score = row["scores"]
                filename = row["ligand"]
                gen = row["generated"]
                score = score.replace("nan", "None")
                score = ast.literal_eval(score)
                scores.append(score)
                filenames.append(filename)
                generated.append(gen)

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
            all_psa = [x[5] for x in all_lipinski_values[0]]

            #Save metrics as csv
            metrics = { f'Smile': smiles, 
                        f'QVinaScore{docking_target}': new_docking_results, 
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
                        'PSA': all_psa,
                        'Lipinski': all_lipinski[0], 
                        'Lipinski_rule1_weight': [x[0] for x in all_lipinski_violations[0]],
                        'Lipinski_rule2_hdonors': [x[1] for x in all_lipinski_violations[0]],
                        'Lipinski_rule3_hacceptors': [x[2] for x in all_lipinski_violations[0]],
                        'Lipinski_rule4_logp': [x[3] for x in all_lipinski_violations[0]],
                        'Lipinski_rule5_rotatablebonds': [x[4] for x in all_lipinski_violations[0]],
                        'Lipinski_rule6_psa': [x[5] for x in all_lipinski_violations[0]],
                        #    'Diversity': per_pocket_diversity[0],
                        }
            metrics_df = pd.DataFrame(metrics)
            metrics_df["ReceptorGen"] = receptor_name_gen
            all_metric_dfs.append(metrics_df)
            all_receptor_names_gen.append(receptor_name_gen)
            all_receptor_names_dock.append(docking_target)

        # Do it for the known inhibitors
        docking_target = receptor_name_gen
        docking_result_path = Path(mol_dir, f'sdf_out_inhibitors/qvina2_scores.csv')
        docking_result = pd.read_csv(docking_result_path)
        filenames = docking_result.iloc[:,2]
        generated = [True if not 'knowninhib' in filename else False for filename in filenames]
        docking_result["generated"] = generated
        filenames = []
        scores = []
        generated = []
        for i, row in docking_result.iterrows():
            score = row["scores"]
            filename = row["ligand"]
            gen = row["generated"]
            score = score.replace("nan", "None")
            score = ast.literal_eval(score)
            scores.append(score)
            filenames.append(filename)
            generated.append(gen)

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
        all_psa = [x[5] for x in all_lipinski_values[0]]

        #Save metrics as csv
        metrics = { f'Smile': smiles, 
                    f'QVinaScore{docking_target}': new_docking_results, 
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
                    'PSA': all_psa,
                    'Lipinski': all_lipinski[0], 
                    'Lipinski_rule1_weight': [x[0] for x in all_lipinski_violations[0]],
                    'Lipinski_rule2_hdonors': [x[1] for x in all_lipinski_violations[0]],
                    'Lipinski_rule3_hacceptors': [x[2] for x in all_lipinski_violations[0]],
                    'Lipinski_rule4_logp': [x[3] for x in all_lipinski_violations[0]],
                    'Lipinski_rule5_rotatablebonds': [x[4] for x in all_lipinski_violations[0]],
                    'Lipinski_rule6_psa': [x[5] for x in all_lipinski_violations[0]],
                    #    'Diversity': per_pocket_diversity[0],
                    }
        
        metrics_df = pd.DataFrame(metrics)
        metrics_df["ReceptorGen"] = receptor_name_gen
        all_metric_dfs.append(metrics_df)
        all_receptor_names_gen.append(receptor_name_gen)
        all_receptor_names_dock.append(docking_target)

    # First merge the dataframes which have molecules generated for the same receptor
    metrics_df = {f'{rec_name}': None for rec_name in np.unique(all_receptor_names_gen)}
    for i, df in enumerate(all_metric_dfs):
        rec_name = all_receptor_names_gen[i]
        if metrics_df[rec_name] is None:
            metrics_df[rec_name] = df
        else:
            # Get common columns
            common_cols = list(set(metrics_df[rec_name].columns).intersection(set(df.columns)))
            # Merge the dataframes based on common columns other than 'Smile' and 'QVinaScore'
            merged_df = pd.merge(metrics_df[rec_name], df, on=common_cols, how='outer')
            metrics_df[rec_name] = merged_df

    # Now append the dataframes together
    final_df = pd.concat([df for df in metrics_df.values()], ignore_index=True)

    #Order the columns
    columns = list(final_df.columns)
    columns.remove("Filename")
    columns.remove("sdfFileIdx")
    columns.remove("Lipinski")
    columns.remove("Smile")
    columns.insert(0, "Filename")
    columns.insert(1, "sdfFileIdx")
    for i, rec in enumerate(np.unique(all_receptor_names_dock)):
        columns.remove(f'QVinaScore{rec}')
        columns.insert(i+2, f'QVinaScore{rec}')
    columns.insert(len(set(all_receptor_names_dock))+2, "Lipinski")
    final_df = final_df[columns]
    #Order by docking score
    final_df = final_df.sort_values(by=[f'QVinaScore'+pdb for pdb in np.unique(all_receptor_names_dock)], ascending=True)

    # Save the dataframe
    final_df.to_csv(args.outdir / 'molecule_metrics.csv', index=False)
