#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import string
from collections import OrderedDict

from configs import DATA_CONFIG
from variables import subfolders, nonstandard_amino_acids, hetatm_non_metal_ion, hetatm_metal_ion
from utils import (
    drop_unnamed_columns,
    parse_csv_bool,
    parse_csv_list,
    get_struct_variant_name,
)


def parse_prep_input_file(inputs, input_dir):
    """Read structure-preparation rows from CSV or accept a dict directly."""
    import pandas as pd

    if isinstance(inputs, dict):
        return inputs

    if isinstance(inputs, str) and '.csv' in inputs:
        inputs_df = pd.read_csv(os.path.join(input_dir, inputs))
        inputs_df = drop_unnamed_columns(inputs_df)
        print(inputs_df)

        prep_inputs = OrderedDict()
        for _, row in inputs_df.iterrows():
            if int(row['process_structure']) != 1:
                continue

            struct_name = row['struct_name']
            prep_inputs[struct_name] = {
                'pdb_id': row['pdb_id'] if 'pdb_id' in row.index else struct_name,
                'ligand_id': row['ligand_id'] if 'ligand_id' in row.index else None,
                'chain_id': parse_csv_list(row['chain_id']) if 'chain_id' in row.index else [],
                'target_chain': parse_csv_list(row['target_chain']) if 'target_chain' in row.index else [],
                'ligand_chain_id': parse_csv_list(row['ligand_chain_id']) if 'ligand_chain_id' in row.index else [],
                'keep_multiple_chains_in_struct': parse_csv_bool(
                    row['keep_multiple_chains_in_struct'] if 'keep_multiple_chains_in_struct' in row.index else None,
                    default=True,
                ),
            }
        return prep_inputs

    raise ValueError(f'Unsupported inputs argument: {inputs}')


def normalize_singleton_list(value):
    """Collapse single-item lists to the item itself for APIs that expect scalars."""
    if isinstance(value, list) and len(value) == 1:
        return value[0]
    return value


def pdb_to_sce(
        pdb_fpath,
        ligand_name,
        chains_to_process=None,
        ligand_chain_id=None,
        keep_ligand=True,
        keep_metal_ion=True,
        split_lig_obj=True,
        save_sce=True,
        skip_chains=None,
        keep_multiple_chains_in_struct=True,
        target_chains=None,
        struct_name=None,
        pdb_name=None,
):
    import yasara

    # initialize yasara
    yasara.info.mode = 'txt'
    yasara.info.licenseshown = 0
    yasara.Console('Off')

    # load pdb
    yasara.Clear()
    yasara.LoadPDB(pdb_fpath)
    yasara.DelWater()
    yasara.CleanAll()
    print('Parsing', os.path.basename(pdb_fpath).split('.')[0], '...')

    # get number of chains
    chain_list = list(set(yasara.NameMol('Obj 1')))
    chain_list.sort()
    if ligand_chain_id is not None:
        non_receptor_chains = [ligand_chain_id]
    else:
        non_receptor_chains = [chain for chain in chain_list if chain not in string.ascii_uppercase]
    receptor_chains = [chain for chain in chain_list if chain not in non_receptor_chains]
    print('Chain list (raw PDB):', len(chain_list), chain_list)
    print('Receptor chains:', receptor_chains)
    print('Non-receptor chains:', non_receptor_chains)

    # if there are overall more than one chain, but only one is based on a protein receptor
    # the 2nd chain must belong to the ligand -- merge this into the same chain as the protein
    if len(chain_list) > 1 and len(non_receptor_chains) > 0 and len(receptor_chains) == 1:
        yasara.NameMol('all', receptor_chains[0])
        chain_list = list(set(yasara.NameMol('Obj 1')))
        chain_list.sort()
        print('Chain list (after merging ligand with receptor chain):', len(chain_list), chain_list)

    if chains_to_process is None:
        chains_to_process = chain_list.copy()
    elif isinstance(chains_to_process, str):
        chains_to_process = list(chains_to_process)
    else:
        chains_to_process = list(chains_to_process)

    if skip_chains is None:
        skip_chains = []

    chains_to_process = [c for c in chains_to_process if c not in skip_chains]
    print('Final Chains to Process:', chains_to_process)

    if not target_chains:
        target_chains = chains_to_process.copy()
    elif isinstance(target_chains, str):
        target_chains = list(target_chains)
    else:
        target_chains = list(target_chains)
    target_chains = [c for c in target_chains if c in chains_to_process]
    print('Target Chains:', target_chains)

    if struct_name is None:
        struct_name = os.path.basename(pdb_fpath).replace('.pdb', '')
    if pdb_name is None:
        pdb_name = struct_name

    sce_fpath_list = []
    chains_to_save = [None] if keep_multiple_chains_in_struct else target_chains
    for chain in chains_to_save:
        if keep_multiple_chains_in_struct:
            print('Processing combined structure...')
        else:
            print('Processing Chain', chain, '...')

        yasara.Clear()
        yasara.LoadPDB(pdb_fpath)
        yasara.DelWater()
        yasara.CleanAll()

        chains_for_scene = chains_to_process if keep_multiple_chains_in_struct else [chain]
        delmol_cmd = ' and '.join([f'not {chain_id}' for chain_id in chains_for_scene])
        if delmol_cmd:
            yasara.DelMol(delmol_cmd)

        chain_list_byreceptor = list(set(yasara.NameMol('Obj 1')))
        print(chain_list_byreceptor)
        resname_notprotein = [r for r in yasara.NameRes('not Protein') if r not in nonstandard_amino_acids]
        print('Non-protein residues (before deletion):', resname_notprotein)

        if ligand_name is None:
            ligand_name = [r for r in resname_notprotein if r not in hetatm_non_metal_ion + hetatm_metal_ion][0]
        print('Ligand ID:', ligand_name)

        hetatm_res_to_keep = list(hetatm_non_metal_ion)
        if keep_metal_ion:
            hetatm_res_to_keep += hetatm_metal_ion
        delres_cmd = "not " + " and not ".join(["Protein"] + hetatm_res_to_keep + nonstandard_amino_acids)
        if keep_ligand:
            delres_cmd += ' and not ' + ligand_name

        print('Deletion cmd:', delres_cmd)
        yasara.DelRes(delres_cmd)
        resname_notprotein = [r for r in yasara.NameRes('not Protein') if r not in nonstandard_amino_acids]
        print('Non-Protein residues (after deletion)', len(resname_notprotein), resname_notprotein)

        if split_lig_obj:
            print('# of obj (before splitting Obj 1):', yasara.CountObj('all'))
            yasara.SplitObj(1)
            print('# of obj (after splitting Obj 1):', yasara.CountObj('all'))
            yasara.JoinObj(f'not Res {ligand_name}', 1)
            print('# of obj (after joining non-ligand objects to Obj 1):', yasara.CountObj('all'))
            obj_list = yasara.ListObj('all')
            obj_list.sort()
            ligand_objs = [obj for obj in obj_list if ligand_name in yasara.NameRes(f'Obj {obj}')]
            if ligand_objs:
                ligand_obj_anchor = ligand_objs[0]
                if ligand_obj_anchor != 2:
                    yasara.SwapObj(ligand_obj_anchor, 2)
                yasara.JoinObj(f'Res {ligand_name}', 2)
            obj_list = yasara.ListObj('all')
            print(obj_list)

        if save_sce:
            pdb_dir = os.path.dirname(pdb_fpath)
            pdb_parent_dir = os.path.dirname(pdb_dir)
            data_subdir = os.path.basename(pdb_dir)
            if data_subdir == 'pdb':
                sce_dir = os.path.join(pdb_parent_dir, 'sce')
            else:
                sce_dir = os.path.join(os.path.dirname(pdb_parent_dir), 'sce', data_subdir)
            struct_variant_name = get_struct_variant_name(
                struct_name,
                pdb_name=pdb_name,
                target_chain=chain,
                keep_multiple_chains_in_struct=keep_multiple_chains_in_struct,
            )
            sce_fname = struct_variant_name + '.sce'
            sce_fpath = os.path.join(sce_dir, sce_fname)
            if not os.path.exists(sce_dir):
                os.makedirs(sce_dir)
                print('Created postOpt sub-directory:', sce_dir)
            yasara.SaveSce(sce_fpath)
            print('Saved .sce file:', sce_fpath)
            sce_fpath_list.append(sce_fpath)

    return sce_fpath_list


class PrepareYasaraSCE:
    def __init__(
        self,
        data_folder,
        data_subfolder,
        input_subfolder='yasara/Input/',
    ):
        self.data_folder = data_folder
        self.data_subfolder = data_subfolder
        self.input_dir = os.path.join(data_folder, input_subfolder)
        self.pdb_dir = os.path.join(data_folder, subfolders['pdb'], data_subfolder)

    def validate_prepared_scene(self, sce_fpath):
        import yasara

        yasara.info.mode = 'txt'
        yasara.Console('Off')
        yasara.Clear()
        yasara.LoadSce(sce_fpath)

        obj_list = yasara.ListObj('all')
        obj_count = len(obj_list)
        issues = []
        if obj_count != 2:
            issues.append(f'expected 2 objects, found {obj_count}')

        return {
            'sce_fpath': sce_fpath,
            'obj_count': obj_count,
            'obj_list': obj_list,
            'issues': issues,
        }

    def prepare_structure(self, struct_name, prep_inputs):
        pdb_name = prep_inputs.get('pdb_id', struct_name)
        pdb_fpath = os.path.join(self.pdb_dir, pdb_name + '.pdb')
        chain_id = prep_inputs.get('chain_id') or None
        ligand_chain_id = normalize_singleton_list(prep_inputs.get('ligand_chain_id') or None)
        return pdb_to_sce(
            pdb_fpath,
            prep_inputs.get('ligand_id'),
            chains_to_process=chain_id,
            ligand_chain_id=ligand_chain_id,
            keep_ligand=True,
            keep_metal_ion=True,
            keep_multiple_chains_in_struct=prep_inputs.get('keep_multiple_chains_in_struct', True),
            target_chains=prep_inputs.get('target_chain'),
            struct_name=struct_name,
            pdb_name=pdb_name,
        )

    def run_pipeline(self, inputs):
        prep_inputs = parse_prep_input_file(inputs, self.input_dir)
        sce_fpaths = []
        scene_validation_issues = []
        for struct_name, struct_inputs in prep_inputs.items():
            print('Preparing', struct_name, '...')
            prepared_sce_fpaths = self.prepare_structure(struct_name, struct_inputs)
            sce_fpaths += prepared_sce_fpaths
            for sce_fpath in prepared_sce_fpaths:
                validation_result = self.validate_prepared_scene(sce_fpath)
                if validation_result['issues']:
                    scene_validation_issues.append(validation_result)
        print('Prepared SCE files:')
        for sce_fpath in sce_fpaths:
            print(sce_fpath)
        if scene_validation_issues:
            print('\nPrepared scenes with issues:')
            for validation_result in scene_validation_issues:
                print(
                    '-',
                    validation_result['sce_fpath'],
                    '>>',
                    '; '.join(validation_result['issues']),
                    f"(objects={validation_result['obj_list']})",
                )
        return sce_fpaths


if __name__ == "__main__":
    config = DATA_CONFIG.copy()
    data_folder = config['data_folder']
    data_subfolder = config['data_subfolder']
    inputs = config['inputs']

    prep = PrepareYasaraSCE(
        data_folder,
        data_subfolder,
    )
    prep.run_pipeline(inputs)
