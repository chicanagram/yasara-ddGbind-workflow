# -*- coding: utf-8 -*-
import yasara
import numpy as np
import os
import traceback
from collections import OrderedDict
from configs import DATA_CONFIG, RUN_CONFIG, RUN_PARAMS
from utils import (
    opsys,
    get_mutstr,
    split_mutation,
    findProcess,
    exit_program,
    save_dict_as_csv,
    combine_csv_files,
    drop_unnamed_columns,
    parse_csv_list,
    parse_csv_bool,
    get_struct_variant_name,
)

# CALCULATE BINDING ENERGY INCLUDING IMPLICIT SOLVATION
# Created by Sebastian Maurer-Stroh, inspired by template macros from Elmar Krieger.
# Edited by Jhoann Miyajima and Charmaine Chia

ENERGY_FEATURE_PREFIXES = ['ebind', 'epot', 'esol', 'esolcol', 'esolvdw', 'esurfacc']
ENERGY_COMPONENT_SUFFIXES = ['Cpx', 'Rtr', 'Lgd']
POSTOPT_SUBDIR = 'postOpt'
FAILED_JOBS_PREFIX = 'failed_jobs'
FAILED_JOB_COLUMNS = [
    'struct_name',
    'ligand_name',
    'target_chain',
    'mutation',
    'proc_num',
    'error_type',
    'error_message',
    'traceback',
]


class SceneValidationError(RuntimeError):
    pass


def parse_input_file(inputs, input_dir, sep):
    """Get struct name, ligand, and mutations to process from inputs file"""
    import pandas as pd
    from utils import mutstr_to_mutations

    def parse_mutation_entries(inputs_df_filt):
        mutations_ = inputs_df_filt['mutations'].tolist()
        mutations = []
        for mut_entry in mutations_:
            mut_list = mut_entry.split(',')
            for idx, mut in enumerate(mut_list):
                if not mut.isalnum():
                    combi_mut, _ = mutstr_to_mutations(mut, sep)
                    mut_list[idx] = combi_mut
            mutations += mut_list
        return mutations

    def parse_sequence_mutations(inputs_df, inputs_df_filt):
        from utils import list_all_mutations

        res_to_mutate = None
        wt_seq = inputs_df_filt.iloc[0]['sequence']
        if 'res' in inputs_df:
            res_to_mutate_ = inputs_df_filt['res'].tolist()
            res_to_mutate = []
            for res_entry in res_to_mutate_:
                res_list = [int(res) for res in res_entry.split(',')]
                res_to_mutate += res_list
            res_to_mutate = sorted(set(res_to_mutate))
        return list_all_mutations(wt_seq, ignore_mutations_to_WT=True, pos_to_mutate=res_to_mutate)

    def parse_structure_metadata(struct_name, inputs_df_filt):
        chain_id, target_chain = [], []
        pdb_id = struct_name
        keep_multiple_chains_in_struct = True
        if 'pdb_id' in inputs_df_filt.columns:
            pdb_id = inputs_df_filt.iloc[0]['pdb_id']
        if 'chain_id' in inputs_df_filt.columns:
            chain_id = parse_csv_list(inputs_df_filt.iloc[0]['chain_id'])
        if 'target_chain' in inputs_df_filt.columns:
            target_chain = parse_csv_list(inputs_df_filt.iloc[0]['target_chain'])
        if 'keep_multiple_chains_in_struct' in inputs_df_filt.columns:
            keep_multiple_chains_in_struct = parse_csv_bool(inputs_df_filt.iloc[0]['keep_multiple_chains_in_struct'], default=True)
        if not target_chain:
            target_chain = chain_id.copy()
        if chain_id and not set(target_chain).issubset(set(chain_id)):
            raise ValueError(f'target_chain must be a subset of chain_id for {(struct_name, ligand_name)}. Got target_chain={target_chain}, chain_id={chain_id}')
        return {
            'pdb_id': pdb_id,
            'chain_id': chain_id,
            'target_chain': target_chain,
            'keep_multiple_chains_in_struct': keep_multiple_chains_in_struct,
        }

    if isinstance(inputs, dict):
        struct_to_mutate_dict = inputs

    elif isinstance(inputs, str) and '.csv' in inputs:
        inputs_df = pd.read_csv(os.path.join(input_dir, inputs))
        inputs_df = drop_unnamed_columns(inputs_df)
        print(inputs_df)

        # Group rows by unique structure / ligand combinations so each
        # prepared scene can be paired with its mutation list.
        struct_to_mutate_dict = OrderedDict()
        struct_ligand_df = inputs_df[['struct_name', 'ligand_name']].drop_duplicates(ignore_index=True)
        struct_ligand_pairs = list(zip(struct_ligand_df['struct_name'].tolist(), struct_ligand_df['ligand_name'].tolist()))

        for struct_name, ligand_name in struct_ligand_pairs:

            # get filtered df
            inputs_df_filt = inputs_df.loc[
                          (inputs_df['struct_name'] == struct_name) &
                          (inputs_df['ligand_name'] == ligand_name) &
                          (inputs_df['process_structure'] == 1), :
                          ]
            if len(inputs_df_filt) == 0:
                print((struct_name, ligand_name), '>> No entries after filtering.')
                continue
            struct_metadata = parse_structure_metadata(struct_name, inputs_df_filt)

            # Either take explicit mutations from the CSV or enumerate
            # all allowed substitutions from a supplied sequence.
            if 'mutations' in inputs_df_filt.columns:
                mutations = parse_mutation_entries(inputs_df_filt)
            # no specific mutations provided
            elif 'mutations' not in inputs_df_filt.columns and 'sequence' in inputs_df_filt.columns:
                mutations = parse_sequence_mutations(inputs_df, inputs_df_filt)
            else:
                print('"mutations", "res", and "sequence" columns not provided in input CSV. Unable to determine mutations for', (struct_name, ligand_name))
                break

            # update struct_to_mutate_dict entry
            struct_to_mutate_dict[(struct_name, ligand_name)] = {
                'ligand_name': ligand_name,
                **struct_metadata,
                'mutations': mutations,
            }
    else:
        raise ValueError(f'Unsupported inputs argument: {inputs}')

    # print out inputs
    for struct_key, struct_dict in struct_to_mutate_dict.items():
        print(struct_key)
        print(struct_dict)

    return struct_to_mutate_dict


def parse_args_process_mutant(args, mutbindDDG, sep):
    mutation_for_log = None
    struct_name = None
    ligand_name = None
    target_chain = None
    proc_num = None
    if len(args) == 13:
        mutation, struct_name, ligand_name, target_chain, keep_multiple_chains_in_struct, move, minimize_energy, resetSce, nrep, ff, mvdist, mvdrug, surfout = args
    elif len(args) == 14:
        mutation, struct_name, ligand_name, target_chain, keep_multiple_chains_in_struct, move, minimize_energy, resetSce, nrep, ff, mvdist, mvdrug, surfout, proc_num = args
    else:
        raise ValueError(f'Unexpected process_mutant arg length: {len(args)}')
    try:
        mutstr, mutation = get_mutstr(mutation, sep=sep)
        mutation_for_log = mutstr
        print(proc_num, struct_name + " Mutation: " + mutstr + "; Target chain: " + str(target_chain) + "; Move: " + move + "; Minimize: " + str(minimize_energy))
        log_fpath_avg, log_fpath_full = mutbindDDG.process_mutant(
            mutation,
            struct_name,
            ligand_name,
            target_chain,
            keep_multiple_chains_in_struct,
            move,
            minimize_energy,
            resetSce,
            nrep,
            ff,
            mvdist,
            mvdrug,
            surfout,
            proc_num,
        )
        return log_fpath_avg, log_fpath_full
    except Exception as exc:
        mutbindDDG.log_failed_job(
            struct_name=struct_name,
            ligand_name=ligand_name,
            target_chain=target_chain,
            mutation=mutation_for_log if mutation_for_log is not None else str(mutation),
            proc_num=proc_num,
            exc=exc,
            traceback_text=traceback.format_exc(),
        )
        print(f'FAILED job for {(struct_name, target_chain, mutation_for_log)}: {exc}')
        return None, None

class MutBindDDG:
    def __init__(
            self,
            data_folder,
            data_subfolder,
            input_subfolder='yasara/Input/',
            output_subfolder='yasara/Output/',
            fix_metal_ion=True,
            cntions=0, # If counterions=1, counter ions will be implicitly considered by setting the net charges to 0
            JToUnit=1.43932620865201e20,
            sep='-',
            exit_program_after_processing=True,
            append_to_existing_output=True,
            energy_calc_method='BoundaryFast'
    ):
        self.data_folder = data_folder
        self.data_subfolder = data_subfolder
        self.sce_dir = os.path.join(data_folder, 'sce', data_subfolder)
        self.input_dir = os.path.join(data_folder, input_subfolder)
        self.output_dir = os.path.join(data_folder, output_subfolder, data_subfolder)
        self.fix_metal_ion = fix_metal_ion
        self.cntions = cntions
        self.JToUnit = JToUnit
        self.sep = sep
        self.exit_program = exit_program_after_processing
        self.append_to_existing_output = append_to_existing_output
        self.energy_calc_method = energy_calc_method

        # results columns
        metadata_list = ['struct', 'target_chain', 'mutations', 'ligname', 'setname', 'minimize_energy', 'resetSce', 'move', 'mvdist', 'mvdrug', 'ff', 'surfcost', 'counterions', 'nrep']
        self.energy_features_list = [
            feature + suffix
            for feature in ENERGY_FEATURE_PREFIXES
            for suffix in ['DDG', 'MT', 'WT'] + [state + element for state in ['MT', 'WT'] for element in ENERGY_COMPONENT_SUFFIXES]
        ]
        self.res_avg_cols = metadata_list + [f+suffix for f in self.energy_features_list for suffix in ['', '_std']]
        self.res_all_cols = metadata_list + ['n'] + self.energy_features_list

        # create output directory if it doesn't exist
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print('Created output directory:', self.output_dir)
        if not os.path.exists(self.postopt_dir):
            os.makedirs(self.postopt_dir)
            print('Created postOpt sub-directory:', self.postopt_dir)

        # initialize yasara
        self.initialize_yasara()

    def prepare_sce_files(self, all_inputs):
        scene_jobs = []

        for (struct_name, ligand_name), struct_inputs_dict in all_inputs.items():
            mutations = struct_inputs_dict['mutations']
            pdb_id = struct_inputs_dict['pdb_id']
            target_chain = struct_inputs_dict['target_chain']
            keep_multiple_chains_in_struct = struct_inputs_dict['keep_multiple_chains_in_struct']
            for chain in target_chain:
                scene_jobs.append({
                    'struct_name': get_struct_variant_name(
                        struct_name,
                        pdb_name=pdb_id,
                        target_chain=chain,
                        keep_multiple_chains_in_struct=keep_multiple_chains_in_struct,
                    ),
                    'ligand_name': ligand_name,
                    'target_chain': chain,
                    'keep_multiple_chains_in_struct': keep_multiple_chains_in_struct,
                    'mutations': mutations,
                })
        return scene_jobs

    def initialize_yasara(self, ff=None):
        yasara.info.mode = 'txt'
        yasara.Console('Off')
        yasara.Clear()
        yasara.Processors(1)
        yasara.EnergyUnit('kcal/mol')
        if ff is not None:
            yasara.ForceField(ff, setpar='yes')

    def save_sce_as_pdb(self, output_dir, struct_subdir='postOpt/', del_sce_files=False):
        # convert sce files to pdb files
        struct_dir = os.path.join(output_dir, struct_subdir.rstrip('/'))
        sce_list = [f for f in os.listdir(struct_dir) if f.endswith('.sce')]
        for sce_fname in sce_list:
            try:
                yasara.Clear()
                sce_fpath = os.path.join(struct_dir, sce_fname)
                yasara.LoadSce(sce_fpath)
                yasara.JoinObj('2', '1', center='No')
                yasara.SavePDB('1', os.path.join(struct_dir, sce_fname.replace('.sce', '.pdb')))
                print('Converted sce to pdb file:', sce_fname)
            except Exception as exc:
                print(f'Failed to export {sce_fname}: {exc}')
        if del_sce_files:
            for sce_fname in sce_list:
                os.remove(os.path.join(struct_dir, sce_fname))

    def get_log_fpath(self, struct_name, target_chain):
        return os.path.join(self.output_dir, f'DDG_{struct_name}_{target_chain}')

    @property
    def postopt_dir(self):
        return os.path.join(self.output_dir, POSTOPT_SUBDIR)

    def get_setname(self, ligname, struct, mutant, ff, resetSce, move, mvdist, mvdrug, surfout, cntions, nrep):
        return '{0}_{1}_{2}_{3}_rs{4}_{5}_d{6}_md{7}_s{8}_c{9}_r{10}'.format(
            ligname,
            struct,
            mutant,
            ff,
            resetSce * 1,
            move,
            mvdist,
            mvdrug,
            surfout,
            cntions,
            nrep,
        )
    def load_structure(self, struct_name):
        yasara.Clear()
        yasara.LoadSce(os.path.join(self.sce_dir, struct_name + '.sce'))
        yasara.CellAuto(extension='10')
        yasara.FixAll()

    def validate_loaded_scene(self, struct_name, target_chain, keep_multiple_chains_in_struct):
        obj_list = yasara.ListObj('all')
        if 1 not in obj_list:
            raise SceneValidationError(f'{struct_name}: missing Obj 1 after loading scene')
        if 2 not in obj_list:
            raise SceneValidationError(f'{struct_name}: missing Obj 2 after loading scene')

        if not keep_multiple_chains_in_struct:
            return

        receptor_selector = f'Obj 1 and Mol {target_chain}'
        ligand_selector = f'Obj 2 and Mol {target_chain}'
        complex_selector = f'Mol {target_chain}'

        if not yasara.NameRes(receptor_selector):
            raise SceneValidationError(f'{struct_name}: no receptor residues found for selector "{receptor_selector}"')
        if not yasara.NameRes(ligand_selector):
            raise SceneValidationError(f'{struct_name}: no ligand residues found for selector "{ligand_selector}"')
        if not yasara.NameRes(complex_selector):
            raise SceneValidationError(f'{struct_name}: no residues found for complex selector "{complex_selector}"')

    def get_failed_jobs_fpath(self, proc_num=None):
        suffix = 'serial' if proc_num is None else str(proc_num)
        return os.path.join(self.output_dir, f'{FAILED_JOBS_PREFIX}_{suffix}.csv')

    @staticmethod
    def sanitize_failed_job_value(value):
        return str(value).replace('\n', ' ').replace(',', ';')

    def log_failed_job(self, struct_name, ligand_name, target_chain, mutation, proc_num, exc, traceback_text):
        failed_jobs_fpath = self.get_failed_jobs_fpath(proc_num)
        file_exists = os.path.exists(failed_jobs_fpath)
        with open(failed_jobs_fpath, 'a') as f:
            if not file_exists:
                f.write(','.join(FAILED_JOB_COLUMNS) + '\n')
            row = [
                struct_name,
                ligand_name,
                target_chain,
                mutation,
                proc_num,
                type(exc).__name__,
                self.sanitize_failed_job_value(exc),
                self.sanitize_failed_job_value(traceback_text.replace('\n', ' | ')),
            ]
            f.write(','.join([self.sanitize_failed_job_value(value) for value in row]) + '\n')

    def combine_failed_job_logs(self, remove_combined_files=False):
        failed_job_fpaths = sorted(
            [
                os.path.join(self.output_dir, fname)
                for fname in os.listdir(self.output_dir)
                if fname.startswith(FAILED_JOBS_PREFIX + '_') and fname.endswith('.csv')
            ]
        )
        if not failed_job_fpaths:
            return None

        import pandas as pd

        df_all = None
        for failed_job_fpath in failed_job_fpaths:
            df = pd.read_csv(failed_job_fpath)
            if df_all is None:
                df_all = df.copy()
            else:
                df_all = pd.concat([df_all, df], axis=0, ignore_index=True)

        combined_fpath = os.path.join(self.output_dir, FAILED_JOBS_PREFIX + '.csv')
        df_all.to_csv(combined_fpath, index=False)

        if remove_combined_files:
            for failed_job_fpath in failed_job_fpaths:
                os.remove(failed_job_fpath)
        return combined_fpath

    def get_res_selector(self, position, target_chain=None):
        if target_chain in [None, '']:
            return f'Res {position}'
        return f'Mol {target_chain} and Res {position}'

    def set_up_sce_for_minimization(self, mutation, target_chain, move, mvdist, mvdrug, continue_processing=True):
        # Free only the local neighborhood around the mutation site
        # unless the run configuration explicitly allows broader motion.
        mutname = ''
        for mutgrp in mutation:
            WT, position, MT = split_mutation(mutgrp)
            # check if current residue position is the same
            res_selector = self.get_res_selector(position, target_chain)
            wt_matches = yasara.NameRes(res_selector)
            if not wt_matches:
                print('No residue found for selection:', res_selector)
                continue_processing = False
                break
            WT_actual = wt_matches[0]
            print('Actual WT residue:', WT_actual, '; Target WT residue:', WT.upper())
            if WT_actual != WT.upper():
                print('Incorrect WT amino acid found at position to mutate >> End processing')
                continue_processing = False
                break
            else:
                print('Correct WT amino acid found at position to mutate >> Continue processing')
            x = '{0}{1}{2}'.format(WT, position, MT)
            mutname += x
            yasara.ShowRes(res_selector)
            yasara.FreeAtom(f'{move} and Mol {target_chain} with distance <{mvdist} from {res_selector}')
        print('mutname:', mutname)

        if mvdrug == 1:  # allow ligand and close by residues to move
            yasara.FreeAtom(f'{move} and Mol {target_chain} with distance <{mvdist} from Obj 2 and Mol {target_chain}')
        if mvdrug == 0:  # fix ligand
            yasara.FixObj("2")
        # fix metal ion
        if self.fix_metal_ion:
            yasara.FixAtom('metal')
        return mutname, continue_processing

    def swapres(self, mutation, WT_or_MT, target_chain):
        # Build either the WT or MT state directly inside the current
        # YASARA scene before minimization / energy evaluation.
        for mutgrp in mutation:
            print(mutgrp)
            WT, position, MT = split_mutation(mutgrp)
            selection = self.get_res_selector(position, target_chain)
            if WT_or_MT == 'MT':
                aa_for_swapres = MT
            else:
                aa_for_swapres = WT
            yasara.SwapRes(selection, aa_for_swapres)
            if aa_for_swapres != 'Gly' and aa_for_swapres != 'Ala':
                yasara.OptimizeRes(aa_for_swapres + ' ' + selection, method='SCWALL')
                yasara.Boundary(Type='Periodic')

    def minimize(self):
        yasara.ExperimentMinimization(convergence=0.01)
        yasara.Experiment('On')
        yasara.Wait('ExpEnd')

    def get_component_energies(self, element, element_idx, method='BoundaryFast'):
        # Legacy object-based energy calculation used for split-chain scenes.
        if element in ['Rtr', 'Lgd']:
            yasara.RemoveObj('not ' + str(element_idx))
        yasara.ChargeObj('All', 0)
        epot_list = yasara.Energy('All')
        epot = sum(epot_list)
        if method == 'PBS':
            esolcol = yasara.SolvEnergy(method='PBS')[0]
            _, esolvdw = yasara.SolvEnergy(method='BoundaryFast')
        elif method == 'BoundaryFast':
            esolcol, esolvdw = yasara.SolvEnergy(method='BoundaryFast')
        yasara.Sim('On')
        surfacc_list = yasara.Surf('Accessible')
        surfacc = surfacc_list.pop()
        esurfacc = surfacc * self.surfcost
        esol = esolcol + esolvdw + esurfacc
        yasara.AddObj('All')
        yasara.Sim('Off')
        return epot, esol, esolcol, esolvdw, esurfacc

    def get_selected_component_energies(self, target_chain, method='BoundaryFast'):
        # Evaluate receptor, target ligand, and their complex directly from
        # chain-aware selections inside the loaded scene.
        print(f'Calculating component energies for Mol {target_chain}...')
        receptor_selector = f'Obj 1 and Mol {target_chain}'
        ligand_selector = f'Obj 2 and Mol {target_chain}'
        complex_selector = f'Mol {target_chain}'

        selections = {
            'Cpx': complex_selector,
            'Rtr': receptor_selector,
            'Lgd': ligand_selector,
        }
        object_map = {
            'Rtr': '1',
            'Lgd': '2',
        }
        results = {}

        # Match the previous workflow behavior by neutralizing the loaded scene
        # before deriving per-selection state energies, while keeping the
        # original remove/re-add object workflow for isolated components.
        for suffix, selector in selections.items():
            print(suffix, selector, end=' >> ')
            if suffix in object_map:
                remove_obj_selector = 'not ' + object_map[suffix]
                yasara.RemoveObj(remove_obj_selector)
                print(f'Removed {remove_obj_selector}')
            yasara.ChargeObj('All', 0)
            epot_list = yasara.EnergyMol(selector, component='All')
            epot = sum(epot_list)
            if method == 'PBS':
                esolcol = yasara.SolvEnergyMol(selector, method='PBS')[0]
                solv_energy_list = yasara.SolvEnergyMol(selector, method='BoundaryFast')
                esolvdw = sum([energy for idx, energy in enumerate(solv_energy_list) if idx % 2 == 1])
            elif method == 'BoundaryFast':
                solv_energy_list = yasara.SolvEnergyMol(selector, method='BoundaryFast')
                esolcol = sum([energy for idx, energy in enumerate(solv_energy_list) if idx % 2 == 0])
                esolvdw = sum([energy for idx, energy in enumerate(solv_energy_list) if idx % 2 == 1])
            yasara.Sim('On')
            surfacc_list = yasara.SurfMol(selector, 'Accessible')
            surfacc = surfacc_list.pop()
            esurfacc = surfacc * self.surfcost
            esol = esolcol + esolvdw + esurfacc
            results['epot' + suffix] = epot
            results['esolcol' + suffix] = esolcol
            results['esolvdw' + suffix] = esolvdw
            results['esurfacc' + suffix] = esurfacc
            results['esol' + suffix] = esol
            results['ebind' + suffix] = epot + esol
            yasara.AddObj('All')
            yasara.Sim('Off')

        for feature_prefix in ENERGY_FEATURE_PREFIXES:
            results[feature_prefix] = results[feature_prefix + 'Cpx'] - results[feature_prefix + 'Rtr'] - results[feature_prefix + 'Lgd']
        return results

    def calculate_component_energies(self, target_chain, keep_multiple_chains_in_struct, method='BoundaryFast'):
        if keep_multiple_chains_in_struct:
            print('Energy mode: selection-based')
            return self.get_selected_component_energies(target_chain, method=method)

        print('Energy mode: object-based')
        object_map = {
            'Cpx': None,
            'Rtr': 1,
            'Lgd': 2,
        }
        results = {}
        for element in ENERGY_COMPONENT_SUFFIXES:
            epot, esol, esolcol, esolvdw, esurfacc = self.get_component_energies(
                element,
                object_map[element],
                method=method,
            )
            results['epot' + element] = epot
            results['esol' + element] = esol
            results['esolcol' + element] = esolcol
            results['esolvdw' + element] = esolvdw
            results['esurfacc' + element] = esurfacc
            results['ebind' + element] = epot + esol

        for feature_prefix in ENERGY_FEATURE_PREFIXES:
            results[feature_prefix] = results[feature_prefix + 'Cpx'] - results[feature_prefix + 'Rtr'] - results[feature_prefix + 'Lgd']
        return results

    def process_mutant(
            self,
            mutation,
            struct_name,
            ligand_name,
            target_chain,
            keep_multiple_chains_in_struct,
            move='!backbone',
            minimize_energy=True,
            resetSce=False,
            nrep=5,
            ff='YASARA2',
            mvdist=4,
            mvdrug=1,
            surfout=0.65,
            multiprocessing_proc_num=None,
    ):
        self.surfout = surfout
        self.surfcost = surfout / 6.02214199e20 * self.JToUnit

        # format mutation / mutant
        mutant, mutation = get_mutstr(mutation, sep=self.sep)
        print('mutation:', mutation)
        print('mutant:', mutant)

        # settings and fpaths
        log_fpath = self.get_log_fpath(struct_name, target_chain)
        setname = self.get_setname(ligand_name, struct_name, mutant, ff, resetSce, move, mvdist, mvdrug, self.surfout, self.cntions, nrep)
        seeds = []
        print(multiprocessing_proc_num, 'Processing for ' + struct_name + ' >> ' + mutant + ' on chain ' + str(target_chain))

        # Collect replicate-level energies first, then summarize them
        # into AVG output once all repeats are complete.
        res = {f: [] for f in self.energy_features_list}
        res_avg = {
            'struct': struct_name,
            'target_chain': target_chain,
            'mutations': mutant,
            'ligname': ligand_name,
            'setname': setname,
            'surfcost': self.surfcost,
            'minimize_energy': minimize_energy,
            'resetSce': resetSce,
            'move': move,
            'mvdist': mvdist,
            'mvdrug': mvdrug,
            'ff': ff,
            'counterions': self.cntions,
            'nrep': nrep,
        }

        # INITIALIZE YASARA #
        self.initialize_yasara(ff)

        # allow mutated positions to move
        continue_processing = True
        # Each replicate reuses or reloads the same starting scene
        # depending on resetSce, then evaluates WT and MT states with
        # the same random-seed schedule for reproducibility.
        for n in range(nrep):
            print('rep#=' + str(n))
            setname_n = setname + '-' + str(n)
            # clear and load scene
            seed = 1234*n
            seeds.append(seed)
            yasara.RandomSeed(seed)
            if n == 0 or resetSce:
                # load structure
                self.load_structure(struct_name)
                self.validate_loaded_scene(struct_name, target_chain, keep_multiple_chains_in_struct)
                # set up scene for minimization
                _, continue_processing = self.set_up_sce_for_minimization(mutation, target_chain, move, mvdist, mvdrug, continue_processing)
            if not continue_processing:
                print(f'ERROR: Could not finish processing {struct_name, mutant}. Moving on to the next simulation...')
                break
            else:
                # perform energy minimization and DDG calculations for WT and MT
                for WT_or_MT in ['WT', 'MT']:
                    # perform mutagenesis
                    self.swapres(mutation, WT_or_MT, target_chain)
                    # perform energy minimization
                    if minimize_energy:
                        self.minimize()
                    # get component energies
                    state_energies = self.calculate_component_energies(
                        target_chain,
                        keep_multiple_chains_in_struct,
                        method=self.energy_calc_method,
                    )
                    for element in ENERGY_COMPONENT_SUFFIXES:
                        res['ebind' + WT_or_MT + element].append(state_energies['ebind' + element])
                        res['epot' + WT_or_MT + element].append(state_energies['epot' + element])
                        res['esol' + WT_or_MT + element].append(state_energies['esol' + element])
                        res['esol' + 'col' + WT_or_MT + element].append(state_energies['esolcol' + element])
                        res['esol' + 'vdw' + WT_or_MT + element].append(state_energies['esolvdw' + element])
                        res['esurfacc' + WT_or_MT + element].append(state_energies['esurfacc' + element])

                    # Reconstruct binding free energy from isolated
                    # component terms after the scene has been minimized.
                    for feature_prefix in ENERGY_FEATURE_PREFIXES:
                        res[feature_prefix + WT_or_MT].append(state_energies[feature_prefix])
                    print('Obtained energies for ' + WT_or_MT)
                    print('ebind' + WT_or_MT + ':', res['ebind'+WT_or_MT])
                    print('seeds:', seeds)

                    # save structure as sce
                    setname_prefix, _, setname_suffix = setname_n.partition('_')
                    setname_mod = setname_prefix + '_postOpt-' + WT_or_MT + '_' + setname_suffix
                    output_postOpt_fpath = os.path.join(self.postopt_dir, setname_mod)
                    print(output_postOpt_fpath)
                    yasara.SaveSce(output_postOpt_fpath + '.sce')

                # calculate overall energy change of mutation
                for feature_prefix in ENERGY_FEATURE_PREFIXES:
                    res[feature_prefix + 'DDG'].append(res[feature_prefix + 'MT'][-1] - res[feature_prefix + 'WT'][-1])

        log_fpath_avg, log_fpath_full = None, None
        if continue_processing:
            # Get average and standard deviation for energies calculated
            for f in self.energy_features_list:
                res_avg[f] = round(np.mean(np.asarray(res[f])), 4)
                res_avg[f + '_std'] = round(np.std(np.asarray(res[f])), 4)
            print('ebindDDG=' + str(round(res_avg['ebindDDG'],2)))

            # save average results as csv
            _, log_fpath_avg, write_mode = save_dict_as_csv(res_avg, self.res_avg_cols, log_fpath, multiprocessing_proc_num=multiprocessing_proc_num)
            print('Saved AVG results for ' + struct_name + '_' + str(target_chain) + '_' + mutant + ' to CSV (mode=' + write_mode + ').')

            # save all results as csv
            res.update({'struct': [struct_name]*nrep, 'target_chain': [target_chain]*nrep, 'mutations': [mutant]*nrep, 'ligname': [ligand_name]*nrep, 'setname': [setname+'-'+str(n) for n in range(nrep)], 'n': [n for n in range(nrep)],
                         'surfcost': [self.surfcost]*nrep, 'minimize_energy': [minimize_energy]*nrep, 'resetSce': [resetSce]*nrep, 'move': [move]*nrep,
                         'mvdist': [mvdist]*nrep, 'mvdrug': [mvdrug]*nrep, 'ff': [ff]*nrep, 'counterions': [self.cntions]*nrep, 'nrep':[nrep]*nrep})
            _, log_fpath_full, write_mode = save_dict_as_csv(res, self.res_all_cols, log_fpath, csv_suffix='_FULL', multiprocessing_proc_num=multiprocessing_proc_num)
            print('Saved FULL results for ' + struct_name + '_' + str(target_chain) + '_' + mutant + ' to CSV (mode=' + write_mode + ').')

        return log_fpath_avg, log_fpath_full

    def run_pipeline(self, inputs, output_fname, params, run_multiprocessing=False, save_minimized_struct=True):

        # Parse the batch specification, then expand it into a flat job
        # list across structures, mutations, and parameter settings.
        all_inputs = parse_input_file(inputs, self.input_dir, self.sep)
        scene_jobs = self.prepare_sce_files(all_inputs)

        # get args
        # Prepare argument list
        args_list = [
            (
                mutation,
                scene_job['struct_name'],
                scene_job['ligand_name'],
                scene_job['target_chain'],
                scene_job['keep_multiple_chains_in_struct'],
                move,
                minimize_energy,
                resetSce,
                nrep,
                ff,
                mvdist,
                mvdrug,
                surfout,
            )
            for scene_job in scene_jobs
            for mutation in scene_job['mutations']
            for nrep in params['nrep']
            for ff in params['ff']
            for resetSce in params['resetSce']
            for minimize_energy in params['minimize_energy']
            for move in params['move']
            for mvdist in params['mvdist']
            for mvdrug in params['mvdrug']
            for surfout in params['surfout']
        ]
        print('\n', '[args list]')
        for i, args in enumerate(args_list): print(i, args)
        print()

        # Execute with or without multiprocessing
        if run_multiprocessing:
            # modify args list with proc_num
            args_list = [args + (i,) for i, args in enumerate(args_list)]
            # get log fpath lists
            log_fpath_list = [self.get_log_fpath(args[1], args[3]) + '_' + str(args[-1]) + '.csv' for args in args_list]
            log_fpath_list_FULL = [self.get_log_fpath(args[1], args[3]) + '_FULL_' + str(args[-1]) + '.csv' for args in args_list]
            # perform multiprocessing
            import multiprocessing as mp
            mp.set_start_method("spawn", force=True)
            from multiprocessing import Pool
            with Pool(processes=run_multiprocessing) as pool:
                pool.starmap(parse_args_process_mutant, [(args, self, self.sep) for args in args_list])

        # without multiprocessing
        else:
            log_fpath_list, log_fpath_list_FULL = [], []
            for args in args_list:
                log_fpath_avg, log_fpath_full = parse_args_process_mutant(args, self, self.sep)
                log_fpath_list.append(log_fpath_avg)
                log_fpath_list_FULL.append(log_fpath_full)
            yasara_pids = findProcess('YASARA.exe' if opsys == 'Windows' else 'yasara')
            if self.exit_program and len(yasara_pids) > 0:
                yasara_pid = yasara_pids[-1]
                exit_program(yasara_pid)

        # combine files spawned
        ## AVG results
        print('log_fpath_list:', log_fpath_list)
        log_fpath_list = [f for f in log_fpath_list if f is not None]
        missing_data = combine_csv_files(list(set(log_fpath_list)), self.output_dir, output_fname, remove_combined_files=True, append_to_existing_output=self.append_to_existing_output)
        print('Saved all AVG results as CSV:' + os.path.join(self.output_dir, output_fname + '.csv'))
        if len(missing_data)>0: print('Missing:', missing_data)

        ## FULL results
        log_fpath_list_FULL = [f for f in log_fpath_list_FULL if f is not None]
        missing_data = combine_csv_files(list(set(log_fpath_list_FULL)), self.output_dir, output_fname + '_FULL')
        print('Saved all FULL results as CSV:' + os.path.join(self.output_dir, output_fname + '_FULL.csv'))
        if len(missing_data)>0: print('Missing:', missing_data)

        if save_minimized_struct:
            self.save_sce_as_pdb(self.output_dir, struct_subdir=POSTOPT_SUBDIR, del_sce_files=True)
        combined_failed_jobs = self.combine_failed_job_logs(remove_combined_files=False)
        if combined_failed_jobs is not None:
            print('Saved failed job log as CSV:' + combined_failed_jobs)


if __name__ == "__main__":

    data_config = DATA_CONFIG.copy()
    run_config = RUN_CONFIG.copy()
    data_folder = data_config['data_folder']
    data_subfolder = data_config['data_subfolder']
    output_fname = data_config['output_fname']
    inputs = data_config['inputs']
    run_multiprocessing = run_config['run_multiprocessing']
    save_minimized_struct = run_config['save_minimized_struct']
    fix_metal_ion = run_config['fix_metal_ion']
    append_to_existing_output = run_config['append_to_existing_output']
    sep = run_config['sep']
    energy_calc_method = run_config['energy_calc_method']
    params = OrderedDict((key, list(value)) for key, value in RUN_PARAMS.items())

    # initialize MutBindDDG object
    mutbindDDG = MutBindDDG(data_folder, data_subfolder, fix_metal_ion=fix_metal_ion,
                            append_to_existing_output=append_to_existing_output, sep=sep, energy_calc_method=energy_calc_method)

    # run pipeline
    mutbindDDG.run_pipeline(inputs, output_fname, params, run_multiprocessing, save_minimized_struct)
