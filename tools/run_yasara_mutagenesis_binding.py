# -*- coding: utf-8 -*-
import yasara
import numpy as np
import os
from collections import OrderedDict
from variables import address_dict, subfolders, aaList
from utils import opsys, get_mutstr, pdb_to_sce, split_mutation, findProcess, exit_program, save_dict_as_csv, combine_csv_files
if opsys == 'Windows':
    process_name = 'YASARA.exe'
else:
    process_name = 'yasara'

# CALCULATE BINDING ENERGY INCLUDING IMPLICIT SOLVATION
# Created by Sebastian Maurer-Stroh, inspired by template macros from Elmar Krieger.
# Edited by Jhoann Miyajima and Charmaine Chia

def parse_input_file(inputs, input_dir, sep):
    """Get struct name, ligand, and mutations to process from inputs file"""
    import pandas as pd
    from utils import mutstr_to_mutations

    if isinstance(inputs,dict):
        struct_to_mutate_dict = inputs

    elif isinstance(inputs,str) and inputs.find('.csv')>-1:
        inputs_df = pd.read_csv(input_dir+inputs)
        print(inputs_df)

        # Group rows by unique structure / ligand combinations so each
        # prepared scene can be paired with its mutation list.
        struct_to_mutate_dict = {}
        struct_ligand_df = inputs_df[['struct_name', 'ligand_name']].drop_duplicates(ignore_index=True)
        struct_ligand_pairs = [(struct_name, ligand_name) for struct_name, ligand_name in zip(struct_ligand_df['struct_name'].tolist(), struct_ligand_df['ligand_name'].tolist())]

        for (struct_name, ligand_name) in struct_ligand_pairs:

            # get filtered df
            inputs_df_filt = inputs_df.loc[
                          (inputs_df['struct_name'] == struct_name) &
                          (inputs_df['ligand_name'] == ligand_name) &
                          (inputs_df['process_structure'] == 1), :
                          ]
            if len(inputs_df_filt)==0:
                print((struct_name, ligand_name), '>> No entries after filtering.')
                continue
            ligand_id, chain_id, ligand_chain_id = None, None, None
            # get ligand ID, if provided
            if 'ligand_id' in inputs_df_filt:
                ligand_id = inputs_df_filt.iloc[0]['ligand_id']
            # get chain ID, if provided
            if 'chain_id' in inputs_df_filt:
                chain_id = inputs_df_filt.iloc[0]['chain_id']
            # get chain ID, if provided
            if 'ligand_chain_id' in inputs_df_filt:
                ligand_chain_id = inputs_df_filt.iloc[0]['ligand_chain_id']

            # Either take explicit mutations from the CSV or enumerate
            # all allowed substitutions from a supplied sequence.
            if 'mutations' in inputs_df_filt:
                mutations_ = inputs_df_filt['mutations'].tolist()
                mutations = []
                # parse individual entries if they contain a list of separate mutations
                for mut_entry in mutations_:
                    mut_list = mut_entry.split(',') # split separate mutants
                    # use '+' as a separator for combi mutations
                    for k, mut in enumerate(mut_list):
                        if not mut.isalnum():
                            combi_mut, sep_orig = mutstr_to_mutations(mut,sep)
                            mut_list[k] = combi_mut
                    mutations += mut_list
            # no specific mutations provided
            elif 'mutations' not in inputs_df_filt and 'sequence' in inputs_df_filt:
                from utils import list_all_mutations
                res_to_mutate = None
                wt_seq = inputs_df_filt.iloc[0]['sequence']
                # get specific residues to process from inputs_df
                if 'res' in inputs_df:
                    res_to_mutate_ = inputs_df_filt['res'].tolist()
                    res_to_mutate = []
                    for res_entry in res_to_mutate_:
                        res_list = [int(res) for res in res_entry.split(',')]
                        res_to_mutate += res_list
                    res_to_mutate = list(set(res_to_mutate))
                    res_to_mutate.sort()
                # elucidate mutations
                mutations = list_all_mutations(wt_seq, ignore_mutations_to_WT=True, pos_to_mutate=res_to_mutate)
            else:
                print('"mutations", "res", and "sequence" columns not provided in input CSV. Unable to determine mutations for', (struct_name, ligand_name))
                break

            # update struct_to_mutate_dict entry
            struct_to_mutate_dict.update({
                struct_name: {
                    'ligand_name':ligand_name,
                    'ligand_id':ligand_id,
                    'chain_id':chain_id,
                    'ligand_chain_id': ligand_chain_id,
                    'mutations':mutations
                }
            })

    # print out inputs
    for struct_name, struct_dict in struct_to_mutate_dict.items():
        print(struct_name)
        print(struct_dict)

    return struct_to_mutate_dict

def parse_args_process_mutant(args, mutbindDDG, sep):
    if len(args)==12:
        mutation, struct_name, ligand_name, move, minimize_energy, resetSce, nrep, save_minimized_struct, ff, mvdist, mvdrug, surfout = args
        proc_num = None
    elif len(args)==13:
        mutation, struct_name, ligand_name, move, minimize_energy, resetSce, nrep, save_minimized_struct, ff, mvdist, mvdrug, surfout, proc_num = args
    mutstr, mutation = get_mutstr(mutation, sep=sep)
    print(proc_num, struct_name + " Mutation: " + mutstr + "; Move: " + move + "; Minimize: " + str(minimize_energy))
    log_fpath_avg, log_fpath_full = mutbindDDG.process_mutant(mutation, struct_name, ligand_name, move, minimize_energy, resetSce, nrep, save_minimized_struct, ff, mvdist, mvdrug, surfout, proc_num)
    return log_fpath_avg, log_fpath_full

class MutBindDDG:
    def __init__(
            self,
            data_folder,
            data_subfolder,
            struct_format='sce',
            input_subfolder='yasara/Input/',
            output_subfolder='yasara/Output/',
            keep_metal_ion=True,
            fix_metal_ion=True,
            skip_add_hyd=False,
            cntions=0, # If counterions=1, counter ions will be implicitly considered by setting the net charges to 0
            JToUnit=1.43932620865201e20,
            structure_dict={'Cpx': 3, 'Rtr': 1, 'Lgd': 2},
            sep='-',
            exit_program_after_processing=True,
            append_to_existing_output=True,
            energy_calc_method='BoundaryFast'
    ):
        self.data_folder = data_folder
        self.data_subfolder = data_subfolder
        self.struct_format = struct_format
        self.struct_subfolder = struct_format + '/'
        self.struct_dir = data_folder + self.struct_subfolder + data_subfolder + '/'
        self.input_dir = data_folder + input_subfolder
        self.output_dir = data_folder + output_subfolder + data_subfolder + '/'
        self.keep_metal_ion = keep_metal_ion
        self.fix_metal_ion = fix_metal_ion
        self.skip_add_hyd = skip_add_hyd
        self.cntions = cntions
        self.structure_dict = structure_dict
        self.JToUnit = JToUnit
        self.sep = sep
        self.exit_program=exit_program_after_processing
        self.append_to_existing_output = append_to_existing_output
        self.energy_calc_method = energy_calc_method

        # results columns
        metadata_list = ['struct', 'mutations', 'ligname', 'setname', 'minimize_energy', 'resetSce', 'move', 'mvdist', 'mvdrug', 'ff', 'surfcost', 'counterions', 'nrep']
        self.energy_features_list = [feature + suffix for feature in ['ebind', 'epot', 'esol', 'esolcol', 'esolvdw', 'esurfacc'] for suffix in ['DDG', 'MT', 'WT'] + [WT_or_MT + el for WT_or_MT in ['MT', 'WT'] for el in ['Cpx', 'Rtr', 'Lgd']]]
        self.res_avg_cols = metadata_list + [f+suffix for f in self.energy_features_list for suffix in ['', '_std']]
        self.res_all_cols = metadata_list + ['n'] + self.energy_features_list

        # create output directory if it doesn't exist
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print('Created output directory:', self.output_dir)
        if not os.path.exists(self.output_dir + 'postOpt/'):
            os.makedirs(self.output_dir + 'postOpt/')
            print('Created postOpt sub-directory:', self.output_dir + 'postOpt/')

        # initialize yasara
        self.initialize_yasara()

    def prepare_sce_files(self, all_inputs):
        struct_to_mutations = {}
        # When starting from PDB input, split receptor / ligand content
        # into YASARA scene files first, then map each generated scene
        # back to the requested mutation set.
        if self.struct_format == 'pdb':
            print('Processing PDB to generate SCE files...')

            for struct_name, struct_inputs_dict in all_inputs.items():
                ligand_name = struct_inputs_dict['ligand_name']
                ligand_id = struct_inputs_dict['ligand_id']
                chain_id = struct_inputs_dict['chain_id']
                ligand_chain_id = struct_inputs_dict['ligand_chain_id']
                if chain_id is not None:
                    chain_id = chain_id.split(',')
                mutations = struct_inputs_dict['mutations']
                pdb_fpath = self.data_folder + subfolders['pdb'] + self.data_subfolder + '/' + struct_name + '.pdb'
                sce_fpaths_bystruct = pdb_to_sce(pdb_fpath, ligand_id, chains_to_process=chain_id, ligand_chain_id=ligand_chain_id, keep_ligand=True, keep_metal_ion=self.keep_metal_ion)
                # update struct_to_mutate_dict_updated
                struct_to_mutations.update({
                    (os.path.basename(sce_fpath).replace('.sce', ''), ligand_name): mutations
                    for sce_fpath in sce_fpaths_bystruct})
        else:
            for struct_name, struct_inputs_dict in all_inputs.items():
                ligand_name = struct_inputs_dict['ligand_name']
                mutations = struct_inputs_dict['mutations']
                struct_to_mutations.update({(struct_name, ligand_name): mutations})
        return struct_to_mutations

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
        sce_list = [f for f in os.listdir(output_dir + struct_subdir) if f.endswith('.sce')]
        for sce_fname in sce_list:
            yasara.Clear()
            yasara.LoadSce(output_dir + struct_subdir + sce_fname)
            yasara.JoinObj('2', '1', center='No')
            yasara.NameMol('Protein', 'R')
            yasara.NameMol('not Protein', 'L')
            yasara.SavePDB('1', output_dir + struct_subdir + sce_fname.replace('.sce', '.pdb'))
            print('Converted sce to pdb file:', sce_fname)
        if del_sce_files:
            for sce_fname in sce_list:
                os.remove(output_dir+struct_subdir+sce_fname)

    def get_setname(self, ligname, struct, mutant, ff, resetSce, move, mvdist, mvdrug, surfout, cntions, nrep):
        return '{0}_{1}_{2}_{3}_rs{4}_{5}_d{6}_md{7}_s{8}_c{9}_r{10}'.format(ligname, struct, mutant, ff,
                                                                                resetSce * 1, move, mvdist, mvdrug,
                                                                                surfout, cntions, nrep)
    def load_structure(self, struct_name):
        yasara.Clear()
        yasara.LoadSce(self.struct_dir.replace('pdb/', 'sce/') + struct_name + '.sce')
        # Standardize protonation / geometry before mutagenesis so each
        # replicate starts from the same cleaned scene state.
        if self.skip_add_hyd:
            yasara.CleanAll(skip='AddHyd')
        else:
            yasara.CleanAll()
        yasara.CellAuto(extension='10')
        yasara.FixAll()

    def set_up_sce_for_minimization(self, mutation, move, mvdist, mvdrug, continue_processing=True):
        # Free only the local neighborhood around the mutation site
        # unless the run configuration explicitly allows broader motion.
        mutname = ''
        for mutgrp in mutation:
            WT, position, MT = split_mutation(mutgrp)
            # check if current residue position is the same
            WT_actual = yasara.NameRes(position)[0]
            print('Actual WT residue:', WT_actual, '; Target WT residue:', WT.upper())
            if WT_actual != WT.upper():
                print('Incorrect WT amino acid found at position to mutate >> End processing')
                continue_processing = False
                break
            else:
                print('Correct WT amino acid found at position to mutate >> Continue processing')
            x = '{0}{1}{2}'.format(WT, position, MT)
            mutname += x
            yasara.ShowRes(str(position))
            yasara.FreeAtom(move + ' with distance <' + str(mvdist) +' from res ' + str(position))
        print('mutname:', mutname)

        if (mvdrug == 1): # allow ligand and close by residues to move
                yasara.FreeAtom(move + ' with distance <' + str(mvdist) +' from Obj 2')
        if (mvdrug == 0): # fix ligand
                yasara.FixObj("2")
        # fix metal ion
        if self.fix_metal_ion:
            yasara.FixAtom('metal')
        return mutname, continue_processing

    def swapres(self, mutation, WT_or_MT):
        # Build either the WT or MT state directly inside the current
        # YASARA scene before minimization / energy evaluation.
        for mutgrp in mutation:
            print(mutgrp)
            WT, position, MT = split_mutation(mutgrp)
            selection = str(position)
            if WT_or_MT == 'MT':
                aa_for_swapres = MT
            elif WT_or_MT == 'WT':
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
        # Evaluate each component in isolation so binding energy can be
        # reconstructed from complex - receptor - ligand terms.
        if element in ['Rtr', 'Lgd']:
            yasara.RemoveObj('not ' + str(element_idx))
        # set charge to 0
        yasara.ChargeObj('All', 0)
        # get potential energy
        epotList = yasara.Energy('All')
        epot = sum(epotList)
        # get component solvation energies
        if method=='PBS':
            esolcol = yasara.SolvEnergy(method='PBS')[0]
            _, esolvdw = yasara.SolvEnergy(method='BoundaryFast')
        elif method=='BoundaryFast':
            esolcol, esolvdw = yasara.SolvEnergy(method='BoundaryFast') # BoundaryFast method
        yasara.Sim('On')
        # get interfacial energy
        surfaccList = yasara.Surf('Accessible')
        surfacc = surfaccList.pop()
        esurfacc = surfacc * self.surfcost
        # get total solvation energy
        esol = esolcol + esolvdw + esurfacc
        yasara.AddObj('All')
        return epot, esol, esolcol, esolvdw, esurfacc

    def process_mutant(
            self,
            mutation,
            struct_name,
            ligand_name,
            move='!backbone',
            minimize_energy=True,
            resetSce=False,
            nrep=5,
            save_minimized_struct=False,
            ff = 'YASARA2',
            mvdist = 4,
            mvdrug = 1,
            surfout = 0.65,
            multiprocessing_proc_num=None,
            del_sce_input_after_processing=False
    ):
        self.surfout = surfout
        self.surfcost = surfout / 6.02214199e20 * self.JToUnit

        # format mutation / mutant
        mutant, mutation = get_mutstr(mutation, sep=self.sep)
        print('mutation:', mutation)
        print('mutant:', mutant)

        # settings and fpaths
        log_fname = 'DDG_' + struct_name
        log_fpath = self.output_dir + log_fname
        setname = self.get_setname(ligand_name, struct_name, mutant, ff, resetSce, move, mvdist, mvdrug, self.surfout, self.cntions, nrep)
        seeds = []
        print(multiprocessing_proc_num, 'Processing for ' + struct_name + ' >> ' + mutant)

        # Collect replicate-level energies first, then summarize them
        # into AVG output once all repeats are complete.
        res = {f: [] for f in self.energy_features_list}
        res_avg = {'struct':struct_name, 'mutations':mutant, 'ligname':ligand_name, 'setname':setname,
                   'surfcost':self.surfcost, 'minimize_energy': minimize_energy, 'resetSce':resetSce, 'move':move, 'mvdist':mvdist, 'mvdrug':mvdrug, 'ff':ff, 'counterions':self.cntions, 'nrep':nrep}

        # INITIALIZE YASARA #
        self.initialize_yasara(ff)

        # allow mutated positions to move
        continue_processing = True
        # Each replicate reuses or reloads the same starting scene
        # depending on resetSce, then evaluates WT and MT states with
        # the same random-seed schedule for reproducibility.
        for n in range(nrep):
            print('rep#='+str(n))
            setname_n = setname + '-' + str(n)
            # clear and load scene
            seed = 1234*n
            seeds.append(seed)
            yasara.RandomSeed(seed)
            if n==0 or resetSce:
                # load structure
                self.load_structure(struct_name)
                # set up scene for minimization
                mutname, continue_processing = self.set_up_sce_for_minimization(mutation, move, mvdist, mvdrug, continue_processing)
            if not continue_processing:
                print(f'ERROR: Could not finish processing {struct_name, mutant}. Moving on to the next simulation...')
                break
            else:
                # perform energy minimization and DDG calculations for WT and MT
                for WT_or_MT in ['WT', 'MT']:
                    # perform mutagenesis
                    self.swapres(mutation, WT_or_MT)
                    # perform energy minimization
                    if minimize_energy:
                        self.minimize()
                    # get component energies
                    for element, element_idx in self.structure_dict.items():
                        epot, esol, esolcol, esolvdw, esurfacc = self.get_component_energies(element, element_idx, method=self.energy_calc_method)
                        res['ebind' + WT_or_MT + element].append(epot + esol)
                        res['epot' + WT_or_MT + element].append(epot)
                        res['esol' + WT_or_MT + element].append(esol)
                        res['esol' + 'col' + WT_or_MT + element].append(esolcol)
                        res['esol' + 'vdw' + WT_or_MT + element].append(esolvdw)
                        res['esurfacc' + WT_or_MT + element].append(esurfacc)

                    # Reconstruct binding free energy from isolated
                    # component terms after the scene has been minimized.
                    for feature_prefix in ['ebind', 'epot', 'esol', 'esolcol', 'esolvdw', 'esurfacc']:
                        res[feature_prefix + WT_or_MT].append(res[feature_prefix + WT_or_MT + 'Cpx'][-1] - res[feature_prefix + WT_or_MT + 'Rtr'][-1] - res[feature_prefix + WT_or_MT + 'Lgd'][-1])
                    yasara.Sim('Off')
                    print('Obtained energies for ' + WT_or_MT)
                    print('ebind' + WT_or_MT + ':', res['ebind'+WT_or_MT])
                    print('seeds:', seeds)

                    # save structure as sce
                    setname_mod = setname_n[:setname_n.find('_')] + '_postOpt-'+WT_or_MT + setname_n[setname_n.find('_'):]
                    output_postOpt_fpath = self.output_dir + 'postOpt/' + setname_mod
                    print(output_postOpt_fpath)
                    yasara.SaveSce(output_postOpt_fpath + '.sce')

                # calculate overall energy change of mutation
                for feature_prefix in ['ebind', 'epot', 'esol', 'esolcol', 'esolvdw', 'esurfacc']:
                    res[feature_prefix + 'DDG'].append(res[feature_prefix + 'MT'][-1] - res[feature_prefix + 'WT'][-1])
                yasara.Sim('Off')

        log_fpath_avg, log_fpath_full = None, None
        if continue_processing:
            # Get average and standard deviation for energies calculated
            for f in self.energy_features_list:
                res_avg[f] = round(np.mean(np.array(res[f])),4)
                res_avg[f+'_std'] = round(np.std(np.array(res[f])),4)
            print('ebindDDG=' + str(round(res_avg['ebindDDG'],2)))

            # save average results as csv
            csv_txt_avg, log_fpath_avg, write_mode = save_dict_as_csv(res_avg, self.res_avg_cols, log_fpath, multiprocessing_proc_num=multiprocessing_proc_num)
            print('Saved AVG results for ' + struct_name + mutant + ' to CSV (mode=' + write_mode + ').')

            # save all results as csv
            res.update({'struct': [struct_name]*nrep, 'mutations': [mutant]*nrep, 'ligname': [ligand_name]*nrep, 'setname': [setname+'-'+str(n) for n in range(nrep)], 'n': [n for n in range(nrep)],
                         'surfcost': [self.surfcost]*nrep, 'minimize_energy': [minimize_energy]*nrep, 'resetSce': [resetSce]*nrep, 'move': [move]*nrep,
                         'mvdist': [mvdist]*nrep, 'mvdrug': [mvdrug]*nrep, 'ff': [ff]*nrep, 'counterions': [self.cntions]*nrep, 'nrep':[nrep]*nrep})
            csv_txt_avg, log_fpath_full, write_mode = save_dict_as_csv(res, self.res_all_cols, log_fpath, csv_suffix='_FULL', multiprocessing_proc_num=multiprocessing_proc_num)
            print('Saved FULL results for ' + struct_name + mutant + ' to CSV (mode=' + write_mode + ').')

        if save_minimized_struct:
            # convert sce to pdb files
            self.save_sce_as_pdb(self.output_dir, struct_subdir='postOpt/', del_sce_files=True)

        return log_fpath_avg, log_fpath_full

    def run_pipeline(self, inputs, output_fname, params, run_multiprocessing=False, save_minimized_struct=True):

        # Parse the batch specification, then expand it into a flat job
        # list across structures, mutations, and parameter settings.
        all_inputs = parse_input_file(inputs, self.input_dir, self.sep)
        struct_to_mutations = self.prepare_sce_files(all_inputs)

        # get args
        # Prepare argument list
        args_list = [
            (mutation, struct_name, ligand_name, move, minimize_energy, resetSce, nrep, save_minimized_struct, ff, mvdist, mvdrug, surfout)
            for (struct_name, ligand_name), mutations in struct_to_mutations.items()
            for mutation in mutations
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
            log_fpath_list = [self.output_dir + 'DDG_' + args[1] + '_' + str(args[-1]) + '.csv' for args in args_list]
            log_fpath_list_FULL = [self.output_dir + 'DDG_' + args[1] + '_FULL_' + str(args[-1]) + '.csv' for args in args_list]
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
            yasara_pid = findProcess('YASARA.exe' if opsys == 'Windows' else 'yasara')[-1]
            if self.exit_program:
                exit_program(yasara_pid)

        # combine files spawned
        ## AVG results
        print('log_fpath_list:', log_fpath_list)
        log_fpath_list = [f for f in log_fpath_list if log_fpath_list is not None]
        missing_data = combine_csv_files(list(set(log_fpath_list)), self.output_dir, output_fname, remove_combined_files=True, append_to_existing_output=self.append_to_existing_output)
        print('Saved all AVG results as CSV:' + self.output_dir + output_fname + '.csv')
        if len(missing_data)>0: print('Missing:', missing_data)

        ## FULL results
        log_fpath_list_FULL = [f for f in log_fpath_list_FULL if log_fpath_list_FULL is not None]
        missing_data = combine_csv_files(list(set(log_fpath_list_FULL)), self.output_dir, output_fname + '_FULL')
        print('Saved all FULL results as CSV:' + self.output_dir + output_fname + '_FULL' + '.csv')
        if len(missing_data)>0: print('Missing:', missing_data)


if __name__ == "__main__":

    # set inputs
    data_folder = address_dict['influenza-resistance'] #
    data_subfolder = 'PA-A37T-I38T' # 'PA_benchmark_mod' # 'IC50benchmark' # 'platinum_influenzaonly'
    output_fname = 'PA-A37T-I38T' # 'PA_benchmark_mod' # 'DDG_IC50benchmark' # 'DDG_platinum_influenzaonly'
    struct_format = 'sce' # 'pdb' #
    inputs = 'PA-A37T-I38T.csv' # 'PA_benchmark_mod.csv' #
    run_multiprocessing = None # 6 # 16 # #
    save_minimized_struct = True # False
    keep_metal_ion = True
    fix_metal_ion = True # False #
    skip_add_hyd = False # True
    append_to_existing_output = True
    sep = '+'
    energy_calc_method = 'BoundaryFast' # 'PBS' #

    # set params
    params = OrderedDict()
    params['nrep'] = [10] # [5] #
    params['minimize_energy'] = [True]
    params['resetSce'] = [False] #
    params['move'] = ['!backbone'] # ['!backbone', 'all'] #
    params['ff'] = ['YASARA2'] # ['AMBER14'] # ['YASARA2', 'AMBER14', 'AMBER15FB']
    params['mvdist'] = [4] # [4,5] #
    params['mvdrug'] = [1]
    params['surfout'] = [0.65]

    # initialize MutBindDDG object
    mutbindDDG = MutBindDDG(data_folder, data_subfolder, struct_format=struct_format, keep_metal_ion=keep_metal_ion, fix_metal_ion=fix_metal_ion, skip_add_hyd=skip_add_hyd,
                            append_to_existing_output=append_to_existing_output, sep=sep, energy_calc_method=energy_calc_method)

    # run pipeline
    mutbindDDG.run_pipeline(inputs, output_fname, params, run_multiprocessing, save_minimized_struct)
