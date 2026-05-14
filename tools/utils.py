#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024

@author: charmainechia
"""
import pandas as pd
import os
import platform
from variables import mapping, aaList

opsys = platform.system()


def mkDir(res, output_dir, remove_existing_dir=True):
    import shutil
    # making new directory
    new_dir = (output_dir + res)
    if os.path.exists(new_dir):
        # remove if directory exists, and make new directory
        if remove_existing_dir:
            shutil.rmtree(new_dir)
            os.makedirs(new_dir)
    else:
        os.makedirs(new_dir)
    return new_dir


def get_mapping(mapping, res):
    # Map residue initial to code
    result = mapping[res]
    return result


def mutstr_to_mutations(mutstr, sep=None):
    """Convert mutation string to list of mutations"""
    mutations = []
    # check if it is a single or combi mutation
    ## single mutation identified
    if mutstr.isalnum():
        mutations.append(mutstr)
        return mutations
    ## combi mutation identified
    else:
        if sep is None:
            sep = list(set([char for char in list(mutstr) if not char.isalnum()]))[0]
        # get mutations by splitting mutstr at separator
        mutations = mutstr.split(sep)
        return mutations, sep
    

def split_mutation(mutation, aa_letter_representation=False):
    # Convert point mutation to wildtype residue, muted residue and mutation position
    mutation = list(mutation)
    WT_res = mutation[0]
    MUT_res = mutation[-1]
    if not aa_letter_representation:
        WT_res = get_mapping(mapping, WT_res)
        MUT_res = get_mapping(mapping, MUT_res)
    MUT_pos = mutation[1:len(mutation)-1]
    MUT_pos = int(''.join(MUT_pos))
    return WT_res, MUT_pos, MUT_res


def list_all_mutations(seq, ignore_mutations_to_WT=True, pos_to_mutate=None):
    if pos_to_mutate is None:
        pos_to_mutate = [i+1 for i in range(len(seq)) if seq[i]!='-']
    mut_all = []
    for pos in pos_to_mutate:
        i = pos-1
        wt_aa = seq[i]
        if ignore_mutations_to_WT:
            mut_pos = [wt_aa + str(pos) + aa for aa in aaList if aa!=wt_aa]
        else:
            mut_pos = [wt_aa + str(pos) + aa for aa in aaList]
        mut_all += mut_pos
    return mut_all


def get_mutstr(mutation, sep='+'):
    """Convert input mutation or list of mutations to mutation string and mutation list"""
    # multiple mutations
    if isinstance(mutation, list):
        mutstr = sep.join(mutation)
    # single mutation
    else:
        mutstr = mutation
        mutation = mutation.split(sep)
    return mutstr, mutation


def pdb_to_sce(pdb_fpath, ligand_name, chains_to_process=None, ligand_chain_id=None, keep_ligand=True, keep_metal_ion=True, split_lig_obj=True, save_sce=True, skip_chains=[]):
    import string
    import yasara
    from variables import nonstandard_amino_acids, hetatm_non_metal_ion, hetatm_metal_ion

    # load pdb
    yasara.Clear()
    yasara.LoadPDB(pdb_fpath)
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
    if len(chain_list) > 1 and len(non_receptor_chains) > 0 and len(receptor_chains)==1:
        # rename all chains with the chain_id of the receptor chain
        yasara.NameMol('all', receptor_chains[0])
        # update chain list
        chain_list = list(set(yasara.NameMol('Obj 1')))
        chain_list.sort()
        print('Chain list (after merging ligand with receptor chain):', len(chain_list), chain_list)

    if chains_to_process is None:
        chains_to_process = chain_list.copy()
    else:
        if isinstance(chains_to_process, str):
            chains_to_process = list(chains_to_process)

    # remove specific chains, if specified
    chains_to_process = [c for c in chains_to_process if c not in skip_chains]
    print('Final Chains to Process:', chains_to_process)

    # retain only chains of interest
    sce_fpath_list = []
    for chain in chains_to_process:
        print('Processing Chain', chain, '...')
        # reload file
        yasara.Clear()
        yasara.LoadPDB(pdb_fpath)

        # get non-protein residues
        receptor_chains_to_delete = [c for c in chain_list if c not in chains_to_process+non_receptor_chains]
        yasara.DelMol(f'not {chain}')

        # if ligands are named on separate non-alphanumeric chains
        # yasara.DelMol(f'not {chain} {" and not ".join(non_receptor_chains)}')
        # del_cmd = f'{" and ".join(receptor_chains_to_delete)}'
        # yasara.DelMol(del_cmd)
        # print('del_cmd:', del_cmd)
        chain_list_byreceptor = list(set(yasara.NameMol('Obj 1')))
        print(chain_list_byreceptor)
        resname_notprotein = [r for r in yasara.NameRes('not Protein') if r not in nonstandard_amino_acids]
        print('Non-protein residues (before deletion):', resname_notprotein)

        # get ligand name
        if ligand_name is None:
            ligand_name = [r for r in resname_notprotein if r not in  hetatm_non_metal_ion + hetatm_metal_ion][0]
        print('Ligand ID:', ligand_name)

        # retain only Protein, Ligand, and optionally metal ions.
        # Copy the imported list so we do not mutate module-level constants
        # across repeated calls to pdb_to_sce.
        hetatm_res_to_keep = list(hetatm_non_metal_ion)
        if keep_metal_ion:
            hetatm_res_to_keep += hetatm_metal_ion
        delres_cmd = "not " + " and not ".join(["Protein"] + hetatm_res_to_keep + nonstandard_amino_acids)
        if keep_ligand:
            delres_cmd += ' and not ' + ligand_name

        # delete unnecessary components
        # print('Deletion cmd:', delres_cmd)
        yasara.DelRes(delres_cmd)
        resname_notprotein = [r for r in yasara.NameRes('not Protein') if r not in nonstandard_amino_acids]
        print('Non-Protein residues (after deletion)', len(resname_notprotein), resname_notprotein)

        # split out ligand into 2nd object
        if split_lig_obj:
            print('# of obj (before splitting Obj 1):', yasara.CountObj('all'))
            yasara.SplitObj(1)
            print('# of obj (after splitting Obj 1):', yasara.CountObj('all'))
            yasara.JoinObj(f'not Res {ligand_name}', 1)
            print('# of obj (after joining non-ligand objects to Obj 1):', yasara.CountObj('all'))
            obj_list = yasara.ListObj('all')
            obj_list.sort()
            yasara.SwapObj(obj_list[-1], 2)
            obj_list = yasara.ListObj('all')
            print(obj_list)

        if save_sce:
            if len(chain_list) == 1:
                chain_suffix = ''
            else:
                chain_suffix = '-' + chain
            sce_dir = os.path.dirname(pdb_fpath).replace('pdb/', 'sce/')
            sce_fname = os.path.basename(pdb_fpath).replace('.pdb', chain_suffix + '.sce')
            sce_fpath = sce_dir + '/' + sce_fname
            if not os.path.exists(sce_dir):
                os.makedirs(sce_dir)
                print('Created postOpt sub-directory:', sce_dir)
            yasara.SaveSce(sce_fpath)
            print('Saved .sce file:', sce_fpath)
            sce_fpath_list.append(sce_fpath)

    return sce_fpath_list


def findProcess(process_name):
    if opsys=='Windows':
        return [int(item.split()[1]) for item in os.popen('tasklist').read().splitlines()[4:] if process_name in item.split()]
    elif opsys=='Linux' or opsys=='Darwin':
        return [int(pid) for pid in os.popen('pidof '+process_name).read().strip(' \n').split(' ') if pid!='']
    

def exit_program(pid):
    import signal
    print("Sending SIGINT to self...")
    os.kill(pid, signal.SIGINT)
    print('Exited program', pid)
def is_float(string):
    if string.replace(".", "").replace("-", "").isnumeric():
        return True
    else:
        return False
    
    
def save_dict_as_csv(datadict, cols, log_fpath, csv_suffix ='', multiprocessing_proc_num=None):
    # save results as CSV
    csv_txt = ''
    # get csv_suffix if running multiprocessing
    if multiprocessing_proc_num is not None:
        csv_suffix += '_' + str(multiprocessing_proc_num)

    # check if file exists yet
    log_fpath_full = log_fpath + csv_suffix + '.csv'
    if not os.path.exists(log_fpath_full):
        # if not, start a new file with headers
        write_mode = 'w'
        csv_txt += ','.join(cols) + '\n'
    else:
        write_mode = 'a'

    # convert dict of lists to list of dicts
    if isinstance(datadict[cols[0]], list):
        num_rows = len(datadict[cols[0]])
        datadict_byrow = []
        for row_idx in range(num_rows):
            row = []
            for col in cols:
                row.append(datadict[col][row_idx])
            datadict_byrow.append(row)
    else:
        row = []
        for col in cols:
            row.append(datadict[col])
        datadict_byrow = [row]

    # add data to csv file
    for row in datadict_byrow:
        csv_txt += ','.join([str(el) for el in row])
        csv_txt += '\n'
    # save the changes
    with open(log_fpath_full, write_mode) as f:
        f.write(csv_txt)
    return csv_txt, log_fpath_full, write_mode

def combine_csv_files(log_fpath_list, output_dir, output_fname, remove_combined_files=True, append_to_existing_output=True):
    # combine files spawned
    missing_data = []
    combined_data = []
    df_all = None
    # fetch logged result
    for i, log_fpath in enumerate(log_fpath_list):
        if os.path.exists(log_fpath):
            df = pd.read_csv(log_fpath, index_col=0)
            if len(combined_data)==0:
                df_all = df.copy()
            else:
                df_all = pd.concat([df_all, df], axis=0)
            combined_data.append(log_fpath)
            print(i, log_fpath, '>> combined')
        else:
            missing_data.append(log_fpath)
            print(i, log_fpath, '>> MISSING')

    # update combined results
    output_csv_fpath = output_dir + output_fname + '.csv'
    if df_all is not None:
        if os.path.exists(output_csv_fpath) and append_to_existing_output:
            df_existing = pd.read_csv(output_csv_fpath, index_col=0)
            df_all = pd.concat([df_existing, df_all], axis=0)
        df_all.to_csv(output_csv_fpath)

    # record missing files
    if len(missing_data)>0:
        with open(output_dir + 'missing_data.txt', 'w') as f:
            missing_txt = '\n'.join(missing_data) + '\n'
            f.write(missing_txt)

    # remove combined files
    if remove_combined_files:
        for log_fpath in [f for f in log_fpath_list if f not in missing_data]:
            os.remove(log_fpath)
    return missing_data
