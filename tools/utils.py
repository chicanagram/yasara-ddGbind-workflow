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


def drop_unnamed_columns(df):
    """Remove CSV index columns accidentally persisted by pandas."""
    return df.loc[:, [col for col in df.columns if not str(col).startswith('Unnamed') and str(col) != '']]


def parse_csv_list(value):
    """Normalize a CSV field into a list of stripped string values."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    if isinstance(value, str):
        value = value.strip()
        if value == '':
            return []
        return [item.strip() for item in value.split(',') if item.strip() != '']
    return [str(value).strip()]


def parse_csv_bool(value, default=False):
    """Normalize CSV booleans stored as 0/1, true/false, or blank."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return default
    if isinstance(value, str):
        value = value.strip().lower()
        if value == '':
            return default
        if value in ['1', 'true', 't', 'yes', 'y']:
            return True
        if value in ['0', 'false', 'f', 'no', 'n']:
            return False
    return bool(value)


def get_struct_variant_name(struct_name, pdb_name=None, target_chain=None, keep_multiple_chains_in_struct=True):
    """Return the saved scene basename for either combined or per-chain scenes."""
    if keep_multiple_chains_in_struct or target_chain in [None, '']:
        return struct_name
    if pdb_name in [None, '']:
        pdb_name = struct_name
    pdb_core = str(pdb_name).split('_', 1)[0]
    if pdb_core not in struct_name:
        return f'{struct_name}-{target_chain}'
    return struct_name.replace(pdb_core, f'{pdb_core}-{target_chain}', 1)


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

def findProcess(process_name):
    if opsys == 'Windows':
        return [int(item.split()[1]) for item in os.popen('tasklist').read().splitlines()[4:] if process_name in item.split()]
    elif opsys == 'Linux' or opsys == 'Darwin':
        return [int(pid) for pid in os.popen('pidof ' + process_name).read().strip(' \n').split(' ') if pid != '']
    

def exit_program(pid):
    import signal
    print("Sending SIGINT to self...")
    os.kill(pid, signal.SIGINT)
    print('Exited program', pid)

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
            df = pd.read_csv(log_fpath)
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
    output_csv_fpath = os.path.join(output_dir, output_fname + '.csv')
    if df_all is not None:
        if os.path.exists(output_csv_fpath) and append_to_existing_output:
            df_existing = pd.read_csv(output_csv_fpath)
            df_all = pd.concat([df_existing, df_all], axis=0)
        df_all.to_csv(output_csv_fpath, index=False)

    # record missing files
    if missing_data:
        with open(os.path.join(output_dir, 'missing_data.txt'), 'w') as f:
            missing_txt = '\n'.join(missing_data) + '\n'
            f.write(missing_txt)

    # remove combined files
    if remove_combined_files:
        for log_fpath in [f for f in log_fpath_list if f not in missing_data]:
            os.remove(log_fpath)
    return missing_data
