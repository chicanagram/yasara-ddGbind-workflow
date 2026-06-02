#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from collections import OrderedDict

from variables import address_dict


YASARA_CONFIG = {
    'version_suffix': '2026',
}


DATA_CONFIG = {
    'data_folder': address_dict['influenza-resistance'],
    'data_subfolder': '',
    'inputs': 'PA-NA_benchmark.csv',
    'output_fname': f'PA-NA_benchmark_yasara{YASARA_CONFIG["version_suffix"]}',
}


RUN_CONFIG = {
    'run_multiprocessing': 16,
    'save_minimized_struct': True,
    'fix_metal_ion': True,
    'append_to_existing_output': True,
    'sep': '+',
    'energy_calc_method': 'BoundaryFast',
    'energy_selection_mode': 'chain-aware',
}


RUN_PARAMS = OrderedDict()
RUN_PARAMS['nrep'] = [10]
RUN_PARAMS['minimize_energy'] = [True]
RUN_PARAMS['resetSce'] = [False]
RUN_PARAMS['move'] = ['!backbone']
RUN_PARAMS['ff'] = ['YASARA2', 'AMBER14', 'AMBER15FB'] # ['YASARA2']
RUN_PARAMS['mvdist'] = [4]
RUN_PARAMS['mvdrug'] = [1]
RUN_PARAMS['surfout'] = [0.65]
