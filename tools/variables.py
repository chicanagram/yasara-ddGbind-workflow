#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024

@author: charmainechia
"""
address_dict = {
    'influenza-resistance': '../influenza-resistance/',
}

subfolders = {
    'sequences': 'sequences/',
    'yasara': 'yasara/',
    'pdb': 'pdb/',
    'sce': 'sce/',
}

mapping = {
    'A': 'Ala',
    'H': 'His',
    'Y': 'Tyr',
    'R': 'Arg',
    'T': 'Thr',
    'K': 'Lys',
    'M': 'Met',
    'D': 'Asp',
    'N': 'Asn',
    'C': 'Cys',
    'Q': 'Gln',
    'E': 'Glu',
    'G': 'Gly',
    'I': 'Ile',
    'L': 'Leu',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'W': 'Trp',
    'V': 'Val'
    }

hetatm_non_metal_ion = ['HEM', 'NAG']
hetatm_metal_ion = [' MG', '  K', ' NA', ' MN', ' CA']

aaList = list("ACDEFGHIKLMNPQRSTVWY")
aaList_with_X = list("ACDEFGHIKLMNPQRSTVWYX")
aa2idx = {aa: i for i, aa in enumerate(aaList)}


nonstandard_amino_acids = [
    # Seleno and special amino acids
    "MSE",  # Selenomethionine
    "SEC",  # Selenocysteine
    "PYL",  # Pyrrolysine
    # Post-translational modifications (PTMs)
    "TPO",  # Phosphothreonine
    "SEP",  # Phosphoserine
    "PTR",  # Phosphotyrosine
    "HYP",  # Hydroxyproline
    "MLY",  # N6-Methyllysine
    "ALY",  # Nε-Acetyllysine
    # Cysteine variants
    "CSS",  # Thioether-cysteine (disulfide-linked)
    "CSO",  # S-hydroxycysteine
    "CSD",  # S-sulfinylcysteine
    "CME",  # S-methylcysteine
    # Others
    "FME",  # N-formylmethionine
    "AIB",  # α-Aminoisobutyric acid
    "ORN",  # Ornithine
    "DPR",  # D-Proline
    "HIC",  # 4-Hydroxyisoleucine
]
