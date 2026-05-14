# YASARA ddG Bind Workflow

This repo stores structures, YASARA scene files, input tables, and output tables for running mutagenesis-based binding energy (`ddG`) calculations with YASARA.

## Repo layout

- `tools/run_yasara_mutagenesis_binding.py`: main workflow driver
- `tools/utils.py`: mutation parsing, PDB to SCE conversion, CSV helpers
- `tools/variables.py`: repo paths and residue/ligand constants
- `influenza-resistance/pdb/`: input PDB structures
- `influenza-resistance/sce/`: input YASARA `.sce` structures
- `influenza-resistance/yasara/Input/`: batch input CSV files
- `influenza-resistance/yasara/Output/`: workflow outputs and post-optimization structures

## Requirements

- Python 3
- YASARA installed locally
- Python packages used by the workflow:
  - `numpy`
  - `pandas`
- `tools/yasara.py` must point to the correct local YASARA installation directory

## Input format

The workflow reads a CSV from `influenza-resistance/yasara/Input/`.

Common columns:

- `struct_name`: structure basename without extension
- `ligand_name`: ligand label used in output naming
- `process_structure`: `1` to run, `0` to skip
- `mutations`: comma-separated mutations, with combined mutations joined by the configured separator
- `chain_id`: receptor chain(s) to keep when starting from PDB
- `ligand_id`: ligand residue name in the structure
- `ligand_chain_id`: optional ligand chain ID

Examples already in the repo include:

- `influenza-resistance/yasara/Input/NA.csv`
- `influenza-resistance/yasara/Input/PA-A37T-I38T.csv`
- `influenza-resistance/yasara/Input/PA_benchmark_mod.csv`

## How to run

1. Open `tools/run_yasara_mutagenesis_binding.py`.
2. In the `__main__` block, set:
   - `data_folder` to `address_dict['influenza-resistance']`
   - `data_subfolder` to the dataset folder you want to run
   - `struct_format` to `sce` or `pdb`
   - `inputs` to a CSV present in `influenza-resistance/yasara/Input/`
   - `output_fname` to the output prefix you want
3. Adjust run parameters if needed:
   - `nrep`
   - `ff`
   - `move`
   - `mvdist`
   - `mvdrug`
   - `surfout`
   
   These variables are defined as lists in the script. The workflow expands
   them into combinations, so you can sweep across multiple parameter values
   in one run. If you want only one value for a given parameter, keep it as a
   single-element list such as `params['ff'] = ['YASARA2']`.
4. Run:

```bash
python3 tools/run_yasara_mutagenesis_binding.py
```

## Workflow summary

The script:

1. Reads the input CSV and collects structure / ligand / mutation jobs.
2. If `struct_format='pdb'`, converts PDB files into YASARA `.sce` files.
3. Loads each structure in YASARA.
4. Applies WT and mutant residue swaps.
5. Minimizes the structure for each replicate.
6. Computes component energies for complex, receptor, and ligand.
7. Writes per-run and aggregated CSV outputs.
8. Optionally exports minimized post-optimization structures.

## Outputs

Results are written under `influenza-resistance/yasara/Output/<data_subfolder>/`.

Typical outputs:

- `DDG_<name>.csv`: averaged results
- `DDG_<name>_FULL.csv`: all replicate-level results
- `postOpt/`: saved post-minimization `.sce` and optional `.pdb` files
- `missing_data.txt`: missing intermediate files during combine step

## Notes

- The workflow assumes YASARA object numbering is consistent with the expected receptor / ligand / complex layout.
- If you start from `pdb`, make sure ligand and chain annotations are correct in the input CSV.
- Large output folders in this repo include archived historical runs as well as current working data.
