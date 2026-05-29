# YASARA ddG Bind Workflow

This repo stores structures, YASARA scene files, input tables, and output tables for running mutagenesis-based binding energy (`ddG`) calculations with YASARA.

## Repo Layout

- `tools/prepare_yasara_sce.py`: PDB -> SCE preparation driver
- `tools/run_yasara_mutagenesis_binding.py`: mutagenesis / ddG workflow driver
- `tools/configs.py`: shared dataset and run configuration
- `tools/utils.py`: mutation parsing and CSV / path helpers
- `tools/variables.py`: repo paths and residue / ligand constants
- `influenza-resistance/pdb/`: raw input PDB structures
- `influenza-resistance/sce/`: prepared YASARA `.sce` structures
- `influenza-resistance/yasara/Input/`: batch input CSV files
- `influenza-resistance/yasara/Output/`: ddG outputs and post-optimization structures

## Requirements

- Python 3
- YASARA installed locally
- Python packages used by the workflow:
  - `numpy`
  - `pandas`
- `tools/yasara.py` must point to the correct local YASARA installation directory

## Workflow Overview

The workflow is now explicitly two-stage:

1. `tools/prepare_yasara_sce.py`
   - reads raw PDB files
   - removes water
   - runs `CleanAll()`
   - keeps the requested receptor chains, ligand, and non-ligand cofactors
   - writes prepared `.sce` files to the parallel `sce/` directory

2. `tools/run_yasara_mutagenesis_binding.py`
   - loads prepared `.sce` files only
   - expands jobs over mutations, parameter combinations, and `target_chain`
   - performs WT / MT swaps and minimization
   - computes ddG outputs
   - optionally exports minimized `.pdb` files from saved post-optimization `.sce` files

## Configuration

Both scripts read shared settings from [tools/configs.py](/Users/charmainechia/Documents/projects/yasara-ddGbind-workflow/tools/configs.py).

### `DATA_CONFIG`

Controls which dataset is read:

- `data_folder`: root data directory
- `data_subfolder`: optional subfolder under `pdb/`, `sce/`, and `yasara/Output/`
- `inputs`: input CSV filename from `yasara/Input/`
- `output_fname`: basename for combined output CSVs

### `RUN_CONFIG`

Controls run-time behavior for the ddG runner:

- `run_multiprocessing`: `None` for serial mode, or an integer worker count
- `save_minimized_struct`: whether to export post-optimization `.pdb` files
- `fix_metal_ion`: whether to keep metal ions fixed during minimization
- `append_to_existing_output`: whether combined CSV output appends to an existing output file
- `sep`: separator for combined mutation strings
- `energy_calc_method`: currently `BoundaryFast` or `PBS`

### `RUN_PARAMS`

Defines the parameter sweep. Each entry is a list, and the runner expands all combinations:

- `nrep`
- `minimize_energy`
- `resetSce`
- `move`
- `ff`
- `mvdist`
- `mvdrug`
- `surfout`

## Input CSV Format

The workflow reads a CSV from `influenza-resistance/yasara/Input/`.

Common columns:

- `struct_name`: base name for the prepared `.sce` and downstream outputs
- `pdb_id`: raw input PDB basename, without `.pdb`
- `ligand_name`: label used in output naming
- `process_structure`: `1` to run, `0` to skip
- `mutations`: comma-separated mutations; combined mutations should use the configured separator
- `chain_id`: comma-separated receptor chains to retain from the source PDB
- `target_chain`: comma-separated subset of `chain_id` to actually mutate / evaluate
- `keep_multiple_chains_in_struct`: `1` to keep all listed receptor chains in one `.sce`, `0` to save one `.sce` per `target_chain`
- `ligand_id`: ligand residue name in the structure
- `ligand_chain_id`: optional ligand chain ID override for prep

### Important naming behavior

`struct_name` and `pdb_id` do different jobs:

- `pdb_id` determines which raw PDB file is loaded
  - example: `pdb_id = 6fs6_Baloxavir` loads `.../pdb/6fs6_Baloxavir.pdb`
- `struct_name` determines the prepared scene / run basename

When `keep_multiple_chains_in_struct=0`, per-chain scene naming is derived from `struct_name`, but the chain suffix is inserted after the PDB core token from `pdb_id`.

Example:

- `pdb_id = 6fs6_Baloxavir`
- `struct_name = PA_pH1N1_6fs6_Baloxavir`
- `target_chain = A`

Prepared scene name:

- `PA_pH1N1_6fs6-A_Baloxavir.sce`

not:

- `PA_pH1N1_6fs6_Baloxavir-A.sce`

## Chain Handling

### `keep_multiple_chains_in_struct=1`

One combined scene is prepared.

- `chain_id` controls which receptor chains remain in the structure
- all retained non-ligand molecules stay in `Obj 1`
- all retained ligand molecules are grouped into `Obj 2`
- the runner expands one job per `target_chain`, but reuses the same combined `.sce`

Use this when you want to keep a multichain assembly intact, such as a homo-oligomer, while still computing one chain-specific run at a time.

### `keep_multiple_chains_in_struct=0`

One chain-specific scene is prepared per requested `target_chain`.

- `chain_id` defines which receptor chains are allowed to be used
- `target_chain` defines which chain-specific scenes are actually written and run
- each saved scene name is derived from `struct_name` plus the chain suffix inserted after the PDB core token

Use this when you want separate SCE files for each chain calculation.

### `target_chain`

`target_chain` is always the chain that the runner mutates and evaluates in a single job.

If the CSV omits `target_chain`, the runner falls back to `chain_id`.

## How To Run

### 1. Prepare `.sce` files from PDB input

Set the dataset in [tools/configs.py](/Users/charmainechia/Documents/projects/yasara-ddGbind-workflow/tools/configs.py), then run:

```bash
python3 tools/prepare_yasara_sce.py
```

This script:

- loads raw PDBs from `.../pdb/<data_subfolder>/` or `.../pdb/`
- always runs `DelWater()` and `CleanAll()`
- writes prepared scenes to the parallel `.../sce/<data_subfolder>/` or `.../sce/`

### 2. Run the ddG workflow from prepared `.sce` files

Set `DATA_CONFIG`, `RUN_CONFIG`, and `RUN_PARAMS` in [tools/configs.py](/Users/charmainechia/Documents/projects/yasara-ddGbind-workflow/tools/configs.py), then run:

```bash
python3 tools/run_yasara_mutagenesis_binding.py
```

The runner:

1. reads the input CSV
2. groups jobs by `(struct_name, ligand_name)`
3. expands each row across all requested `target_chain` values
4. expands parameter combinations from `RUN_PARAMS`
5. loads the prepared `.sce`
6. applies WT / MT residue swaps only on the current `target_chain`
7. minimizes the scene for each replicate
8. computes chain-specific receptor / ligand / complex energies
9. combines replicate outputs into AVG and FULL CSV files
10. optionally converts post-optimization `.sce` files to `.pdb`

## Outputs

Results are written under:

- `influenza-resistance/yasara/Output/<data_subfolder>/`

Typical outputs:

- `<output_fname>.csv`: combined averaged results
- `<output_fname>_FULL.csv`: combined replicate-level results
- `postOpt/`: saved post-minimization `.sce` files and optional exported `.pdb` files
- `missing_data.txt`: missing intermediate files during the combine step

Worker-level temporary CSVs use internal `DDG_<struct>_<target_chain>...` naming and are combined at the end of the run.

## Multiprocessing Notes

The runner supports multiprocessing through `RUN_CONFIG['run_multiprocessing']`.

Current behavior:

- serial mode: set `run_multiprocessing = None`
- parallel mode: set `run_multiprocessing` to an integer worker count

Notes:

- the helper `parse_args_process_mutant(...)` is retained specifically for multiprocessing
- temporary worker CSV naming includes `target_chain`, so parallel outputs do not collide
- `.sce -> .pdb` export happens once after all workers finish, rather than inside each worker
- if a worker raises an uncaught exception, the whole run will stop
- YASARA multiprocessing also depends on your local YASARA installation and license behavior

## Notes

- The runner assumes prepared `.sce` files exist before execution.
- The workflow assumes YASARA object numbering remains compatible with the prep layout:
  - `Obj 1`: receptor plus retained non-ligand molecules
  - `Obj 2`: ligand molecules
- If a structure has missing loops near the mutation site or ligand, reconstructing that region before running the workflow is usually safer than leaving a broken local structure in place.
