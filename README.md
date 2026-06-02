# YASARA ddG Bind Workflow

This repo stores structures, YASARA scene files, input tables, and output tables for running mutagenesis-based binding energy (`ddG`) calculations with YASARA.

## Quickstart

1. Set the active YASARA version in [tools/configs.py](/Users/charmainechia/Documents/projects/yasara-ddGbind-workflow/tools/configs.py):
   - `YASARA_CONFIG['version_suffix']`
   <br>This suffix determines the YASARA folder name to access (e.g. `yasara2026`), and can also be set to an empty string if applicable (i.e. `""`)
2. Set the dataset in `DATA_CONFIG`:
   - `data_folder`
   - `data_subfolder`
   - `inputs`
   - `output_fname`
3. Prepare `.sce` files from raw PDBs:

```bash
python3 tools/prepare_yasara_sce.py
```

4. Set runner options in `RUN_CONFIG` and `RUN_PARAMS`.
5. Run the ddG workflow:

```bash
python3 tools/run_yasara_mutagenesis_binding.py
```

The workflow is two-stage:
- `prepare_yasara_sce.py` converts raw PDB inputs into prepared `.sce` scenes
- `run_yasara_mutagenesis_binding.py` loads those `.sce` files and performs WT/MT ddG calculations

## Repo Layout

- `tools/prepare_yasara_sce.py`: PDB -> SCE preparation driver
- `tools/run_yasara_mutagenesis_binding.py`: mutagenesis / ddG workflow driver
- `tools/configs.py`: shared YASARA, dataset, and run configuration
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
- [tools/yasara.py](/Users/charmainechia/Documents/projects/yasara-ddGbind-workflow/tools/yasara.py) must point to a valid local YASARA installation layout

## Configuration

Both scripts read shared settings from [tools/configs.py](/Users/charmainechia/Documents/projects/yasara-ddGbind-workflow/tools/configs.py).

### `YASARA_CONFIG`

- `version_suffix`: used by [tools/yasara.py](/Users/charmainechia/Documents/projects/yasara-ddGbind-workflow/tools/yasara.py) to build the local YASARA install path

### `DATA_CONFIG`

- `data_folder`: root data directory
- `data_subfolder`: optional subfolder under `pdb/`, `sce/`, and `yasara/Output/`
- `inputs`: input CSV filename from `yasara/Input/`
- `output_fname`: basename for combined output CSVs

### `RUN_CONFIG`

- `run_multiprocessing`: `None` for serial mode, or an integer worker count
- `save_minimized_struct`: whether to export post-optimization `.pdb` files
- `fix_metal_ion`: whether to keep metal ions fixed during minimization
- `append_to_existing_output`: whether combined CSV output appends to an existing output file
- `sep`: separator for combined mutation strings
- `energy_calc_method`: currently `BoundaryFast` or `PBS`

### `RUN_PARAMS`

These lists define the parameter sweep:

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

### Naming behavior

`pdb_id` selects the raw PDB to load, while `struct_name` controls the prepared scene / run basename.

When `keep_multiple_chains_in_struct=0`, the chain suffix is inserted after the PDB core token derived from `pdb_id`.

Example:

- `pdb_id = 6fs6_Baloxavir`
- `struct_name = PA_pH1N1_6fs6_Baloxavir`
- `target_chain = A`

Prepared scene name:

- `PA_pH1N1_6fs6-A_Baloxavir.sce`

## Chain Handling

### `keep_multiple_chains_in_struct=1`

- one combined scene is prepared
- `chain_id` controls which receptor chains remain in the structure
- non-ligand molecules remain in `Obj 1`
- retained ligand molecules are grouped into `Obj 2`
- the runner expands one job per `target_chain`, reusing the same combined `.sce`

Use this when you want to keep a multichain assembly intact while still calculating one chain-specific run at a time.

### `keep_multiple_chains_in_struct=0`

- one chain-specific scene is prepared per requested `target_chain`
- `chain_id` defines the allowable receptor chains
- `target_chain` defines which per-chain scenes and jobs are created

Use this when you want separate SCE files for separate chain calculations.

If `target_chain` is omitted, the runner falls back to `chain_id`.

## Prep Step

Run:

```bash
python3 tools/prepare_yasara_sce.py
```

The prep script:

- loads raw PDBs from `.../pdb/<data_subfolder>/` or `.../pdb/`
- runs `DelWater()` and `CleanAll()`
- keeps the requested receptor chains, ligand, and retained non-ligand molecules
- writes prepared scenes to the parallel `.../sce/<data_subfolder>/` or `.../sce/`
- prints a footer summary for any prepared scenes with unexpected object counts

## Runner Step

Run:

```bash
python3 tools/run_yasara_mutagenesis_binding.py
```

The runner:

1. reads the input CSV
2. groups jobs by `(struct_name, ligand_name)`
3. expands each row across requested `target_chain` values
4. expands all parameter combinations from `RUN_PARAMS`
5. loads the prepared `.sce`
6. applies WT / MT swaps on the current `target_chain`
7. minimizes the scene for each replicate
8. calculates receptor / ligand / complex energies
9. combines replicate outputs into AVG and FULL CSVs
10. optionally converts post-optimization `.sce` files to `.pdb`

## Outputs

Results are written under:

- `influenza-resistance/yasara/Output/<data_subfolder>/`

Typical outputs:

- `<output_fname>.csv`: combined averaged results
- `<output_fname>_FULL.csv`: combined replicate-level results
- `postOpt/`: saved post-minimization `.sce` files and optional exported `.pdb` files
- `failed_jobs.csv`: jobs that failed validation or runtime execution
- `missing_data.txt`: missing intermediate files during combine

Worker-level temporary CSVs use internal `DDG_<struct>_<target_chain>...` naming and are combined at the end of the run.

## Multiprocessing Notes

The runner supports multiprocessing through `RUN_CONFIG['run_multiprocessing']`.

- serial mode: set `run_multiprocessing = None`
- parallel mode: set `run_multiprocessing` to an integer worker count

Notes:

- `parse_args_process_mutant(...)` is retained specifically for multiprocessing
- temporary worker CSV naming includes `target_chain`, so parallel outputs do not collide
- `.sce -> .pdb` export happens once after all workers finish
- failed jobs are logged and the rest of the batch continues
- YASARA multiprocessing still depends on local YASARA installation and license behavior

## Notes

- The runner assumes prepared `.sce` files already exist before execution.
- The workflow assumes YASARA object numbering remains compatible with the prep layout:
  - `Obj 1`: receptor plus retained non-ligand molecules
  - `Obj 2`: ligand molecules
- If a structure has missing loops near the mutation site or ligand, reconstructing that region before running the workflow is usually safer than leaving a broken local structure in place.
