# Organometallics_Fragmentation
This Python script processes molecular structures from an .extxyz file, identifying and removing CO ligands bonded to metal centers. It uses ASE and OpenBabel for ligand identification and modification, saving the updated molecules to an output file. 

# CO Ligand Removal Script

This Python script, `co_ligand_removal.py`, is designed to process molecular structures stored in an `.extxyz` file, specifically identifying and removing CO ligands from metal complexes. The script uses ASE to handle molecular geometries and OpenBabel for ligand identification and processing. 

## Features
- **CO Ligand Identification**: Detects CO ligands in molecules where a carbon atom is bonded to both a metal atom and an oxygen atom.
- **Ligand Removal**: Removes one CO ligand (a carbon and oxygen pair) from each molecule that contains such a ligand.
- **File Processing**: Reads molecular geometries from an input file, processes them, and saves the modified structures to an output file. Any failed molecules are logged separately.
- **Progress Tracking**: Uses `tqdm` to display the progress of processing molecules in real-time.

## Requirements
The script requires the following Python packages:
- `ase`
- `numpy`
- `openbabel`
- `tqdm`

Install these dependencies via pip:
```bash
pip install ase numpy openbabel tqdm
```

## Input and Output Files

- **Input File**: 
  The input file (`tmQM_neutral.extxyz`) is expected to be an extended XYZ file format containing molecular structures.
  
- **Output File**: 
  The script generates an output file (`tmQM_neutral_CO_diss_subset.extxyz`) containing the modified molecular structures with one CO ligand removed (if present).

- **Failed Molecules File**: 
  Molecules that cannot be processed are logged in the `failed_molecules.txt` file.

## Usage

Run the script from the command line after setting up the necessary input files. The script processes each molecule in the input file, modifies it if a CO ligand is present, and writes the results to the output file.

```bash
python co_ligand_removal.py
```

## Functions Overview

1. **suppress_kekulization_warnings()**: Suppresses OpenBabel kekulization warnings for cleaner output.
2. **atoms_to_xyz_string(atoms)**: Converts an ASE Atoms object into an XYZ string format with the molecule's CSD code as a comment.
3. **identify_co_ligands(ob_mol, metal_elements)**: Identifies CO ligands in a molecule by searching for carbon atoms bonded to both a metal atom and an oxygen atom.
4. **process_molecule(csd_code, num_atoms, atoms, outfile, failfile)**: Processes each molecule to detect and remove one CO ligand if present.

## File Structure

- `tmQM_neutral.extxyz`: Input file with molecular geometries.
- `tmQM_neutral_CO_diss_subset.extxyz`: Output file containing molecules with one CO ligand removed (if detected).
- `failed_molecules.txt`: Log file for molecules that failed during processing.

## Progress Tracking

The script uses `tqdm` to provide a real-time progress bar during molecule processing.

## License

This project is open-source and available under the MIT License.
