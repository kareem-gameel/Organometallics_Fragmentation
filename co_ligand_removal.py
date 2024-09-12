from ase.io import read, write
from ase import Atoms
import numpy as np
from openbabel import openbabel
import warnings
from tqdm import tqdm

# Define the input and output files
input_file = '../geometries/tmQM_neutral.extxyz'
output_file = '../geometries/tmQM_neutral_CO_diss_subset.extxyz'
failed_file = '../geometries/failed_molecules.txt'

# Suppress specific OpenBabel warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning, module="openbabel")

# Suppress Open Babel kekulization warnings
def suppress_kekulization_warnings():
    openbabel.obErrorLog.StopLogging()
    openbabel.obErrorLog.SetOutputLevel(0)

suppress_kekulization_warnings()

# List of metal elements symbols
metal_elements = [
    'Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Cs', 'Ba',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
    'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U',
    'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
    'Ds', 'Rg', 'Cn'
]

# Function to convert ASE Atoms object to XYZ string with CSD code
def atoms_to_xyz_string(atoms):
    # Start with the number of atoms and CSD code as the comment line
    xyz_str = f"{len(atoms)}\n{atoms.info.get('CSD_code', '')}\n"
    # Add each atom's symbol and position to the string
    for symbol, position in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        xyz_str += f"{symbol} {position[0]:.8f} {position[1]:.8f} {position[2]:.8f}\n"
    return xyz_str

# Function to identify CO ligands in the molecule
def identify_co_ligands(ob_mol, metal_elements):
    co_ligands = []
    for ob_atom in openbabel.OBMolAtomIter(ob_mol):
        if ob_atom.GetAtomicNum() == 6:  # Check if the atom is carbon
            neighbors = list(openbabel.OBAtomAtomIter(ob_atom))
            #print(f"Carbon atom {ob_atom.GetIdx()} has {len(neighbors)} neighbors")
            if len(neighbors) == 2:  # Ensure the carbon is bonded only to two atoms
                metal_center = None
                oxygen_atom = None
                for neighbor in neighbors:
                    atomic_num = neighbor.GetAtomicNum()
                    element_symbol = neighbor.GetType()[0:2].strip()
                    #print(f"  Neighbor {neighbor.GetIdx()} has atomic number {atomic_num} and element symbol {element_symbol}")
                    if element_symbol in metal_elements:  # Check if the neighbor is a metal
                        metal_center = neighbor
                        #print(f"  Metal center found: {neighbor.GetIdx()} ({element_symbol})")
                    elif atomic_num == 8 and neighbor.GetTotalDegree() == 1:  # Check if the neighbor is an oxygen bonded only to the carbon
                        oxygen_atom = neighbor
                        #print(f"  Oxygen atom found: {neighbor.GetIdx()} ({atomic_num})")
                if metal_center and oxygen_atom:  # Ensure both metal and oxygen are present
                    #print(f"  CO ligand identified: C {ob_atom.GetIdx()} - O {oxygen_atom.GetIdx()} - Metal {metal_center.GetIdx()}")
                    co_ligands.append((ob_atom.GetIdx(), oxygen_atom.GetIdx()))
    return co_ligands

# Function to process each molecule
def process_molecule(csd_code, num_atoms, atoms, outfile, failfile):
    try:
        # Create ASE Atoms object
        symbols = [atom[0] for atom in atoms]
        positions = np.array([atom[1:] for atom in atoms])
        molecule = Atoms(symbols=symbols, positions=positions)
        molecule.info['CSD_code'] = csd_code

        # Convert ASE Atoms object to XYZ string
        xyz_string = atoms_to_xyz_string(molecule)

        # Convert XYZ string to OpenBabel OBMol object
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("xyz", "mol")
        ob_mol = openbabel.OBMol()
        if not obConversion.ReadString(ob_mol, xyz_string):
            raise ValueError(f"{csd_code}: Failed to read XYZ string with OpenBabel")

        # Find CO ligands
        co_ligands = identify_co_ligands(ob_mol, metal_elements)

        # Remove one CO ligand from the molecule
        if co_ligands:
            # Select the first CO ligand to remove
            atom_idx, oxygen_idx = co_ligands[0]

            # Get the positions of the atoms to be removed
            positions = molecule.get_positions()
            elements = molecule.get_chemical_symbols()

            # Remove the atoms from ASE Atoms object
            new_positions = np.delete(positions, [atom_idx - 1, oxygen_idx - 1], axis=0)
            new_elements = np.delete(elements, [atom_idx - 1, oxygen_idx - 1])

            # Create a new ASE Atoms object with the remaining atoms
            new_molecule = Atoms(symbols=new_elements, positions=new_positions)
            new_molecule.info['CSD_code'] = molecule.info.get('CSD_code', '')

            # Save the modified molecule to the output file
            outfile.write(atoms_to_xyz_string(new_molecule))
            return 'success'
        else:
            return 'no_co'

    except Exception as e:
        failfile.write(f"{csd_code}: {str(e)}\n")
        return 'failed'

# Read the input file and process each molecule
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile, open(failed_file, 'w') as failfile:
    lines = infile.readlines()
    num_lines = len(lines)
    i = 0
    failed_count = 0
    success_count = 0
    no_co_count = 0
    total_processed = 0

    with tqdm(total=num_lines, desc="Processing molecules") as pbar:
        while i < num_lines:
            try:
                num_atoms = int(lines[i].strip())
                csd_code_line = lines[i + 1].strip()
                csd_code = csd_code_line.split()[0] if "CSD_code=" not in csd_code_line else csd_code_line.split("CSD_code=")[1].split()[0]
                atoms = []
                for j in range(num_atoms):
                    parts = lines[i + 2 + j].split()
                    atoms.append([parts[0]] + [float(x) for x in parts[1:4]])
                i += num_atoms + 2

                result = process_molecule(csd_code, num_atoms, atoms, outfile, failfile)
                if result == 'success':
                    success_count += 1
                elif result == 'no_co':
                    no_co_count += 1
                else:
                    failed_count += 1

                total_processed += 1
                pbar.update(num_atoms + 2)

            except Exception as e:
                failfile.write(f"{csd_code}: {str(e)}\n")
                failed_count += 1
                total_processed += 1
                i += 1  # move to the next line
                pbar.update(1)

    print(f"Total molecules processed: {total_processed}")
    print(f"Total successful CO dissociations: {success_count}")
    print(f"Total failed molecule conversions: {failed_count}")
    print(f"Total molecules without CO detected: {no_co_count}")
    print(f"Modified molecules with one CO ligand removed (if present) have been saved to {output_file}")
    print(f"List of failed molecules saved to {failed_file}")
