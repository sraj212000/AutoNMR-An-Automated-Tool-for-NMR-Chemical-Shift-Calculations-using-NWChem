from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
import os
import re
import subprocess
import time
import sys
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from rdkit import Chem

def generate_conformers(smiles, num_conformers):
    """
    Generate conformers for a molecule from SMILES and return them with energies.

    Args:
        smiles (str): SMILES string of the molecule.
        num_conformers (int): Number of conformers to generate.

    Returns:
        tuple: (sdf_data, energies)
            sdf_data (str): SDF file content containing all conformers.
            energies (list): List of tuples containing (conformer_id, energy).
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # Add hydrogen atoms

    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol)
    
    # Optimize geometry using UFF
    AllChem.UFFOptimizeMolecule(mol)
    
    # Generate conformers
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers)
    
    # Optimize each conformer using UFF and calculate energies
    energies = []
    for cid in cids:
        AllChem.UFFOptimizeMolecule(mol, confId=cid)
        energy = AllChem.UFFGetMoleculeForceField(mol, confId=cid).CalcEnergy()
        energies.append((cid, energy))
    
    # Sort conformers by energy
    energies.sort(key=lambda x: x[1])
    
    # Write conformers to an SDF file in a temporary directory
    with tempfile.TemporaryDirectory() as tmp_dir:
        sdf_file = os.path.join(tmp_dir, 'conformers.sdf')
        writer = Chem.SDWriter(sdf_file)
        for cid, energy in energies:
            mol.SetProp('_Name', f'Conformer_{cid}_Energy_{energy:.2f}')
            writer.write(mol, confId=cid)
        writer.close()
        
        # Read the SDF file content
        with open(sdf_file, 'rb') as f:
            sdf_data = f.read().decode('utf-8')
    
    return sdf_data, energies

def generate_3d_coordinates(smiles, sdf_data, xc_functional, output_file_prefix, case):
    """
    Generates 3D geometry input files for NWChem based on optimized conformers.

    Args:
        smiles (str): The SMILES string of the molecule (for information purposes).
        sdf_data (str): SDF file content containing all optimized conformers.
        xc_functional (str): The exchange-correlation functional to use in the DFT calculations.
        output_file_prefix (str): The prefix for the output NWChem input files.

    Returns:
        list: A list of generated NWChem input file paths.
    """
    atomic_number_to_symbol = {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B',
        6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
        11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
        16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
        # Extend this mapping as needed for other elements
    }

    atomic_numbers_to_include = {
        1: 'Hydrogen',  # Case 1: Hydrogen atoms only
        6: 'Carbon'    # Case 2: Carbon atoms only
    }

    if case not in atomic_numbers_to_include:
        raise ValueError("Invalid case selection. Choose 1 for Hydrogen or 2 for Carbon.")

    target_atomic_number = case
    target_element = atomic_numbers_to_include[target_atomic_number]

    suppl = Chem.SDMolSupplier('conformers.sdf', removeHs=False)

    nwchem_input_files = []

    for mol_idx, mol in enumerate(suppl):
        if mol is None:
            print(f"Error: Could not read molecule from SDF data.")
            continue

        for conf_idx in range(mol.GetNumConformers()):
            conformer = mol.GetConformer(conf_idx)

            nwchem_input_file = f"{output_file_prefix}_mol_{mol_idx}_conformer_{conf_idx}.nw"
            nwchem_input_files.append(nwchem_input_file)

            with open(nwchem_input_file, 'w') as f:
                f.write(f"start\n\ntitle \"Molecule conformer {conf_idx}\"\n\nmemory total 4000 mb\n\necho\ngeometry units angstroms\n")
                
                atom_indices = []

                for atom_idx in range(mol.GetNumAtoms()):
                    atom = mol.GetAtomWithIdx(atom_idx)
                    atomic_num = atom.GetAtomicNum()
                    symbol = atomic_number_to_symbol.get(atomic_num, 'Xx')  # 'Xx' for unknown symbols
                    pos = conformer.GetAtomPosition(atom_idx)
                    x, y, z = pos.x, pos.y, pos.z
                    f.write(f"   {symbol:<2s} {x:>10.4f} {y:>10.4f} {z:>10.4f}\n")
                    # Collect indices of target atoms
                    if atomic_num == target_atomic_number:
                        atom_indices.append(atom_idx + 1)  # 1-based index for NWChem

                f.write(
                    "end\n\n\nbasis spherical\n* library 6-311G\nend\n\n"
                    "driver\n"
                    "  maxiter 500\n  xyz final\n"
                    "end\n\n"
                    "relativistic\n"
                    "  zora on\n"
                    "  zora:cutoff_NMR 1d-8\n"
                    "  zora:cutoff 1d-30\n"
                    "end\n\n"
                    "dft\n"
                    "  direct\n"
                    "  grid fine\n"
                    f"  xc {xc_functional}\n"
                    "  mult 1\n"
                    "  maxiter 1500\n"
                    '  noprint "final vectors analysis" multipole\n'
                    "end\n\n"
                    "task dft optimize\n\n"
                    "task dft freq\n\n"
                    "property\n"
                )

                # Specify target atoms for shielding tensor calculation
                if atom_indices:
                    f.write(f"  shielding {len(atom_indices)} {' '.join(map(str, atom_indices))}\n")
                
                f.write(
                    "end\n\n"
                    "cosmo\n"
                    "   solvent cdcl3\n"
                    "end\n\n"
                    "task dft property\n"
                )

    return nwchem_input_files
    
    
   
def run_nwchem(commands):
    """
    Executes an NWChem command using the system shell and measures its execution time.

    This function:
    1. Prints a message indicating that NWChem is running.
    2. Records the start time of the NWChem command.
    3. Executes the provided NWChem command using the system shell.
    4. Records the end time of the NWChem command.
    5. Calculates and prints the execution time of the NWChem command.

    Args:
        commands (str): The NWChem command to be executed.

    Notes:
        - The function uses `subprocess.run` to execute the command.
        - The `shell=True` argument allows the command to be executed through the shell, which is necessary for complex commands.
        - The execution time is printed to the standard error stream.

    Example:
        run_nwchem("mpiexec --use-hwthread-cpus -np 40 nwchem input_file.nw > output_file.out")
    """
    print("Running NWChem...")
    start_time = time.time()  # Record start time
    subprocess.run(commands, shell=True)  # Execute the command
    end_time = time.time()  # Record end time
    execution_time = end_time - start_time  # Calculate execution time
    print("NWChem execution time: {:.2f} seconds".format(execution_time), file=sys.stderr)  # Print execution time

def run(smiles):
    """
    Generates NWChem input and output file names based on a given SMILES string.

    This function:
    1. Checks if the given SMILES string contains any special characters.
    2. If special characters are found, the SMILES string is enclosed in double quotes.
    3. Generates NWChem input and output file names using the (possibly modified) SMILES string.

    Args:
        smiles (str): A SMILES string representing a molecule.

    Returns:
        tuple: A tuple containing the names of the NWChem input file and output file.

    Global Variables:
        output_file_prefix (str): The prefix for the NWChem input files.
        output_file_name (str): The name of the NWChem output file.    atomic_number_to_symbol = {


    Notes:
        - Special characters in file names can cause issues in shell commands. Enclosing the SMILES string in 
          double quotes helps prevent these issues.
        - The special characters checked for are: @, _, !, #, $, %, ^, &, *, (, ), <, >, ?, /, \, |, }, {, ~, and :.
    """
    # Make own character set and pass this as argument in compile method
    regex = re.compile(r'[@_!#$%^&*()<>?/\|}{~:]')

    # Check for special characters in the string
    if regex.search(smiles) is not None:
        smiles = f'"{smiles}"'

    # Generate NWChem input and output file names
    output_file_prefix = f"{smiles}"
    output_file_name = f"{smiles}.out"
    
    return output_file_prefix, output_file_name
    
    
def free_energy_output(output_file):
 
 
    with open(output_file, 'r') as f:
        content = f.read()
        #  Extract the final DFT energy
        energy_matches = re.findall(r'Total DFT energy =\s+([-]?\d+\.\d+)', content)
        if energy_matches:
            e_elec = float(energy_matches[-1]) * 627.5095
                                          
        # Extract thermal correction to energy
        thermal_corr_match = re.search(r'Thermal correction to Energy\s+=\s+(\d+\.\d+)\s+kcal/mol', content)
        if thermal_corr_match:
            thermal_corr = float(thermal_corr_match.group(1))
                           

        # Extract total entropy
        entropy_match = re.search(r'Total Entropy\s+=\s+(\d+\.\d+)\s+cal/mol-K', content)
        if entropy_match:
            entropy = float(entropy_match.group(1)) / 1000  # Convert to kcal/mol-K
                                                        
        # Calculate free energy
        T = 298.15  # Temperature in K
        G = e_elec + thermal_corr - T * entropy

    return G

def compute_boltzmann_weights(relative_energies, temperature):
    """
    Compute Boltzmann weights based on relative free energies.

    Args:
    - relative_energies (list or numpy array): List of relative free energies for each conformer.
    - temperature (float): Temperature in Kelvin.

    Returns:
    - numpy array: Boltzmann weights for each conformer.
    """
    # Normalize relative energies to prevent overflow
    min_energy = np.min(relative_energies)
    normalized_energies = relative_energies - min_energy

    beta = 1 / (0.0083145 * temperature)  # Boltzmann constant in kJ/mol/K
    exp_terms = np.exp(-beta * normalized_energies)
    boltz_weights = exp_terms / np.sum(exp_terms)
    return boltz_weights

def apply_boltzmann_weights_to_tensors(tensors, boltz_weights):
    """
    Apply Boltzmann weights to NMR shielding tensors.

    Args:
    - tensors (dict): Dictionary containing NMR shielding tensors for each atom.
                      Example: {'H1': np.array([...]), ...}
    - boltz_weights (numpy array): Boltzmann weights for each conformer.

    Returns:
    - dict: Boltzmann-weighted average NMR shielding tensors for each atom.
    """
    weighted_tensors = {}

    for atom, tensor_data in tensors.items():
        weighted_tensor = np.zeros(len(tensor_data[0]))   


        for i in range(len(boltz_weights)):
            weighted_tensor += boltz_weights[i] * tensor_data[i]

        weighted_tensors[atom] = weighted_tensor

    return weighted_tensors
    
    
def parse_shielding_tensor_from_output(output_file):
    """
    Parse isotropic shielding tensor values from NWChem output file.

    Args:
        output_file (str): Path to the NWChem output file.

    Returns:
        dict: A dictionary with combined (atom name + atom index) as keys and isotropic shielding values as values.
    """
    shielding_values = {}

    with open(output_file, 'r') as file:
        data = file.read()

        # Find all matches of atom information using regular expressions
        matches = re.findall(r'Atom:\s+(\d+)\s+(\w+)\s+.*?isotropic\s+=\s+([\d.-]+)', data, re.DOTALL)

        # Iterate over the matches and store the isotropic values
        for match in matches:
            atom_number = int(match[0])
            atom_name = match[1]
            isotropic_value = float(match[2])
            key = f"{atom_name}{atom_number}"
            shielding_values[key] = isotropic_value

    return shielding_values

    

def write_chemical_shifts_to_excel(output_dir, output_excel):
    """
    Write the boltzmann weighted average NMR shielding tensor for each conformer to an Excel sheet.

    Args:
        output_dir (str): Directory containing the NWChem output files.
        output_excel (str): Path to the output Excel file.

    Returns:
        None
    """
    data = []
    free_energies = {}
    tensors = {}

    # Read NWChem output files and extract free energies and shielding values
    for output_file in os.listdir(output_dir):
        if output_file.endswith('.out'):
            conformer_id = os.path.splitext(output_file)[0]
            output_path = os.path.join(output_dir, output_file)

            free_energies[conformer_id] = free_energy_output(output_path)
            shielding_values = parse_shielding_tensor_from_output(output_path)

            if conformer_id not in tensors:
                tensors[conformer_id] = {}

            for atom, shielding in shielding_values.items():
                if atom in tensors[conformer_id]:
                    tensors[conformer_id][atom].append(shielding)
                else:
                    tensors[conformer_id][atom] = [shielding] 

    # Calculate Boltzmann weights
    relative_energies = np.array(list(free_energies.values()))
    boltz_weights = compute_boltzmann_weights(relative_energies, temperature=298.15)

    # Combine tensors for Boltzmann-weighted average
    combined_tensors = {}
    unique_atoms = set(atom for conformer in tensors.values() for atom in conformer.keys())
    for atom in unique_atoms:
        atom_tensors = np.array([tensors[conformer][atom] for conformer in tensors if atom in tensors[conformer]])
        combined_tensors[atom] = apply_boltzmann_weights_to_tensors({atom: atom_tensors}, boltz_weights)[atom]

    print(f"Boltzmann-weighted average NMR shielding tensors: {combined_tensors}")

    

    for atom, shielding_tensor in combined_tensors.items():
        # Convert numpy.ndarray to float if it contains a single value
        if isinstance(shielding_tensor, np.ndarray) and shielding_tensor.size == 1:
            shielding_tensor = shielding_tensor.item()  # or shielding_tensor = shielding_tensor[0]
        
        data.append([atom,shielding_tensor])

    # Write results to Excel
    df = pd.DataFrame(data, columns=['Atom', 'Shielding Tensor'])
    df.to_excel(output_excel, index=False)
    
def main():


    """
    Main function to generate conformers, perform NWChem calculations, and write chemical shifts to an Excel file.

    Steps:
    1. Parse command-line arguments to obtain SMILES string, exchange-correlation functional, case, number of conformers, and output Excel file name.
    2. Generate conformers for the given SMILES string and print their energies for verification.
    3. Write SDF data to a file if needed.
    4. Generate 3D coordinates and create NWChem input files for each conformer.
    5. Run NWChem calculations for each input file using MPI execution.
    6. Write the calculated chemical shifts to an Excel file.
    """

    parser = argparse.ArgumentParser(description="Generate conformers and calculate chemical shifts.")
    parser.add_argument("smiles", type=str, help="SMILES string of the molecule.")
    parser.add_argument("xc_functional", type=str, help="Exchange-correlation functional for NWChem calculations.")
    parser.add_argument("--case", type=int, help="Case 1: 1H Spectrum , Case 2: 13C Spectrum")
    parser.add_argument("--num_conformers", type=int, default=10, help="Number of conformers to generate (default: 10).")
    parser.add_argument("--output_excel", type=str, default="chemical_shifts.xlsx", help="Output Excel file (default: chemical_shifts.xlsx).")

    args = parser.parse_args()

    # Generate conformers
    sdf_data, energies = generate_conformers(args.smiles, args.num_conformers)

    # Print energies for verification
    print("Conformer Energies:")
    for cid, energy in energies:
        print(f"Conformer {cid}: Energy = {energy:.2f}")

    # Example: Write sdf_data to a file if needed
    with open('conformers.sdf', 'w') as f:
        f.write(sdf_data)

    # Generate 3D coordinates and create NWChem input files
    nwchem_input_files = generate_3d_coordinates(args.smiles, sdf_data, args.xc_functional, args.smiles, args.case)
    
    for nwchem_input_file in nwchem_input_files:
        nwchem_input_file, nwchem_output_file = run(nwchem_input_file)
        nwchem_command = f"nwchem {nwchem_input_file} > {nwchem_output_file}"
        final_command = f"mpiexec --use-hwthread-cpus -np 40 {nwchem_command}"
        run_nwchem(final_command)

    output_dir = './'  # Directory containing NWChem output files

    # Write chemical shifts to Excel
    write_chemical_shifts_to_excel(output_dir, args.output_excel)


if __name__ == "__main__":
    main()
