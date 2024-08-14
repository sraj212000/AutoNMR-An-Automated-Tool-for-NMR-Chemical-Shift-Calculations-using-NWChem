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
#from NMR_Simulation import get_hydrogen_neighbors,nmr_splitting_patterns,plot_nmr_spectrum,extract_chemical_shifts_from_excel

    
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

    
def compute_chemical_shifts(shielding_values, reference_shielding):
    """
    Compute chemical shifts by subtracting the reference shielding tensor value.

    Args:
        shielding_values (dict): Dictionary of isotropic shielding values for a conformer.
        reference_shielding (dict): Dictionary of reference shielding tensor values.

    Returns:
        dict: Chemical shifts for each atom.
    """
    chemical_shifts = {}
    for atom, shielding in shielding_values.items():
        atom_name = ''.join([i for i in atom if not i.isdigit()])
        if atom_name in reference_shielding:
            chemical_shifts[atom] =  shielding - reference_shielding[atom_name] 
        else:
            chemical_shifts[atom] = None  # No reference value available

    return chemical_shifts

def write_chemical_shifts_to_excel(output_dir, reference_shielding, output_excel):
    """
    Write the chemical shifts for each conformer to an Excel sheet.

    Args:
        output_dir (str): Directory containing the NWChem output files.
        reference_shielding (dict): Dictionary of reference shielding tensor values.
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

    # Compute chemical shifts using the combined Boltzmann-weighted average shielding tensor values
    chemical_shifts = compute_chemical_shifts(combined_tensors, reference_shielding)

    for atom, shielding_tensor in combined_tensors.items():
        # Convert numpy.ndarray to float if it contains a single value
        if isinstance(shielding_tensor, np.ndarray) and shielding_tensor.size == 1:
            shielding_tensor = shielding_tensor.item()  # or shielding_tensor = shielding_tensor[0]
        
        # Compute the corresponding chemical shift
        shift = chemical_shifts.get(atom, None)
        if isinstance(shift, np.ndarray) and shift.size == 1:
            shift = shift.item()  # or shift = shift[0]
        
        data.append([atom, shift, shielding_tensor])

    # Write results to Excel
    df = pd.DataFrame(data, columns=['Atom', 'Chemical Shift', 'Shielding Tensor'])
    df.to_excel(output_excel, index=False)
    
def main():
   
    # Parse reference shielding into a dictionary
    reference_shielding = {
        'H': 32.8684825,
        'C': 190.9843647
    } 
    print(reference_shielding)
  
    output_dir = './'  # Directory containing NWChem output files

    # Write chemical shifts to Excel
    write_chemical_shifts_to_excel(output_dir, reference_shielding,'chemical_shifts.xlsx')
    #chemical_shifts = extract_chemical_shifts_from_excel(args.output_excel)
    #atom_h_neighbors = get_hydrogen_neighbors(args.smiles, spectrum_type="1H")
    #splitting_patterns = nmr_splitting_patterns(chemical_shifts, atom_h_neighbors)
    #plot_nmr_spectrum(splitting_patterns, spectrum_type="1H")


if __name__ == "__main__":
    main()
