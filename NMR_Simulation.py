import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
import pandas as pd
import argparse

def get_hydrogen_neighbors(smiles, spectrum_type):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    
    atom_h_neighbors = {}
    atom_index=0
    for atom in mol.GetAtoms():
        if spectrum_type == "1H" and atom.GetAtomicNum() == 1:
            atom_index = atom_index + 1
            
            # Get the carbon atom this hydrogen is attached to
            carbon_atom = next((nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6), None)
            
            neighboring_hydrogens = 0
            equivalent_hydrogens = 0
        
            if carbon_atom is not None:
                    # Count the hydrogens attached to the same carbon atom (equivalent hydrogens)
                    equivalent_hydrogens = sum(1 for nbr in carbon_atom.GetNeighbors() if nbr.GetAtomicNum() == 1)
                
                    # Count the hydrogens attached to neighboring carbon atoms (for splitting pattern)
                    
                    for neighbor in carbon_atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 6:  # Ensure it is a carbon neighbor
                            neighboring_hydrogens += sum(1 for nbr in neighbor.GetNeighbors() if nbr.GetAtomicNum() == 1)
            
            atom_h_neighbors[atom_index] = (equivalent_hydrogens, neighboring_hydrogens)

    print(atom_h_neighbors)
    return atom_h_neighbors

def nmr_splitting_patterns(chemical_shifts, atom_h_neighbors, J_coupling=7.0, spectrometer_frequency=90000000):
    splitting_patterns = {}

    for atom_index, (equivalent_hydrogens, neighboring_hydrogens) in atom_h_neighbors.items():
        atom_label = f"H{atom_index}"
        shift = chemical_shifts.get(atom_label, None)
        if shift is not None:
            num_peaks = neighboring_hydrogens + 1
            coupling_ppm = J_coupling / spectrometer_frequency
            splitting_pattern = (num_peaks, coupling_ppm * equivalent_hydrogens, equivalent_hydrogens)
            splitting_patterns[atom_label] = (shift, splitting_pattern)

    print(splitting_patterns)
    return splitting_patterns
    

def lorentzian(x, x0, gamma):
    return (gamma / np.pi) / ((x - x0) ** 2 + gamma ** 2)

def pascals_triangle_row(n):
    row = [1]
    for k in range(1, n + 1):
        row.append(row[k - 1] * (n - k + 1) // k)
    return row

def plot_nmr_spectrum(splitting_patterns, spectrum_type, line_width=0.02, points=10000000):
    fig, ax = plt.subplots()

    x_min = min(shift - coupling / 2 for shift, (num_peaks, coupling, _) in splitting_patterns.values()) - 1
    x_max = max(shift + coupling / 2 for shift, (num_peaks, coupling, _) in splitting_patterns.values()) + 1
    x = np.linspace(x_min, x_max, points)
    y = np.zeros_like(x)

    for atom_label, (shift, (num_peaks, coupling, intensity)) in splitting_patterns.items():
        peak_positions = np.linspace(shift - coupling / 2, shift + coupling / 2, num_peaks)
        intensities = pascals_triangle_row(num_peaks - 1)
        total_intensity = sum(intensities)
       
        for peak_position, relative_intensity in zip(peak_positions, intensities):
            y +=  intensity * lorentzian(x, peak_position, line_width) 
           
   
    ax.plot(x, y, color='blue')
    print (y)
    ax.set_xlabel('Chemical Shift (ppm)')
    ax.set_ylabel('Intensity')
    ax.invert_xaxis()
    plt.title(f'Simulated {spectrum_type} NMR Spectrum')
    plt.show()


def extract_chemical_shifts_from_excel(file_path):
    # Read the Excel file
    df = pd.read_excel(file_path)
    # Drop rows with missing values
    df = df.dropna()
    
    # Initialize a dictionary to store chemical shifts
    chemical_shifts = {}
    
    # Iterate over the rows of the dataframe
    for _, row in df.iterrows():
        # Get the first string before the space as the atom label
        atom_label = row['Atoms'].split()[0]
        # Get the shift value
        shift = row['Shift Values']
        # Store the shift value in the dictionary
        chemical_shifts[atom_label] = shift
    
    # Print the dictionary of chemical shifts
    print(chemical_shifts)
    
    # Return the dictionary
    return chemical_shifts

# Example usage:
# file_path = 'path_to_your_excel_file.xlsx'
# chemical_shifts = extract_chemical_shifts_from_excel(file_path)




def main():
	
    parser = argparse.ArgumentParser(description="Generate NMR Spectrum")
    parser.add_argument("smiles", type=str, help="SMILES string of the molecule.")
    args = parser.parse_args()
   
        
    chemical_shifts = extract_chemical_shifts_from_excel('chemical_shifts.xlsx')
    atom_h_neighbors = get_hydrogen_neighbors(args.smiles, spectrum_type="1H")
    splitting_patterns = nmr_splitting_patterns(chemical_shifts, atom_h_neighbors)
    plot_nmr_spectrum(splitting_patterns, spectrum_type="1H")

if __name__ == "__main__":
    main()
