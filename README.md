# AutoNMR-An-Automated-Tool-for-NMR-Chemical-Shift-Calculations-using-NWChem
## Introduction
AutoNMR is a novel Python package designed to streamline and automate the workflow of NMR chemical shift calculations using NWChem. By accepting a molecular structure in the form of a SMILES code, AutoNMR automates the processes of structure generation, conformer generation, geometry optimization, frequency calculations, free energy calculations, and NMR shielding tensor calculations. This package also incorporates Boltzmann weighting to provide accurate chemical shift predictions.
## Methodology
AutoNMR Tool performs all the calclutions using NWChem. Below is the sequence wise algorithm which has been implemented in the tool to give the better results on the automation process.
### Structure & Conformer Generation

AutoNMR uses RDKit to convert the SMILES code into a 3D molecular structure. Conformers are generated using RDKit's ETKDG method, which efficiently samples the conformational space of the molecule. Multiple conformers are generated to ensure comprehensive sampling.

### Generates NWChem Input File

NWChem input file (.nw) is generated for each conformer for geometry optimization, frequency calculations, and shielding tensor calculations using DFT. The sample input file can be seen in the GitHub repository for NMR chemical shift calculations.

### NMR Shielding Tensor Calculation

The NWChem input file for each conformer runs on the terminal for generating the output file for each conformer, which contains information about the shielding tensor values and thermodynamic values useful for further calculations. Necessary information is extracted from the output files of each conformer.

### Free Energy Calculation

The free energy of each conformer is calculated from the vibrational frequencies and the electronic energy obtained from geometry optimization. These free energies are essential for Boltzmann weighting. These free energy values for each conformer are used to determine the Boltzmann weight.

### Boltzmann Weighting

Boltzmann weighting is applied to the calculated shielding tensors based on the free energy of each conformer. This step accounts for the population distribution of different conformers. This concept results in the Boltzmann weighted shielding tensor value for each atom of the molecule.

### Chemical Shift Calculation

The Boltzmann weighted shielding tensor values are averaged to get the cluster midpoint values of hydrogen within the same environment. These values are then used with the Linear Regression Model to get the exact chemical shift value.

### Scaling using Linear Regression

Running the test on molecules and comparing it with the experimental known chemical shift values builds a linear regression model with the same functionals and DFT specifications. Use the intercept and the slope obtained to scale the computational calculated value of shielding tensors.

### NMR Spectrum

The scaled values of chemical shifts are then extracted to plot the NMR spectrum using the Lorentzian function. The concept of J-Coupling and Pascal’s ratio for the intensity of peaks are also included to determine a relevant spectrum.

## Requirements
The following are the requirements to run the AutoNMR package:

- **Python**: Ensure you have Python installed. AutoNMR is compatible with Python 3.6 and above.
- **NWChem**: NWChem must be installed and properly configured on your system to perform quantum chemical calculations.
- **Python Distribution**: The following Python distribuiton are required:
  - **Anaconda 4.3**: It is recommended to use Anaconda for managing the Python environment and dependencies. Anaconda includes many of the necessary packages and makes installation easier.

  Additional packages that need to be installed can be found in the `requirements.txt` file provided in the repository.
