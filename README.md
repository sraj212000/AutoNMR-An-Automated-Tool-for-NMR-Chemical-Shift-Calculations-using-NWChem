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

### Operating System
AutoNMR has been tested on the following operating systems:
- **Linux**: Preferred due to better support for NWChem and scientific computing.
- **macOS**: Compatible, but ensure NWChem is properly configured.
- **Windows**: Can be used with appropriate configuration, but might require additional setup for NWChem.
### Software
- **Python**: Ensure you have Python installed. AutoNMR is compatible with Python 3.6 and above.
- **NWChem**: NWChem must be installed and properly configured on your system to perform quantum chemical calculations.
- **Python Distribution**: The following Python distribuiton are required:
  - **Anaconda 4.3**: It is recommended to use Anaconda for managing the Python environment and dependencies. Anaconda includes many of the necessary packages and makes installation easier.

  Additional packages that need to be installed can be found in the `requirements.txt` file provided in the repository.
## Installation Guide for NWChem
The following tutorial provides instructions for installing NWChem on Linux. If you are using a different operating system, such as Windows or macOS, please visit the [NWChem website](https://nwchemgit.github.io/Download.html#nwchem-availability-in-linux-distributions) for specific instructions tailored to your OS.
### Installation Steps
1. **Open the terminal.**

2. **Update your package list by running the following command:**

    ```sh
    sudo apt update
    ```

3. **Install NWChem by running the following command:**

    ```sh
    sudo apt install nwchem
    ```

4. **Verify the installation by running the following command to check the NWChem command is available in your system's PATH**

    ```sh
    command -v nwchem
    ```

    If NWChem is installed and properly configured, these commands will return the path to the nwchem executable. If NWChem is not installed or not in your PATH, these commands will not return any output.

Congratulations! You have successfully installed NWChem on your system.
## Usage of AutoNMR - A Case Study of 2-Chloropropene
### Brief Introduction
Let's begin to understand how to use AutoNMR tool step-by-step for getting the NMR Chemical Shifts and the NMR Spectrum of the respective atoms.Here we have taken the example of **2-Chloropropene** whose canonical smiles code is `CC(=C)Cl`

The AutoNMR tool contains four important python script which can be downlaoded from this repository. The scripts are mentioned below:

1. [**AutoNMR.py**](AutoNMR.py)- This script performs all the task related to NMR Shift calculation using NWChem. This returns the boltzmann calculated shielding tensor values of atoms for further study. Further details regarding the script can be found later in this tutorial.
2. [**Linear_Regression.py**](Linear_Regression.py) - This script perfoms linear regression calclution on a set of data ( Experimental and calculated ) to build a model which returns the slope and intercept required to scale the shielding tensor values obtained from [**AutoNMR.py**](AutoNMR.py).
3. [**Scaling.py**](Scaling.py)- This script gives the final calculated NMR Chemcial shift values after scaling the resulted output of shielding tensor values obtained from [**AutoNMR.py**](AutoNMR.py) using the slope and the intercept values obtained from [**Linear_Regression.py**](Linear_Regression.py).
4. [**NMR_Simulation.py**](NMR_Simulation.py) - This script gives you the NMR Spectrum of the resultant chemical shift obtained from the [**Scaling.py**](Scaling.py).

After Downloading these script create a working directory where you can save all these script together, say `AutoNMR` which contains all the four script mentioned above.

### Beiginning with [**AutoNMR.py**](AutoNMR.py) - First Step
This script performs the following task using NWChem module : 
- Structure and Conformer Generation
- Generates NWChem Input File
- NMR Shielding Tensor Calculation
- Free Energy Calculation
- Boltzmann Weighting
- Returns the Boltzmann weighted shielding tensor values

**Step 1:** Navigate the terminal to the working directory say `AutoNMR` which contains all these files.

**Step 2:** Get the canonical smiles of the molecule here the example is of **2-Chloropropene** whose canonical smiles code is `CC(=C)Cl`

**Step 3:** Run the following command on the terminal.
```
python3 AutoNMR.py "CC(=C)Cl" "B3LYP" --num_conformers 10 --case 1
```
Now here the smiles code is take for example is `CC(=C)Cl` and the DFT level functional used in the calculation is `B3LYP` you may also use other DFT level functionals as available in the [NWChem Documentation](https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html) 

The `num_conformers` which describes the number of conformers the user wants to consider for studying the resepective molecule whose default value is 10.

The last command line argument `case` explains that which spectrum user wants to study in the respective molecule. It is **1 for Hydrogen spectrum** and **6 for carbon spectrum** so it calculates the shielding tensor values of the specified atoms in the molecule.

**Note** Here the calcultion are performed on No. of processors = 40 . In case your machine is compatible with less or more number of processor you need to mention the number of processors in the line
