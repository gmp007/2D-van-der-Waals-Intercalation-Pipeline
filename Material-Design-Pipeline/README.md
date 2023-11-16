## Instructions:
How to generate structures and other VASP input files of organic molecules, MX2 bilayers, and intercalated hybrids (MX2/organic)? Follow these steps to generate structures using the provided workflow files:

(If you are here for the pre-generated structures of the systems that were used in the active learning workflow, go inside "Structures-VASP-format".)

### 1. Setup
1. Copy the contents of the "VDW_workflow_files" directory to a local directory on your machine.

2. Edit the `prep_VASP.py` file:
   - Replace the "VDW_PATH" variable within `prep_VASP.py` with the path to the "VDW_workflow_files" directory on your system.

3. Add VASP POTCAR files:
   - Place VASP POTCAR files in "VDW_workflow_files\POTCAR_files\PBE_52" using the naming convention POTCAR_"element-symbol" (ex: POTCAR_H, POTCAR_C). Ensure inclusion of POTCAR files for all investigated atoms.

4. Replace the submission script:
   - Replace the contents of "VDW_workflow_files\other_DFT_files\submit_vasp_gam.sh" with your own HPC submission script.

5. Create Conda Environment and Install Packages:
   - Then create a conda environment, name - VDW_wrkflow (can be anything) - install rkdit:
     ```
     conda create -c conda-forge -n VDW_wrkflow rdkit
     ```
   - Now activate the environment:
     ```
     conda activate VDW_wrkflow
     ```
   - Finally, install the following packages (make sure you wait before each package is installed before executing the subsequent commands. this might take a while):
     1) ```
        conda install numpy
        ```
     2) ```
        conda install -c conda-forge openbabel
        ```
     3) ```
        conda install -c conda-forge ase
        ```
     4) ```
        conda install -c anaconda scipy
        ```
   6. Alternatively, you can replicate the conda environment using the configuration file "VDW_wrkflow.yml"
        ```
        conda env create -f VDW_wrkflow.yml -n your-env-name
        ```

### 2. Execution
Within the `run.py` script:

1. Update the input files:
   - Update the "molecules.txt" file with PubChem IDs and SMILES of molecules for intercalation.
   - Update the "Bilayers.txt" file with the chemical composition of bilayer transition metal dichalcogenides.

2. Generate VASP input files:
   - Loop over combinations of molecules and bilayers:
     - Execute lines 50-54 in the script to generate directories with VASP input files.
     - Do not execute lines 55, 59-61 at this stage. Bilayer formation should occur only upon bulk convergence.

3. Bilayer formation:
   - Upon convergence of bulk structures:
     - Execute line 55 in place of line 54 in the script to create bilayer structures.

4. Intercalated structures:
   - Upon convergence of isolated gas molecules and bilayers:
     - Execute lines 59-61 in the script to generate intercalated structures.
     - Do not re-execute lines 50-56.



