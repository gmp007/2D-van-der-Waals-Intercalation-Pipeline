This directory contains 

1. The optimized structures (CONTCARs) of 24 MX2 bilayers (see "MX2-Bilayers")
2. The POSCARs and CONTCARs of around 1300 hybrid materials studied during active learning are in "CalculatedHybirds.tar"
3. The POSCARs and CONTCARs of the top 50 hybrid materials (result of active learning)

To generate structures of organic molecules in lone-gas phase follow the steps below: 

1. Un-comment lines 50 and 51, and comment-out lines below 51 in the 'run.py' script 

2. Upadate "molecules.txt" with the PubChem IDs and SMILES of the molecules. This information is available in "Active-Learning-Workflow/DataSource/filtered_planar_mols.txt" or "Active-Learning-Workflow/DataSource/filtered_planar_mol_by_fit.txt" or "Active-Learning-Workflow/DataSource/planar_mols_output_final.txt"
