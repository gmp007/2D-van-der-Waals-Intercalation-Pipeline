# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 12:01:15 2020

@author: chr218
"""

import sys#For break points.
import numpy as np
import subprocess#For executing BASH commands.
import openbabel#For cif manip

from openbabel import pybel as pb#For SMILES 2 mole coord
from ase.io import read, write#For reading xyz
from ase.io import vasp, cif#For reading VASP & cif.
from ase import Atoms
from ase.build import mx2#For structure db
from ase.constraints import FixAtoms, FixCom, FixedPlane
from ase.visualize import view
from scipy import constants,optimize


from rdkit import Chem as CH
from rdkit.Chem import Descriptors, Descriptors3D
from rdkit.Chem import AllChem

import prep_VASP
import gen_strucs

if __name__ == "__main__":

    #collect molecule SMILES
    SMILES_list = []
    mol_names = []
    f=open('molecules.txt')
    for line in f:
        split_line = line.strip().split()
        SMILES_list.append(split_line[1])
        mol_names.append(split_line[0])
    f.close()

    #collect Bilayer formula
    Bilayer_list = []
    f=open('Bilayers.txt')
    for line in f:
        Bilayer_list.append(line.strip())
    f.close()

    #loop over combinations of SMILES and Bilayers, generating dirs.
    #for i in range(len(SMILES_list)):
        #prep_VASP.make_gas_DFT(SMILES_list[i],mol_names[i])

    #for j in range(len(Bilayer_list)):
    #    prep_VASP.make_Bulk(Bilayer_list[j])
    #    prep_VASP.make_Bilayer(Bilayer_list[j])
 

 #uncomment and run after gas & frmwrk relaxed. MAKE SURE TO COMMENT-OUT THE OTHER TWO LOOPS!!!
    for i in range(len(SMILES_list)):
        for j in range(len(Bilayer_list)):
            gen_strucs.insert_mol(mol_names[i],Bilayer_list[j])
