# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:39:09 2022
The purpose of this code is to find potential intercalates
by screening the PubChem database.
@author: chr218
"""

import numpy as np
import sys
from rdkit import Chem as CH
from rdkit.Chem import Descriptors, Descriptors3D
from rdkit.Chem import AllChem


'''
constants & thresholds
'''

planarity_threshold = 1.0e-5
#Cyclohexane PBF = 0.293, sphl=0.125856
#PBF of toluene is: 1.16E-06

'''
Functions and Subroutines
'''

def get_smiles():
    smi_list = []
    f=open('Pubchem_atoms_3_final.txt', "r")
    for line in f:
        strip_line = line.strip()
        smi_list.append(strip_line.split())
    f.close()
    return smi_list

def filter_planar(smi):
    f = open('planar_mols.txt',"w")
    planar_mol_list = []
    for row in smi:
        smi_str = row[-1]
        CID = row[0]
        try:
            m = CH.MolFromSmiles(smi_str)#convert smiles string to Rdkit object.
            mh = CH.AddHs(m)#Add hydrogens.
            AllChem.EmbedMolecule(mh)#Get 3D coordinates.
            Num_Rings = CH.rdMolDescriptors.CalcNumRings(mh)#Calculate numebr of Rings in molecule.
            if Num_Rings > 0:
                try:
                    PBF = CH.rdMolDescriptors.CalcPBF(mh)#plane of best fit,or reference, the PBF of toluene is:1.16E-06
                except:
                    PBF = 1
                #sphI = Descriptors3D.SpherocityIndex(mh)#Additionally, we can look at "sphereocity", toluene:2.18E-12
                if PBF < planarity_threshold:
                    planar_mol_list.append(smi)
                    f.write('%s %s\n' % (str(CID),str(smi_str)))
        except:
            pass
            
    f.close()
    return planar_mol_list
            
            
        

'''
Program Begins Here
'''

smi_list = get_smiles()
planar_list=filter_planar(smi_list)

#Save smiles to text file
#f = open('planar_mols.txt')
#for smi in planar_list:
#    f.write('%s\n' % (str(smi)))
#f.close()


