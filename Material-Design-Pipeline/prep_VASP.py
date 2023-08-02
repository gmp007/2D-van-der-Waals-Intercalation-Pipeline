'''
This file includes all functions related to
making directories for running VASP calculations.
'''

import sys#For break points.
import os #Get envir vars
import numpy as np
import subprocess#For executing BASH commands.

from ase.io import read, write#For reading xyz
from ase.io import vasp, cif#For reading VASP & cif.
from ase import Atoms, formula
from ase.build import mx2#For structure db

from rdkit import Chem as CH
from rdkit.Chem import Descriptors, Descriptors3D
from rdkit.Chem import AllChem

import gen_strucs 

#Get VDW environ var
echo_arg = os.environ["USER"]
cmd = subprocess.Popen(["echo", echo_arg], stdout=subprocess.PIPE)
VDW_PATH='/home/'+cmd.stdout.read().rstrip().decode("utf-8")+'/bin/VDW_workflow_files'

'''
Make POTCAR for POSCAR

Input: ASE atoms object.

'''
def make_POTCAR(atoms):
        symbol_list = []
        for atom in atoms:
            symbol_list.append(atom.symbol)

        symbol_arr = np.asarray(symbol_list)
        POTCAR_symbols = np.unique(symbol_arr)

        print('Generating POTCAR for atoms: \n')
        for sym in POTCAR_symbols:
                print(sym)
                subprocess.call(['cat '+VDW_PATH+'/POTCAR_files/PBE_52/POTCAR_'+sym+' >> POTCAR'],shell=True)

'''
Collects INCAR, submission script and KPOINTS to
directory.

Input: ASE atoms object

'''
def prep_DFT_files(atms_):
    make_POTCAR(atms_)#Generate POTCAR file
    subprocess.call("cp "+VDW_PATH+"/other_DFT_files/INCAR .",shell=True)#Generate INCAR file. (Just pull it from bin for now- FIX LATER)
    subprocess.call("cp "+VDW_PATH+"/other_DFT_files/submit_vasp_gam.sh .",shell=True)#Generate submission script, again just pull from bin for now.
    subprocess.call("cp "+VDW_PATH+"/other_DFT_files/KPOINTS .",shell=True)#Generate KPOINTS file, again just pull from bin for now.


'''
Makes directories and lone gas-phase intercalant structures
with necessary VASP files.

Input: name of directory, unit cell parameters. (default=None)

'''
def make_gas_DFT(smiles,nme,ucell=None):

    #define directory name
    dir_name = 'gas_'+nme

    #Make 3D molecule from SMILES using Rdkit 
    gen_strucs.SMILES2coord(smiles)

    #Assume Rdkit output will always be 'mol_3D.xyz'
    atoms_=read('mol_3D.xyz')

    #get moment of inertia matrix.
    I = get_PMI_vec(atoms_)

    if ucell == None:
        atoms_.center(vacuum=15.0,about=(0,0,0))#Add vaccum along each side & center mol.
    else:
        atoms_.cell[:] = ucell#Set unit cell to input array. Note: must be 2D array!

    I_inv = np.linalg.inv(I)#calculate inverse of moment of inertia matrix.
    atoms_.set_positions(np.dot(atoms_.get_positions(),I_inv))#tranform coordinates of atoms to princple axes.

    ucell_vec = atoms_.cell.cellpar()[:3]
    new_pos = np.zeros(np.shape(atoms_.get_positions()))
    for i,atm in enumerate(atoms_):
        for j,pos in enumerate(atm.position):
            new_pos[i,j] += pos + ucell_vec[j]/2

    atoms_.set_positions(new_pos)
    vasp.write_vasp('POSCAR',atoms_,sort=True)
    vasp.read_vasp('POSCAR')
    prep_DFT_files(atoms_)
    subprocess.call("mkdir "+dir_name,shell=True)
    subprocess.call("mv KPOINTS POTCAR submit_vasp_gam.sh INCAR POSCAR "+dir_name,shell=True)

'''
Calculates Principle Moments of Inertia of lone gas-phase intercalant structure.
Necessary for orienting gas-phase structure such that its plane faces along z-axis.

Input: ASE atoms object

Return: Inertia Tensor
'''
def get_PMI_vec(atms_):
    #translate mol to COM
    com = atms_.get_center_of_mass(scaled=False) 
    for atom_ in atms_: 
        atom_.position = atom_.position - com

    #caclulate moments of inertia tensor 
    eig_val, eig_vec = atms_.get_moments_of_inertia(vectors=True)

    return eig_vec

'''
Makes directories and bulk structures with necessary
VASP files.

Input: Bulk structure formula as string., phase (defualt 2H for now-make sure directory is within bin/$VDW_WORKFLOW_PATH)

'''
def make_Bulk(frmlu,phase_='2H'):

    #Define dir
    dir_name = 'struc_Bulk_'+frmlu

    frmlu_obj=formula.Formula(frmlu)#Initialize ASE formula object.
    frmlu_sorted = frmlu.format('metal')#Order formula such that metals come first.
    frmlu_obj=formula.Formula(frmlu)#Initialize ASE formula object.
    frmlu_dict=frmlu_obj.count()
    metal,chalcogen = list(frmlu_dict.keys())

    #Get reference MoS2 framework.
    frmwrk_ref=vasp.read_vasp(VDW_PATH+'/ref_frmwrks/'+phase_+'/MoS2/POSCAR')
     
    #Change frmwrk symbols to match frmlu
    sym_arr=frmwrk_ref.get_chemical_symbols()
    new_sym = []
    for i,sym in enumerate(sym_arr):
        if sym == 'Mo':
            new_sym.append(metal)
        elif sym =='S':
            new_sym.append(chalcogen)

    #Set framework symbols.
    frmwrk_ref.set_chemical_symbols(new_sym)
    frmwrk = frmwrk_ref.copy()

    #Prep DFT files
    vasp.write_vasp('POSCAR',frmwrk,sort=True)
    prep_DFT_files(frmwrk)
    subprocess.call("mv INCAR "+"INCAR_org",shell=True)
    out_file = open('INCAR', "w")
    sub = subprocess.call(['sed', 's/.*ISIF.*/ISIF=3/g', 'INCAR_org'], stdout=out_file)
    subprocess.call("rm INCAR_org",shell=True)
    subprocess.call("mkdir "+dir_name,shell=True)
    subprocess.call("mv KPOINTS POTCAR submit_vasp_gam.sh INCAR POSCAR "+dir_name,shell=True)

'''
Makes directories and Bilayer structures with necessary
VASP files.

Input: Bilayer structure formula as string., phase (defualt 2H for now-make sure directory is within bin/$VDW_WORKFLOW_PATH)


'''
def make_Bilayer(frmlu):

    #Define dir to extract CONTCAR
    dir_name = 'struc_Bulk_'+frmlu

    #Collect relaxed Bulk CONTCAR.
    frmwrk = vasp.read_vasp(dir_name+'/CONTCAR')
    ucv = frmwrk.get_cell()#Extract "unit cell vector"
    ucv[2][-1] += 15#Expand unit cell vector in z-direction to add vaccum between bilayers.
    frmwrk.set_cell(ucv)#Set unit cell vector.

    #Define dir to write Bilayer POSCAR
    dir_name = 'struc_Bilayer_'+frmlu

    #Prep DFT files
    vasp.write_vasp('POSCAR',frmwrk,sort=True)
    prep_DFT_files(frmwrk)
    subprocess.call("mv INCAR "+"INCAR_org",shell=True)
    out_file = open('INCAR', "w")
    sub = subprocess.call(['sed', 's/.*ISIF.*/ISIF=2/g', 'INCAR_org'], stdout=out_file)
    subprocess.call("rm INCAR_org",shell=True)
    subprocess.call("mkdir "+dir_name,shell=True)
    subprocess.call("mv KPOINTS POTCAR submit_vasp_gam.sh INCAR POSCAR "+dir_name,shell=True)


if __name__ == '__main__':

    #Generate DFT relaxation of gas-phase thiophene from SMILEs string.
    smiles_ = 'c1ccsc1'
    name_ = 'thiophene'
    make_gas_DFT(smiles_,name_)

    #Generate of DFT relaxation of Bulk WTe2.
    make_Bulk('WTe2')

    #Generate Bilayer of Bulk strucutre.
    make_Bilayer('WTe2')

 

