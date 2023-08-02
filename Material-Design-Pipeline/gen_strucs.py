'''
This file includes all functions related 
to the manipulation of atoms objects. This
includes but is not limited to the seperation
of bilayers, and inserting intercalants into the
bilayers.
'''

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

'''

Input: ASE atoms object of: DFT converged intercalant, DFT converged seperated Bilayer, Ceneter of mass vector of the Bilayer.

Return: ASE atoms object of constrained, intercalated, bilayer.

'''

def apply_constraints(mol_atms_tmp,mol_frm_tmp,frmwrk_com):

    #Sort indices between intercalant atoms and Bilayer atoms.
    constraint_list = []
    mol_indices = []
    top_layer_indices = []
    bottom_layer_indices = []
   
 
    for atm in mol_frm_tmp:
            if atm.position in mol_atms_tmp.positions:
                mol_indices.append(atm.index)
            else:
                if atm.position[-1] > frmwrk_com[-1]:
                    top_layer_indices.append(atm.index)
                else:
                    bottom_layer_indices.append(atm.index) 
   
    #Fix bottom layer atoms.
    bottom_layer = FixAtoms(bottom_layer_indices)  
    constraint_list.append(bottom_layer)
    
    #Set constraints
    mol_frm_tmp.set_constraint(constraint_list)
    
    return mol_frm_tmp

'''
Inserts converged intercalant into Bilayer. User has chouice of inclusing random orientations.

Input: ASE atoms objects of  gas-phase intercalant and Bulk structure, framework name, maximum and minimum ucell vector along z axis,
number of random orientations (defualt is 1), Boolean of whether to perform random rotations and translations, Boolean of whether to fix bottom layer.

'''    
def gen_orient(mol,frm,mol_nme,frmwrk_nme,zmax,zmin,N_images=1,rand_rotate=False,rand_translate=True,set_fixed=True):#turn-on random trans/rot here. Set rand_translate=True etc.
    #Define dir
    dir_name = mol_nme+'_'+frmwrk_nme
    
    i = 0
    num_attempts = 0
    while i < N_images:        
        #np.random.seed(i)

        #Create copy of atoms object 'mol'.
        mol_cp = mol.copy()        

        #Randomly generate rotation about COM axis for gas mol placement.
        if rand_rotate:
            rot_gen = np.multiply(np.random.random_sample([1]),np.array([180]))
    
            #mol_cp.euler_rotate(phi=rot_gen[0],theta=rot_gen[1],psi=rot_gen[2],center=('COM'))
            mol_cp.euler_rotate(psi=rot_gen[0],center=('COM'))

        if rand_translate:
            #Randomly generate center of mass position for gas mol to be placed.
            com_gen_XY = np.random.uniform(0,0,[2,1])
            #com_gen_Z = np.random.uniform(zmin,zmax)
            com_gen_Z = zmin + (zmax-zmin)/2 #place exaclty in center of bilayer.
            com_gen = np.append(com_gen_XY,com_gen_Z)

            #Replicate com_gen vector along axis= 0 to form array.
            com_gen = np.tile(com_gen,(len(mol_cp.get_chemical_symbols()),1))
            #Translate gas mol to new com position and wrap atoms outstanding from unit cell.
            mol_cp.set_scaled_positions(np.add(mol_cp.get_scaled_positions(),com_gen))
            mol_cp.wrap()


        #Create system 'frm_mol' by appending mol atoms object to framework atoms object.
        mol_frm = mol_cp+frm
        mol_frm.wrap()
        mol_atms = mol_frm[0:mol.get_global_number_of_atoms()]

        #Ensure new com position does not overlap with framework atoms.
        atoms_far_enough_criteria = 0
        for j,atom in enumerate(mol_atms):
            dist = mol_frm.get_distances(j,range(mol.get_global_number_of_atoms(),mol_frm.get_global_number_of_atoms()))
            if np.amin(dist) > 0.1:#Adjust this for molecule placement. Helps avoid overlap!
                atoms_far_enough_criteria += 1
             

        #Ensure molecule is within bilayer
        atoms_within_bilayer_criteria = 1
        for j,pos in enumerate(mol_atms.get_scaled_positions(wrap=True)):
            if pos[-1] > zmax or pos[-1] < zmin:
                atoms_within_bilayer_criteria = 0
                break


        if atoms_far_enough_criteria == mol.get_global_number_of_atoms() and atoms_within_bilayer_criteria: 
            if set_fixed:
                mol_frm = apply_constraints(mol_atms,mol_frm,frm.get_center_of_mass())# REMINDER: include options to decided which constraints to implement!

            mol_frm.wrap()
            f_out_dir = (dir_name+'_%s' % str(i))
            f_out = ('POSCAR')
            mol_frm.wrap()
            vasp.write_vasp(f_out,(mol_frm),sort=True)
            #view(mol_frm) #uncomment to view formed CONTCARs within ASE gui.i

            prep_VASP.prep_DFT_files(mol_cp+frm)
            subprocess.call("mkdir "+f_out_dir,shell=True)
            subprocess.call("mv KPOINTS POTCAR submit_vasp_gam.sh INCAR POSCAR "+f_out_dir,shell=True)
            i += 1
            continue

        elif num_attempts > 100:
            print('Insertion failed, check structures:',mol_nme,frmwrk_nme)
            break
        else:
            num_attempts += 1
            continue
        

'''
Converts SMILES string into 3D relaxed structure
and writes to .xyz file for further interpretation.

Input: SMILES string of intercalant.

'''
def SMILES2coord(smi_):

    m = CH.MolFromSmiles(smi_)#convert smiles string to Rdkit object.
    mh = CH.AddHs(m)#Add hydrogens.
    AllChem.EmbedMolecule(mh)#Get 3D coordinates.
    CH.rdmolfiles.MolToXYZFile(mh,'mol_3D.xyz')
    
    #Use openbabel to convert SMILEs to 3D config and write as xyz
    #mol_ = pb.readstring("smi", smi_)#OpenBabel to read input as SMILES string.
    #mol_.draw(show=False, update=True)#Openbabel to draw 2D config from SMILES.
    #mol_.addh()#Must add hydrogen atoms to 2D config. (See man)
    #mol_.write("xyz","mol_2D.xyz",overwrite=True)#Write 2D struc to xyz file.
    #mol_.make3D(forcefield='mmff94',steps=50)#Relax SMILES mol using mmff94 ff for 50 steps to gen 3D config.
    #mol_.write("xyz","mol_3D.xyz",overwrite=True)#Write 3D struc to xyz file.
    #NOTE: OpenBabel VASP capabilites are Read-only!

'''
Inserts relaxed inercalant into relaxed Bilayer.

Input: name of intercalant (string), Bilayer structure formula
'''
def insert_mol(mol_nme,struc_nme):
    #Generate atoms objects for molecule and framework..
    gas_mol = vasp.read_vasp('./gas_'+mol_nme+'/CONTCAR')
    #Get framework atoms object.
    frmwrk = vasp.read_vasp('./struc_Bilayer_'+struc_nme+'/CONTCAR')
    
    #Obtain reference atom type and number of bulk, single layer.
    #frmwrk_ref = vasp.read_vasp('./struc_Bulk_'+struc_nme+'/CONTCAR')
    frmwrk_ref = mx2(formula=struc_nme, a=3.18, thickness=3.19, size=(3, 3, 1),vacuum=1.5)#Also change Bulk dimensions here.-But always keep as (N,N,1). Need number of atoms of monolayer!

    #Seperate Bilayer by iterating through Bulk atoms, and sorting Bilayer atoms by z- distance
    btm_layer_atm_idcs = []
    for i,ref_atm in enumerate(frmwrk_ref):#Loop over Bulk, single, layer as reference.-Want to find atoms in Bilayer which are closest to these!
        dist_list = []
        idx_list = []
        for j, atm in enumerate(frmwrk):#iterate over bilayer atoms
            if ref_atm.symbol == atm.symbol and j not in btm_layer_atm_idcs:
                dist_list.append(atm.position[-1])
                idx_list.append(j)

        btm_layer_atm_idcs.append(idx_list[np.argsort(np.asarray(dist_list))[0]])#find which atom had smallest 'z' distance.

    #Translate upper layer of Bilayer by "x" [ang].
    upper_layer_z_dist = []
    lower_layer_z_dist = []

    upper_layer_z_dist_cart = []
    lower_layer_z_dist_cart = []
    for i,atm in enumerate(frmwrk):
        if i not in btm_layer_atm_idcs:
            frmwrk[i].position[-1] += 3.0#CHANGE THIS VALUE TO CHANGE INTERLAYER DISTANCE.
            upper_layer_z_dist.append(frmwrk.get_scaled_positions(wrap=True)[i][-1])
            upper_layer_z_dist_cart.append(frmwrk.get_positions(wrap=True)[i][-1])
 
        else:
            lower_layer_z_dist.append(frmwrk.get_scaled_positions(wrap=True)[i][-1])
            lower_layer_z_dist_cart.append(frmwrk.get_positions(wrap=True)[i][-1])



    #obtain largest/smallest "z" bilayer distances to randomly place molecule within bilayer.
    z_max = np.sort(np.asarray(upper_layer_z_dist))[0]
    z_min = np.sort(np.asarray(lower_layer_z_dist))[-1]

    #Set molecule's unit cell to be that of framework's unit cell.
    gas_mol.set_cell(frmwrk.cell.cellpar())
    gas_mol.center()

    #Determine angle to rotate the intercalant, such that it's PMI is 
    #aligned paralell to the hexagonal unit-cell diagonal.
    la_,lb_,lc_,alpha_,beta_,gamma_ = gas_mol.cell.cellpar()[:] 
    #print(a_,b_,c_,alpha_,beta_,gamma_)
    a_,b_,c_ = gas_mol.get_cell(complete=True)
    #print(a_,b_,c_)
    #diagonal vector is sum of a, b vectors.
    d_ = -a_ + b_ #We want negative 'a' because 'a' points in positive direction.
    ld_ = np.sqrt(np.sum(np.square(d_)))
    #Normalize unit cell vector.
    d_ /= ld_

    I = prep_VASP.get_PMI_vec(gas_mol)
    I0 = I[0,:]#Get smallest moment of inertia vector. (corresponds to diaganol)
    
    gas_mol.rotate(tuple(I0),tuple(d_),center='COM')#rotate PMI_0 of intercalant into diagnol vector.

    #Ensure molecule periodic images do not overlap:
    gas_mol_tmp = gas_mol.repeat((2,2,1))#Repeat gas molecule after converting unit cell to hexagonal.

    #Calculate minimum distance excluding the molecule's own bonds.
    dist_tmp = []
    dist_inds = []
    dist_pos = []
    for i, mol_atm in enumerate(gas_mol):
        for j in range(gas_mol.get_global_number_of_atoms(),gas_mol_tmp.get_global_number_of_atoms()):
            dist_tmp.append(np.linalg.norm(gas_mol_tmp[i].position-gas_mol_tmp[j].position))
            dist_inds.append([i,gas_mol_tmp[j].index])
            dist_pos.append([gas_mol_tmp[i].position,gas_mol_tmp[j].position])

    dist_arr = np.asarray(dist_tmp,dtype=float)
    min_idx = np.argmin(dist_arr)
    if dist_arr[min_idx] < 1.3:#Set overlap between periodic images threshold here.
        print('Warning, overlap between periodic images. Check Strcuture')
        print('Distance [ang]: ',dist_arr[min_idx])
        print('Atom indices: ',dist_inds[min_idx])
        print('Atom positions [ang]:',dist_pos[min_idx])
        print(mol_nme,struc_nme)
        
    else:
        print('Maximum overlap between periodic images [ang]:',dist_arr[min_idx])
        gas_mol.wrap()
 
        #insert gas-molecule into Bilayer.
        gen_orient(gas_mol,frmwrk,mol_nme,struc_nme,z_max,z_min) 

if __name__ == '__main__':
    #Generate ASE atoms object by reading converged gas-phase intercalant.
    mol_smiles_ = 'c1ccsc1'
    mol_name_ = 'thiophene'
    struc_frmlu_ ='WTe2'

    insert_mol(mol_name_,struc_frmlu_)


