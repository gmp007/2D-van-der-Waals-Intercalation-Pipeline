# -*- coding: utf-8 -*-
"""
November 10 2022

@author: smk520
"""
import numpy as np
import pandas as pd


'''
Functions
'''
'''
The function 'IC_calculator' will calculate the intercalation energy given by the formula E_I = E_system - E_2D - E_gas. 'sys_data' is a .txt file containing two columns, 'System' and 'IC' with name of the system in the form 12606301_CrS2_0, and E_system value. 'gas_data' is a two column .txt file with first column 'Gas' containing gas_CID, e.g, gas_12606301 and second column 'Energy' being the DFT optimized energy of the gas molecule. 'bilayer_data' is another two column .txt file with first column 'Bilayer' containing the name of the structure in the form 'struc_Bilayer_CrO2' and second column with the DFT optimized total energy of the bilayer.  
'''

def IC_calculator(sys_data, gas_data, bilayer_data):
    gasCID = []
    MaterialName = []
    for i in sys_data['System']:
        strip_text = str(i).split('_', 1)
        MatName = strip_text[1].strip('_0')
        gasCID.append(strip_text[0])
        MaterialName.append(MatName)
    
    sys_data['Name'] = pd.DataFrame(MaterialName)
    sys_data['CID'] = pd.DataFrame(gasCID)
    
    gasCID = []
    for i in gas_data['Gas']:
        strip_text = str(i).strip('gas_')
        gasCID.append(strip_text)
    gas_data['CID'] = gasCID
    
    Material = []
    for i in bilayer_data['Bilayer']:
        MatName = str(i).strip('struc_Bilayer_')
        Material.append(MatName)
    
    bilayer_data['Name'] = pd.DataFrame(Material)
    
    res = sys_data.merge(gas_data, on=["CID"], how="right")
    res = res.dropna()
    ICE_data = res.merge(bilayer_data, on=["Name"], how="left")
    
    ICE_data['IC_ENERGY'] = ICE_data['IC'] - ICE_data['Energy'] - ICE_data['BilayerEnergy']
    ICE_data.drop(['Bilayer', 'Gas', 'IC', 'Energy', 'BilayerEnergy'], axis = 1, inplace=True)
    return(ICE_data)


'''
Program Begins Here
'''
if __name__ == '__main__':
    
    Gas = pd.read_csv('DataSource/Energy_gas.txt' , delimiter=',')
    Bilayer = pd.read_csv('DataSource/BilayerEnergy.txt', delimiter=',')
    IC_1 = pd.read_csv('DataSource/IC_energy1.txt', delimiter =',')
    IC = IC_calculator(IC_1, Gas, Bilayer)
    print(IC)   
