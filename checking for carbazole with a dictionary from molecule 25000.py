#!/usr/bin/env python
# coding: utf-8

# In[1]:


# This programme follows the same logic as the file titled, 'Checking for carbazole in the database'
# Due to the large number of molecules, the programme had to be run seperately and in a different file for the 
# second half of the database.
# Molecules 25000 - 34999
from rdkit import Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem import rdMMPA
import pandas as pd
import numpy as np

def kekulize(mol):
    """
    Kekulize aromatic bonds in the given molecule.
    """
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.AROMATIC:
            bond.SetBondType(Chem.rdchem.BondType.SINGLE)
            bond.GetBeginAtom().SetNumExplicitHs(1)
            bond.GetEndAtom().SetNumExplicitHs(1)
    Chem.Kekulize(mol)

    
data = pd.read_csv('dissertation data.csv')
SMILES_data = data['SMILES'].iloc[25000:35000]
carbazole_list = [] 
cannot_kekulize_list = []

data_frame = pd.DataFrame({'ID': data['ID'],
                           'doi': data ['doi'],
                           'formula': data ['formula'],
                           'NAts': data ['NAts'],
                           'SMILES': data ['SMILES'],
                           'HOMO' : data ['HOMO'],
                           'LUMO' : data ['LUMO'],
                           'E(S1)' : data ['E(S1)'],
                           'f(S1)' : data ['f(S1)'],
                           'E(S2)' : data ['E(S2)'],
                           'f(S2)' : data ['f(S2)'],
                           'E(S3)' : data ['E(S3)'],
                           'f(S3)' : data ['f(S3)'],
                           'E(T1)' : data ['E(T1)'],
                           'E(T2)' : data ['E(T2)'],
                           'E(T3)' : data ['E(T3)'],
                           'S1-T1 gap' : data ['S1-T1 gap'],
                           'HOMO positive' : data ['HOMO positive'], 
                           'LUMO positive' : data ['LUMO positive'], 
                           'HOMO-LUMO Gap' : data ['HOMO-LUMO Gap']})
    
for SMILE in SMILES_data:
    preFragList = []
    fragList = []
    finlFragLstHs = []
    # mol1 = Chem.MolFromSmiles(Chem.MolToSmiles(mol1))
    # mol1 = Chem.MolFromSmiles(SMILE)
    try:
        mol1 = AllChem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(SMILE)))
        kekulize(mol1)
        
        frag = rdMMPA.FragmentMol(mol1)

        for individFrags in frag:
            if individFrags[0] != None:
                frag1 = Chem.MolToSmiles(individFrags[0])
                frag2 = Chem.MolToSmiles(individFrags[1])
            else:
                frag1 = ""
                frag2 = Chem.MolToSmiles(individFrags[1])

            
            res = frag1+frag2
            res = res.replace('([*:1])','')
            res = res.replace('[*:1]','')
            res = res.replace('([*:2])','')
            res = res.replace('[*:2]','')
            
            preFragList.append(res)
      

        for string in preFragList: #separate fragments if there's a '.' between them
            subGroup = string.split('.')
            fragList.append(subGroup)

        

        finFragLst = sum(fragList, []) #convert list of list to just one list
        finFragLst = list(dict.fromkeys(finFragLst)) #remove duplicate fragments


        convertedFragLst = []

        for frag in finFragLst:
            try:
                frag = Chem.MolToSmiles(Chem.MolFromSmiles(frag))
                convertedFragLst.append(frag)
            except:
                print('no')

        carbazoleArr = ['C1=CC=C2C(=C1)C3=CC=CC=C3N2']
        lists_subgroup = []
        for j in carbazoleArr:
            lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


        carbazole_SMILES = lists_subgroup[0]


        if carbazole_SMILES in convertedFragLst: 
            carbazole_list.append(SMILE)
            

    except:
        cannot_kekulize_list.append(SMILE)
        

print("Carbazole derivatvies: ", carbazole_list)
print("Molecules that couldn't be kekulized: ", cannot_kekulize_list)

carbazole_data_frame = data_frame[data_frame['SMILES'].isin(carbazole_list)]

with open('carbazole derivatives from using a dictionary.csv', 'a') as f:
    carbazole_data_frame.to_csv(f, header=False, index=False)
    

cannot_kekulize_data_frame = data_frame[data_frame['SMILES'].isin(cannot_kekulize_list)]

with open('cannot kekulize list.csv', 'a') as f2:
    cannot_kekulize_data_frame.to_csv(f2, header=True, index=False)


# In[2]:


# Molecules 35000 - 44999

def kekulize(mol):
    """
    Kekulize aromatic bonds in the given molecule.
    """
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.AROMATIC:
            bond.SetBondType(Chem.rdchem.BondType.SINGLE)
            bond.GetBeginAtom().SetNumExplicitHs(1)
            bond.GetEndAtom().SetNumExplicitHs(1)
    Chem.Kekulize(mol)

    
data = pd.read_csv('dissertation data.csv')
SMILES_data = data['SMILES'].iloc[35000:45000]
carbazole_list = [] 
cannot_kekulize_list = []

data_frame = pd.DataFrame({'ID': data['ID'],
                           'doi': data ['doi'],
                           'formula': data ['formula'],
                           'NAts': data ['NAts'],
                           'SMILES': data ['SMILES'],
                           'HOMO' : data ['HOMO'],
                           'LUMO' : data ['LUMO'],
                           'E(S1)' : data ['E(S1)'],
                           'f(S1)' : data ['f(S1)'],
                           'E(S2)' : data ['E(S2)'],
                           'f(S2)' : data ['f(S2)'],
                           'E(S3)' : data ['E(S3)'],
                           'f(S3)' : data ['f(S3)'],
                           'E(T1)' : data ['E(T1)'],
                           'E(T2)' : data ['E(T2)'],
                           'E(T3)' : data ['E(T3)'],
                           'S1-T1 gap' : data ['S1-T1 gap'],
                           'HOMO positive' : data ['HOMO positive'], 
                           'LUMO positive' : data ['LUMO positive'], 
                           'HOMO-LUMO Gap' : data ['HOMO-LUMO Gap']})
    
for SMILE in SMILES_data:
    preFragList = []
    fragList = []
    finlFragLstHs = []
    # mol1 = Chem.MolFromSmiles(Chem.MolToSmiles(mol1))
    # mol1 = Chem.MolFromSmiles(SMILE)
    try:
        mol1 = AllChem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(SMILE)))
        kekulize(mol1)
        
        frag = rdMMPA.FragmentMol(mol1)

        for individFrags in frag:
            if individFrags[0] != None:
                frag1 = Chem.MolToSmiles(individFrags[0])
                frag2 = Chem.MolToSmiles(individFrags[1])
            else:
                frag1 = ""
                frag2 = Chem.MolToSmiles(individFrags[1])

            
            res = frag1+frag2
            res = res.replace('([*:1])','')
            res = res.replace('[*:1]','')
            res = res.replace('([*:2])','')
            res = res.replace('[*:2]','')
            
            preFragList.append(res)
      

        for string in preFragList: #separate fragments if there's a '.' between them
            subGroup = string.split('.')
            fragList.append(subGroup)

        

        finFragLst = sum(fragList, []) #convert list of list to just one list
        finFragLst = list(dict.fromkeys(finFragLst)) #remove duplicate fragments


        convertedFragLst = []

        for frag in finFragLst:
            try:
                frag = Chem.MolToSmiles(Chem.MolFromSmiles(frag))
                convertedFragLst.append(frag)
            except:
                print('no')

        carbazoleArr = ['C1=CC=C2C(=C1)C3=CC=CC=C3N2']
        lists_subgroup = []
        for j in carbazoleArr:
            lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


        carbazole_SMILES = lists_subgroup[0]


        if carbazole_SMILES in convertedFragLst: 
            carbazole_list.append(SMILE)
            

    except:
        cannot_kekulize_list.append(SMILE)
        

print("Carbazole derivatvies: ", carbazole_list)
print("Molecules that couldn't be kekulized: ", cannot_kekulize_list)

carbazole_data_frame = data_frame[data_frame['SMILES'].isin(carbazole_list)]

with open('carbazole derivatives from using a dictionary.csv', 'a') as f:
    carbazole_data_frame.to_csv(f, header=False, index=False)
    


# In[3]:


# Molecules 45000 - 49000
from rdkit import Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem import rdMMPA
import pandas as pd
import numpy as np

def kekulize(mol):
    """
    Kekulize aromatic bonds in the given molecule.
    """
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.AROMATIC:
            bond.SetBondType(Chem.rdchem.BondType.SINGLE)
            bond.GetBeginAtom().SetNumExplicitHs(1)
            bond.GetEndAtom().SetNumExplicitHs(1)
    Chem.Kekulize(mol)

    
data = pd.read_csv('dissertation data.csv')
SMILES_data = data['SMILES'].iloc[45000:50000]
carbazole_list = [] 
cannot_kekulize_list = []

data_frame = pd.DataFrame({'ID': data['ID'],
                           'doi': data ['doi'],
                           'formula': data ['formula'],
                           'NAts': data ['NAts'],
                           'SMILES': data ['SMILES'],
                           'HOMO' : data ['HOMO'],
                           'LUMO' : data ['LUMO'],
                           'E(S1)' : data ['E(S1)'],
                           'f(S1)' : data ['f(S1)'],
                           'E(S2)' : data ['E(S2)'],
                           'f(S2)' : data ['f(S2)'],
                           'E(S3)' : data ['E(S3)'],
                           'f(S3)' : data ['f(S3)'],
                           'E(T1)' : data ['E(T1)'],
                           'E(T2)' : data ['E(T2)'],
                           'E(T3)' : data ['E(T3)'],
                           'S1-T1 gap' : data ['S1-T1 gap'],
                           'HOMO positive' : data ['HOMO positive'], 
                           'LUMO positive' : data ['LUMO positive'], 
                           'HOMO-LUMO Gap' : data ['HOMO-LUMO Gap']})
    
for SMILE in SMILES_data:
    preFragList = []
    fragList = []
    finlFragLstHs = []
    # mol1 = Chem.MolFromSmiles(Chem.MolToSmiles(mol1))
    # mol1 = Chem.MolFromSmiles(SMILE)
    try:
        mol1 = AllChem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(SMILE)))
        kekulize(mol1)
        
        frag = rdMMPA.FragmentMol(mol1)

        for individFrags in frag:
            if individFrags[0] != None:
                frag1 = Chem.MolToSmiles(individFrags[0])
                frag2 = Chem.MolToSmiles(individFrags[1])
            else:
                frag1 = ""
                frag2 = Chem.MolToSmiles(individFrags[1])

            
            res = frag1+frag2
            res = res.replace('([*:1])','')
            res = res.replace('[*:1]','')
            res = res.replace('([*:2])','')
            res = res.replace('[*:2]','')
            
            preFragList.append(res)
      

        for string in preFragList: #separate fragments if there's a '.' between them
            subGroup = string.split('.')
            fragList.append(subGroup)

        

        finFragLst = sum(fragList, []) #convert list of list to just one list
        finFragLst = list(dict.fromkeys(finFragLst)) #remove duplicate fragments


        convertedFragLst = []

        for frag in finFragLst:
            try:
                frag = Chem.MolToSmiles(Chem.MolFromSmiles(frag))
                convertedFragLst.append(frag)
            except:
                print('no')


        carbazoleArr = ['C1=CC=C2C(=C1)C3=CC=CC=C3N2']
        lists_subgroup = []
        for j in carbazoleArr:
            lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


        carbazole_SMILES = lists_subgroup[0]


        if carbazole_SMILES in convertedFragLst: 
            carbazole_list.append(SMILE)
            

    except:
        cannot_kekulize_list.append(SMILE)
        

print("Carbazole derivatvies: ", carbazole_list)
print("Molecules that couldn't be kekulized: ", cannot_kekulize_list)

carbazole_data_frame = data_frame[data_frame['SMILES'].isin(carbazole_list)]

with open('carbazole derivatives from using a dictionary.csv', 'a') as f:
    carbazole_data_frame.to_csv(f, header=False, index=False)
    


# In[ ]:




