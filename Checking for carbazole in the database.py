#!/usr/bin/env python
# coding: utf-8

# In[7]:


# This programme checks for carbazole derivatives in the data base of semiconductors 
# If a derivative is found it is added to a new database, if the molecule cannot be kekulized, it is added to a list
# This programme was ran a number of times as it required a lot of processing power so the molecules in the database 
# were tested in batches
# Molecules 0 - 1999
# RDKit is imported along with its submodules and pandas and numpy
from rdkit import Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem import rdMMPA
import pandas as pd
import numpy as np

# A function is defined to kekulize aromatic bonds in a molecule
# This is necessary for the RDKit to be able to fragment the molecule
# The function takes a molecule as an argument and returns a kekulized molecule
# The SetBondType() function is used to change the bond type from aromatic to single
# The SetNumExplicitHs() function is used to change the number of explicit hydrogens from 0 to 1
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

    
# The data is read from the csv file and stored in a pandas dataframe
# The SMILES column is stored in a variable
# Two empty lists are created to store the SMILES of molecules that contain carbazole derivatives and molecules that couldn't be kekulized
# A dataframe is created to store the data of molecules that contain carbazole derivatives 
data = pd.read_csv('dissertation data.csv')
SMILES_data = data['SMILES'].iloc[0:2000]
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
  


# The SMILES column is iterated over
# Three empty lists are created to store the fragments of the molecules
# The molecule is converted to a SMILES string and then to a molecule
# The molecule is fragmented using the FragmentMol() function 
# The fragments are iterated over and stored in a list 
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
        
# The fragments are iterated over again and stored in a list 
# The list of lists is converted to a list 
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

# The list is iterated over and the fragments are converted to SMILES strings
# The SMILES strings are iterated over and the carbazole derivatives are checked for
# If the molecule contains a carbazole derivative, the SMILES string is stored in the carbazole list
# If the molecule couldn't be kekulized, the SMILES string is stored in the cannot kekulize list 
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

# The SMILES strings of molecules that contain carbazole derivatives are stored in a dataframe 
# The dataframe is exported to a csv file
carbazole_data_frame = data_frame[data_frame['SMILES'].isin(carbazole_list)]
carbazole_data_frame.to_csv('carbazole derivatives from using a dictionary.csv', header=True, index=False)


cannot_kekulize_data_frame = data_frame[data_frame['SMILES'].isin(cannot_kekulize_list)]
cannot_kekulize_data_frame.to_csv('cannot kekulize list.csv', header=True, index=False)


# In[9]:


# Molecules 2000 - 4999
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
SMILES_data = data['SMILES'].iloc[2000:5000]
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
    





# In[10]:


# Molecules 5000 - 14999
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
SMILES_data = data['SMILES'].iloc[5000:15000]
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
    




# In[11]:


# Molecules 15000 - 24999
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
SMILES_data = data['SMILES'].iloc[15000:25000]
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
    cannot_kekulize_data_frame.to_csv(f2, header=False, index=False)


# In[16]:


# This programme used rdkit to visualise one of the derivatives found as a test to see if it contains carbazole
# O=C1N(C(=O)c2ccccc12)c1ccc2Nc3ccccc3c2c1

from rdkit import Chem

subArr = ['COC(=O)c1ccccc1N1c2ccccc2c2ccccc12']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[ ]:




