#!/usr/bin/env python
# coding: utf-8

# This programme analyses the database containing the carbazole derivatives
# The neccessary modules an csv. data are imported 
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

data = pd.read_csv('carbazole derivatives from using a dictionary.csv')
print(data)


# The HOMO and LUMO data is plotted as a histogram using the seaborn and matplotlab modules
x = data['HOMO']
x2 = data['LUMO']

import seaborn as sns
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

sns.histplot(x, bins=range(-12, 3), ax=ax, kde=True, legend=True, color="#12A6D5")
sns.histplot(x2, bins=range(-12, 3), ax=ax, kde=True, legend=True, color="#1256D5")

ax.set_xlim(-12, 2)

ax.set(title='Distribution of HOMO and LUMO Energies in Dataset', xlabel="Energy (eV)", ylabel="Frequency")
ax.legend(loc='upper right', labels=['HOMO', 'LUMO'])


# The histogram is normalised by setting the stat to density so that it can be compared to the histogram for the 
# Whole databse
import seaborn as sns
import matplotlib.pyplot as plt

x = data['HOMO']
x2 = data['LUMO']

fig, ax = plt.subplots()

sns.histplot(x, bins=range(-12, 3), ax=ax, kde=True, legend=True, color="#12A6D5", stat="density")
sns.histplot(x2, bins=range(-12, 3), ax=ax, kde=True, legend=True, color="#1256D5", stat="density")

ax.set_xlim(-12, 2)

ax.set(title='Distribution of HOMO and LUMO Energies in Carbazole Derivatives', xlabel="Energy (eV)", ylabel="Density")
ax.legend(loc='upper right', labels=['HOMO', 'LUMO'])

plt.show()


# The columns containing the energy levels of the three singlet states are imported and stored as variables.
# The variables containing the singlet energy data are plotted on a histogram and normalised 
s1 = data['E(S1)']
s2 = data['E(S2)']
s3 = data ['E(S3)']

print(s1)
print(s2)
print(s3)

fig, ax=plt.subplots()
for singlet1 in [s1]:
    sns.histplot(singlet1, bins = range(-1, 6, 1), ax = ax, kde = True, legend = True, color = "green", stat = "density").set(title='Distribution of Singlet Energy Levels in Carbazole Derivatives', xlabel = "Energy (eV)", ylabel = "Density")

for singlet2 in [s2]:
    sns.histplot(singlet2, bins = range(-1, 6, 1), ax = ax, kde = True, legend = True, color = "red", stat = "density")
    
for singlet3 in [s1]:
    sns.histplot(singlet3, bins = range(-1, 6, 1), ax = ax, kde = True, legend = True, color = "orange", stat = "density")

ax.legend(loc = 'upper right', labels=['E(S1)', 'E(S2)', 'E(S3)'], title = 'Singlet Energy Level')


# The columns containing the oscillating strengths for each of the singlet states are imported and stored as variables.
# A historgram is plotted for the oscillator strength and normalised 
f1 = data['f(S1)']
f2 = data['f(S2)']
f3 = data ['f(S3)']

print(f1)
print(f2)
print(f3)

fig, ax=plt.subplots()
for field1 in [f1]:
    sns.histplot(field1, bins = range(0, 6, 1), ax = ax, kde = True, legend = True, color = "blue", stat="density").set(title='Distribution of Oscillator Strength in Carbazole Derivatives', xlabel = "Oscillator Strength", ylabel = "Density")

for field2 in [f2]:
    sns.histplot(field2, bins = range(0, 6, 1), ax = ax, kde = True, legend = True, color = "pink", stat="density")
    
for field3 in [f1]:
    sns.histplot(field3, bins = range(0, 6, 1), ax = ax, kde = True, legend = True, color = "purple", stat="density")

ax.legend(loc = 'upper right', labels=['f(S1)', 'f(S2)', 'f(S3)'], title = 'Singlet Energy Level')


# The columns containing the energy levels of the three triplet states are imported and stored as variables.
# The three triplet energies T1, T2 and T3 are then plotted on the same histogram which is normalised
t1 = data['E(T1)']
t2 = data['E(T2)']
t3 = data ['E(T3)']

import seaborn as sns
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

ax.set_xlim(-2, 6)
ax.set_ylim(0, 3.4)


sns.histplot(t1, bins=range(0, 6), ax=ax, kde=True, legend=True, color="#003CF6", stat = "density")
sns.histplot(t2, bins=range(0, 6), ax=ax, kde=True, legend=True, color="#CE1FF1", stat = "density")
sns.histplot(t3, bins=range(0, 6), ax=ax, kde=True, legend=True, color="#00FFAA", stat = "density")


ax.set(title='Distribution of Triplet Energy Levels in Carbazole Derivatives', xlabel="Energy (eV)", ylabel="Density")
ax.legend(loc='upper right', labels=['E(T1)', 'E(T2)', 'E(T3)'], title='Triplet Energy Level')


# S1-T2 gap

T2_gap = data['S1-T2 gap']
fig, ax=plt.subplots()
sns.histplot(T2_gap, bins = range(0, 4, 1), ax = ax, kde = True, color = "#9A6FF5", stat = "density").set(title='Distribution of \u0394E(S1-T2) in Database', xlabel = "Energy (eV)", ylabel = "Density")
    
ax.set_ylim(0, 2.0)
ax.set_xlim(-0.15, 3.00)

# The values for the S1-T1 gap are plotted on a histogram and normalised
energy_gap = data['S1-T1 gap']
print(energy_gap)

fig, ax=plt.subplots()

sns.histplot(energy_gap, bins = range(0, 4, 1), ax = ax, kde = True, color = "#75FF94", stat = "density").set(title='Distribution of \u0394E(S1-T1) in Carbazole Derivatives', xlabel = "Energy (eV)", ylabel = "Density")
    
    
# The values for the HOMO-LUMO gap are plotted on a histogram and normalised
HOMO_LUMO = data['HOMO-LUMO Gap']
print(HOMO_LUMO)

fig, ax=plt.subplots()
sns.histplot(HOMO_LUMO, bins = range(0, 9, 1), ax = ax, kde = True, color = "orange", stat="density").set(title='Distribution of HOMO-LUMO gap in Carbazole Derivatives', xlabel = "Energy (eV)", ylabel = "Density")
    

# The singlet fission rule is now tested on the carbazole derivatives
# The carbazole S1 and T1 values are stored in variables and the deviation is set to 1
# An empty list is created that will store all of the molecules that satisfy the equatiion within this deviation
carbazole_s1 = 4.02
carbazole_t1 = 3.21
deviation = 1
smiles_list = []

# Two lists are created from the T1 column and SMILES column in the database
t1_list = data['E(T1)']
smiles_list_all = data['SMILES']

# the T1 list is looped over to check if it statisfies the equation with a small enough deviation
for i, t1 in enumerate(t1_list):
    sum_t1 = carbazole_t1 + t1

# If the difference between the sum of the triplets and the carbazole S1 is smaller than the deviation, the molecule
# SMILES is appended to smiles_list
    if abs(sum_t1 - carbazole_s1) <= deviation:
        smiles_list.append(smiles_list_all[i])
    

# The smiles_list is printed
print(smiles_list)


# As these SMILES were edited to create the .xyz files and to visualise the molecule the full smile is needed, 
# The edited SMILES are added to the Cz derivatives database so they can be matched to their full SMILES after 
# the geometry optimisation calculation


# The neccessary modules an csv. data are imported 
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

data = pd.read_csv('carbazole derivatives from using a dictionary.csv')


# This programme is designed to take a list of SMILES and convert them into .xyz files
# The SMILES are taken from a .csv file and an empty list is created to store the SMILES
SMILES_data = data['SMILES']
SMILES_data_list = []
res_list = []


# The SMILES are then appended to the empty list
for SMILES in SMILES_data: 
    SMILES_data_list.append(SMILES)

# The list is then printed to check that it has worked   
print(SMILES_data_list)
pass_list = []

# The modules are imported
import numpy as np
#import pytest
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdmolops
from urllib.request import urlopen
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import sys
import tkinter as tk
  
    
# This function takes the SMILES and converts them into .xyz files
# It also prints the SMILES that it has converted
# It returns the .xyz file and the name of the file is the SMILES with the characters removed
# The smiles that didn't work are printed
def generate_structure_from_smiles(smiles): 
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        # status = AllChem.UFFOptimizeMolecule(mol)  #Optimises Structure
        # conformer = mol.GetConformer()
        # coordinates = conformer.GetPositions()
        # coordinates = np.array(coordinates)
        # atoms = get_atoms(mol)
        print(smiles)
        res = ids
        res = res.replace('#','')
        res = res.replace('(','')
        res = res.replace(')','')
        res = res.replace('[','')
        res = res.replace(']','')
        res = res.replace('=','')
        res = res.replace('+','')
        res = res.replace('-','')
        res_list.append(res) # add res to res_list
        return Chem.MolToXYZFile(mol,res + ".xyz")
        #Returns named .xyz file
    except: 
        print(f"SMILES that didn't work: {smiles}")

# for each SMILES in the list, the function is called
# The SMILES that have been converted are appended to the pass list    
for ids in SMILES_data:
    generate_structure_from_smiles(ids)
    pass_list.append(ids)
    print(len(res_list))
    
data['res'] = res_list
data.to_csv('carbazole derivatives from using a dictionary.csv', index=False)


# After running the above programme, 214 .xyz files were genereated 
# However, there were 285 carbazole derivatives in the database 
# This programme checks for duplicateds in the database to account for this difference
SMILES_data = data['SMILES']
df = pd.DataFrame({'SMILES': SMILES_data})
duplicates = 2 

# Count the number of times each value appears in the 'SMILES' column.
# If a value appears twice, it will have a count of 2.
duplicate_counts = df.value_counts()
exact_duplicates = duplicate_counts[duplicate_counts == 2]

# print the results
if len(exact_duplicates) > 0:
    print("There are", len(exact_duplicates), "values with exactly two duplicates in the 'SMILES' column:")
    print(exact_duplicates)
else:
    print(f"There are no values with exactly {duplicates} duplicates in the 'SMILES' column.")


# There are also two molecules that appear five times in the database, as shown here
SMILES_data = data['SMILES']
df = pd.DataFrame({'SMILES': SMILES_data})
duplicates = 5

# Count the number of times each value appears in the 'SMILES' column.
# If a value appears five times, it will have a count of 5
duplicate_counts = df.value_counts()
exact_duplicates = duplicate_counts[duplicate_counts == 5]


if len(exact_duplicates) > 0:
    print("There are", len(exact_duplicates), "values with exactly five duplicates in the 'SMILES' column:")
    print(exact_duplicates)
else:
    print(f"There are no values with exactly {duplicates} duplicates in the 'SMILES' column.")


# The moleucles that passed the code creating the .xyz files are printed 
print('Items that pass:', pass_list)

print(len(pass_list))


# This programme creates a folder for each file in a directory and moves the file into the new folder
# A folder called carbazole_structures was created on the local computer and the .xyz files were moved there
# This cell was then run
# This is necessary for the gaussian calculation so that each molecule can be accessed via its own directory
# The modules are imported and the folder is defined
# The glob module is used to find all files in the folder
# The os module is used to create a new directory for each file
# The shutil module is used to move the file into the new directory
import glob, os, shutil

folder = '/Users/artemiswebster/Downloads/carbazole_structures'
for file_path in glob.glob(os.path.join(folder, '*.*')):
    new_dir = file_path.rsplit('.', 1)[0]
    os.mkdir(os.path.join(folder, new_dir))
    shutil.move(file_path, os.path.join(new_dir, os.path.basename(file_path))

