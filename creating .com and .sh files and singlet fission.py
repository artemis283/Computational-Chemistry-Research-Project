#!/usr/bin/env python
# coding: utf-8

# In[32]:


# This programme uses the singlet fission rule to check for guest-host matches in the database
# Pandas, numpy, seaborn and matplotlib are imported and the database is read using pandas
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

data = pd.read_csv('dissertation data.csv')
print(data)

# The carbazole S1 and T1 values are stored in variables and the tolerance is set to 0.0002
# An empty list is created that will store all of the molecules that satisfy the equatiion within this tolerance
carbazole_s1 = 4.02
carbazole_t1 = 3.21
tolerance = 0.005
smiles_list = []

# Two lists are created from the T1 column and SMILES column in the database
t1_list = data['E(T1)']
smiles_list_all = data['SMILES']

# the T1 list is looped over to check if it statisfies the equation with a small enough tolerance
for i, t1 in enumerate(t1_list):
    sum_t1 = carbazole_t1 + t1

# If the difference between the sum of the triplets and the carbazole S1 is smaller than the tolerance, the molecule
# SMILES is appended to smiles_list
    if abs(sum_t1 - carbazole_s1) <= tolerance:
        smiles_list.append(smiles_list_all[i])


# In[33]:


# smiles_list is printed
print(smiles_list)


# In[3]:


# rdkit is imported and this programme takes a molecule's SMILES as an input to visualise the molecule
# this molecule satisfies the singlet fission rule with the smallest tolerance of 0.0002
from rdkit import Chem

subArr = ['CC(C)[Si](C#CC1=C2Oc3cccc(C#N)c3N=C2C(=C2Oc3cccc(C#N)c3N=C12)C#C[Si](C(C)C)(C(C)C)C(C)C)(C(C)C)C(C)C']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[9]:


# FC(F)(F)c1c2cc3ccccc3cc2c(c2cc3ccccc3cc12)C(F)(F)F

from rdkit import Chem

subArr = ['FC(F)(F)c1c2cc3ccccc3cc2c(c2cc3ccccc3cc12)C(F)(F)F']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[10]:


# CCSC1=C2C(=S=C(SCC)C2=C(S1)SC(C)(C)C)SC(C)(C)C

from rdkit import Chem

subArr = ['CCSC1=C2C(=S=C(SCC)C2=C(S1)SC(C)(C)C)SC(C)(C)C']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[11]:


from rdkit import Chem

subArr = ['CC(C)[Si](C#CC1=C2Oc3cccc(C#N)c3N=C2C(=C2Oc3cccc(C#N)c3N=C12)C#C[Si](C(C)C)(C(C)C)C(C)C)(C(C)C)C(C)C']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[12]:


from rdkit import Chem

subArr = ['c1ccc(cc1)c1c2ccccc2c(c2ccccc2)c2nc3c(nc12)c(c1ccccc1)c1ccccc1c3c1ccccc1']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[13]:


from rdkit import Chem

subArr = ['CCCCN1C(Sc2ccccc12)=CC1=C([O-])C(=CC2=[N+](CCCC)c3ccccc3S2)C1=C(C#N)C#N']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[14]:


from rdkit import Chem

subArr = ['CN(C)c1ccc(cc1)N=[NH+]c1ccc(cc1)S(=O)(=O)[O-]']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[15]:


from rdkit import Chem

subArr = ['CCCCCCCCN1C(=[Se]c2ccccc12)C=C1C(=O)C(=CC2=[Se]c3ccccc3N2CCCCCCCC)C1=C(C#N)C#N']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[16]:


from rdkit import Chem

subArr = ['COc1cc2CCN3N(NC(=C3c2cc1OC)C#N)c1ccccc1']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[17]:


from rdkit import Chem

subArr = ['CN1C(=C(c2ccccc12)C1=C(C2=C(C)SC(=C2)c2cccc(c2)C(F)(F)F)C(F)(F)C(F)(F)C1(F)F)C']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[8]:


# This script generates .com files and submission scripts for Gaussian calculations from .xyz files
# It adds the .com files and submission scripts to the same directory as the .xyz files
import glob, os, shutil

xyz_dir = '/Users/artemiswebster/Downloads/carbazole_structures'

# set the path to the directory where you want to save the .com files and submission scripts
gau_dir = '/data/scratch/ca20056/carbazole_calculations/carbazole_structures/'

# loop over all the .xyz files in the directory
for file in os.listdir(xyz_dir):
    if not file.endswith('.sh') and not file.endswith('.DS_Store'):
        file_path = f"/Users/artemiswebster/Downloads/carbazole_structures/{file}/{file}.xyz"
        print(f"file path: {file_path}")
        # generate the name of the .com file and submission script based on the name of the .xyz file
        basename = os.path.splitext(file)[0]
        basename = os.path.splitext(file)[0]
        file_dir = os.path.join(xyz_dir, file)
        if os.path.isdir(file_dir):
            opt_com_path = os.path.join(gau_dir, f'{basename}/' + 'Opt_T1.com')
            # opt_com_file = os.path.join(gau_dir, f'{basename}/' + 'Opt_T1.com') # change this
            single_com_path = os.path.join(gau_dir, f'{basename}/' + 'SP_T1.com') # change this
            sub_script = os.path.join(gau_dir, basename + '.sh')
            print(f"sub_script: {sub_script}")
            # rest of your code
        else:
            print(f"{file} is not a directory, skipping...")

        opt_com_file = os.path.join(xyz_dir, f'{file}/' + 'Opt_T1.com')
        single_com_file = os.path.join(xyz_dir, f'{file}/' + 'SP_T1.com')
        sub_script = os.path.join(xyz_dir, basename + '.sh')
        print(f"sub_script: {sub_script}")

        # read the contents of the .xyz file and generate the contents of the .com file
        print(f"Opening file {file_path}")
        # read the contents of the .xyz file and generate the contents of the .com file
        with open(file_path, 'r') as f:
            #Add XYZ code to end of this file
            xyzBlock = f.readlines()[2:]
            xyz_contents = '{}\n'.format(''.join(xyzBlock))
            print(f"xyz contents: {xyz_contents}") 
            opt_com_contents = f"""%mem=64GB 
%nproc=8
%chk={basename}.chk
# tda=(nstates=3,triplets,root=1) opt=maxcycle=100 M062X def2SVP symmetry=none density=(current)

{basename}

0 1
{xyz_contents}
"""
            
            SP_com_contents = f"""%mem=64GB
%nproc=8
%chk={basename}.chk
# tda=(nstates=5,50-50,root=1) M062X def2SVP symmetry=none density=(current) Geom=Check Guess=Read

mh

0 1
""" 
            
            with open(opt_com_file, 'w') as f:
                f.write(opt_com_contents)
            with open(single_com_file, 'w') as f:
                f.write(SP_com_contents)


        # generate the contents of the submission script
        sub_contents = f"""#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=72:0:0
#$ -l h_vmem=8G

dos2unix {opt_com_path}
dos2unix {single_com_path}
module load gaussian/g16
g16 {opt_com_path}
g16 {single_com_path}

"""

    # write the contents of the submission script to disk
    with open(sub_script, 'w') as f:
        f.write(sub_contents)

    print(f"sub contents: {sub_contents}")


# In[9]:


# This programme moves the .sh files into the correct folder
# It also removes the .sh file if it already exists in the folder
# This is so that each calculation can be run on guassian in its own directory for each molecule

import glob, os, shutil

for filename in os.listdir('/Users/artemiswebster/Downloads/carbazole_structures'):
    print(f"Filename: {filename}")
    if filename.endswith('.sh') and filename != '.DS_Store.sh':
        SMILES = filename.strip('.sh')
        print(f"SMILES: {SMILES}")
        shutil.move(f'/Users/artemiswebster/Downloads/carbazole_structures/{filename}', f'/Users/artemiswebster/Downloads/carbazole_structures/{SMILES}/{SMILES}.sh')
    else: 
        foldername = filename
        print(f"Foldername: {foldername}")


# In[ ]:




