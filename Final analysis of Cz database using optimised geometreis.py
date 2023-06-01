#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Importing necessary modules and S1 data, printing data and checking the length

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

S1_data = pd.read_csv('S1_output copy.csv')
print(S1_data)

S1_values = S1_data['Optimized S1 (eV)']
S1_SMILES = S1_data['SMILE']

len(S1_values)


# In[2]:


# Importing T1 data, printing data and checking the length
T1_data = pd.read_csv('T1_output copy.csv')
print(T1_data)

T1_values = T1_data['Optimized T1 (eV)']
T1_SMILES = T1_data['SMILE']

len(T1_values)


# In[3]:


# Plotting a line graph for all of the derivatives in the S1 and T1 databases
from matplotlib.ticker import FormatStrFormatter

fig, ax = plt.subplots(1, 1, figsize=(8, 4))
sns.scatterplot(x=np.arange(len(T1_values)), y=T1_values , ax=ax,marker="_",alpha=1,s=700,color="pink", label='T1')
sns.scatterplot(x=np.arange(len(S1_values)), y=S1_values, ax=ax,marker="_",alpha=1,s=700,color="purple", label='S1')

plt.ylabel('Energy (eV)', fontsize=16)
plt.xlabel('Molecule number', fontsize = 16)
plt.tick_params(axis='both', which='major', labelsize=16)

ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
sns.set_style('white', {'axes.linewidth': 0.5})
ax.tick_params(bottom=True, left=True)

plt.plot([], linestyle='-')
plt.show()


# In[4]:


# Testing the singlet fission rule on optimised geometries
# The carbazole S1 and T1 values are stored in variables and the tolerance is set to 0.65
# An empty list is created that will store all of the molecules that satisfy the equatiion within this tolerance
carbazole_s1 = 4.02
carbazole_t1 = 3.21
tolerance = 0.65
smiles_list = []

# Two lists are created from the T1 column and SMILES column in the database
t1_list = T1_data['Optimized T1 (eV)']
smiles_list_all = T1_data['SMILE']

# the T1 list is looped over to check if it statisfies the equation with a small enough tolerance
for i, t1 in enumerate(t1_list):
    sum_t1 = carbazole_t1 + t1

# If the difference between the sum of the triplets and the carbazole S1 is smaller than the tolerance, the molecule
# SMILES is appended to smiles_list
    if abs(carbazole_s1 - carbazole_t1 - t1) <= tolerance:
        smiles_list.append(smiles_list_all[i])


# In[5]:


# The smiles_list is printed, to show the molecules that satisfy the critera
print(smiles_list)


# In[6]:


# First moleucle that staisfies heterofission mechanism: Oc1cccc2cccCCCNN3c4ccccc4c4ccccc34nc12

from rdkit import Chem

subArr = ['Oc1cccc2ccc(C=C(C#N)N3c4ccccc4c4ccccc34)nc12']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[7]:


# Second moleucle that staisfies heterofission mechanism: CC1(C)OC(=C(C#N)C#N)C(=C1C=Cc1ccc(cc1)N1c2ccccc2c2ccccc12)C#N

from rdkit import Chem

subArr = ['CC1(C)OC(=C(C#N)C#N)C(=C1C=Cc1ccc(cc1)N1c2ccccc2c2ccccc12)C#N']

lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[8]:


# Third moleucle that staisfies heterofission mechanism: N#CC(=Cc1cccc2ccccc12)N1c2ccccc2c2ccccc12
from rdkit import Chem

subArr = ['N#CC(=Cc1cccc2ccccc12)N1c2ccccc2c2ccccc12']

lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[9]:


# Fourth moleucle that staisfies heterofission mechanism: BrC1=CC=C(C=CN2c3ccccc3c3ccccc23)C2=NSN=C12
from rdkit import Chem

subArr = ['BrC1=CC=C(C=CN2c3ccccc3c3ccccc23)C2=NSN=C12']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[10]:


# Fifth moleucle that staisfies heterofission mechanism: COC(=O)C(C#Cc1ccc(cc1)N1c2ccccc2c2ccccc12)=C(C#Cc1ccc(cc1)N1c2ccccc2c2ccccc12)C(=O)OC
from rdkit import Chem

subArr = ['COC(=O)C(C#Cc1ccc(cc1)N1c2ccccc2c2ccccc12)=C(C#Cc1ccc(cc1)N1c2ccccc2c2ccccc12)C(=O)OC']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[11]:


# Sixth moleucle that staisfies heterofission mechanism: COc1ccc(C=Nc2ccc(cc2)N2c3ccccc3c3ccccc23)c(O)c1

subArr = ['COc1ccc(C=Nc2ccc(cc2)N2c3ccccc3c3ccccc23)c(O)c1']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[12]:


# Seventh moleucle that staisfies heterofission mechanism: c1ccc(cc1)C(=Cc1ccc(cc1)N1c2ccccc2c2ccccc12)c1ccncc1

subArr = ['c1ccc(cc1)C(=Cc1ccc(cc1)N1c2ccccc2c2ccccc12)c1ccncc1']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[13]:


# Eighth moleucle that staisfies heterofission mechanism: COc1ccc(cc1)C1=NN=C(S1)N=Cc1ccc(cc1)N1c2ccccc2c2ccccc12

subArr = ['COc1ccc(cc1)C1=NN=C(S1)N=Cc1ccc(cc1)N1c2ccccc2c2ccccc12']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[55]:


# Ninth moleucle that staisfies heterofission mechanism: c1ccc(cc1)C=CC(=Cc1ccccc1)N1c2ccccc2c2ccccc12

subArr = ['c1ccc(cc1)C=CC(=Cc1ccccc1)N1c2ccccc2c2ccccc12']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[56]:


# Tenth moleucle that staisfies heterofission mechanism: N#CC(C#N)=C(N1c2ccccc2c2ccccc12)C(c1ccccc1)=C(C#N)C#N

subArr = ['N#CC(C#N)=C(N1c2ccccc2c2ccccc12)C(c1ccccc1)=C(C#N)C#N']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[57]:


# Eleventh moleucle that staisfies heterofission mechanism: c1ccc(C=Nc2ccc(C=Cc3ccc(cc3)N3c4ccccc4c4ccccc34)cc2)nc1

subArr = ['c1ccc(C=Nc2ccc(C=Cc3ccc(cc3)N3c4ccccc4c4ccccc34)cc2)nc1']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[58]:


# The carbazole S1 and T1 values are stored in variables and the tolerance is set to 0.2, to test a smaller deviation
# An empty list is created that will store all of the molecules that satisfy the equatiion within this tolerance
carbazole_s1 = 4.02
carbazole_t1 = 3.21
tolerance = 0.2
smiles_list = []

# Two lists are created from the T1 column and SMILES column in the database
t1_list = T1_data['Optimized T1 (eV)']
smiles_list_all = T1_data['SMILE']

# the T1 list is looped over to check if it statisfies the equation with a small enough tolerance
for i, t1 in enumerate(t1_list):
    sum_t1 = carbazole_t1 + t1

# If the difference between the sum of the triplets and the carbazole S1 is smaller than the tolerance, the molecule
# SMILES is appended to smiles_list
    if abs(sum_t1 - carbazole_s1) <= tolerance:
        smiles_list.append(smiles_list_all[i])


# In[59]:


print(smiles_list)


# In[60]:


# Molecules that satisty the heterofission mechanism with a 0.2 deviation:
# COCOCCCc1ccccc1N1c2ccccc2c2ccccc12CCCc1ccccc1N1c2ccccc2c2ccccc12COOC

subArr = ['COC(=O)C(C#Cc1ccc(cc1)N1c2ccccc2c2ccccc12)=C(C#Cc1ccc(cc1)N1c2ccccc2c2ccccc12)C(=O)OC']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[61]:


# N#CC(=Cc1cccc2ccccc12)N1c2ccccc2c2ccccc12

subArr = ['N#CC(=Cc1cccc2ccccc12)N1c2ccccc2c2ccccc12']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[62]:


# N#CC(C#N)=C(N1c2ccccc2c2ccccc12)C(c1ccccc1)=C(C#N)C#N
# all three moleucles that satisfy the singlet fission rule +/- 0.2, contain an alkyne.

subArr = ['N#CC(C#N)=C(N1c2ccccc2c2ccccc12)C(c1ccccc1)=C(C#N)C#N']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))


print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# In[14]:


# Plotting a line graph for the molecuels that satisfy the heterofission mechanism with a 0.65 deviation
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter

smiles_list_data = pd.read_csv('smiles_list.csv')
smiles_list_T1_values = [0.5524, 1.3962, 0.9177, 1.4199, 0.7682, 0.4211, 0.4988, 1.4367, 1.3297, 0.9452, 0.1732, 3.21]
smiles_list_S1_values = [2.2273, 2.6815, 2.5334, 2.7164, 2.2991, 0.7156, 2.2229, 2.787, 3.0701, 2.1793, 0.5656, 4.02]
x_labels = ['12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23']

fig, ax = plt.subplots(1, 1, figsize=(8, 4))

sns.scatterplot(x=np.arange(len(smiles_list_T1_values)), y=smiles_list_T1_values , ax=ax,marker="_",alpha=1,s=700,color="pink", label='T1')
sns.scatterplot(x=np.arange(len(smiles_list_S1_values)), y=smiles_list_S1_values, ax=ax,marker="_",alpha=1,s=700,color="purple", label='S1')

plt.ylabel('Energy (eV)', fontsize=16)
plt.xlabel('Molecule number', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)

ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
sns.set_style('white', {'axes.linewidth': 0.5})
ax.tick_params(bottom=True, left=True)

ax.legend(fontsize=14, loc='best')

# set the x-axis ticks and labels
ax.set_xticks(np.arange(len(smiles_list_T1_values)))
ax.set_xticklabels(x_labels)

# set the x-axis limit to include all the data points
ax.set_xlim(-0.5, len(smiles_list_T1_values)-0.5)

plt.plot([], linestyle='-')
plt.show()



# In[64]:


# Host and guest T1 energy scatter graph with the molecuels that satisfy the heterofission mechanism with a 0.65 devation
# Carbazole is shown in red 

import matplotlib.pyplot as plt
import numpy as np

t1_values = [0.5524, 1.3962, 0.9177, 1.4199, 0.7682, 0.4211, 0.4988, 1.4367, 1.3297, 0.9452, 0.1732, 3.21]
x_labels = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21']


x = np.arange(len(t1_values))  # create an array of x-values
colors = ['black' if val != 3.21 else 'red' for val in t1_values]
plt.scatter(x_labels, t1_values, c=colors)
ax.set_xticklabels(x_labels)
plt.xlabel('Molecule Number', fontsize=16)
plt.ylabel('T1 Excited State (eV)', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.show()


# In[65]:


# Host and guest T1 energy scatter graph with all of the Cz derivatives that the calculations have been run on
# Carbazole is shown in red

import matplotlib.pyplot as plt
import numpy as np


x = np.arange(len(T1_values))  # create an array of x-values
colors = ['black' if val != 3.21 else 'red' for val in T1_values]
plt.scatter(x, T1_values, c=colors)
plt.xlabel('Molecule', fontsize=16)
plt.ylabel('T1 Excited State (eV)', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.show()

print(x)
print(T1_values)


# In[149]:


# Finding the difference of the S1 Host - T1 host - T1 guest so it can be used in a table

carbazole_s1 = 4.02
carbazole_t1 = 3.21
difference_list = []

# Data is read and the triplet values are stored in a variable
SMILES_data = pd.read_csv('smiles_list.csv')
t1_list = SMILES_data['Optimized T1 (eV)']
SMILES = SMILES_data['SMILE']


# the T1 list is looped over to calculate the difference
for i, t1 in enumerate(t1_list):
    difference = carbazole_s1 - carbazole_t1 - t1
    round_difference = round(difference, 2)
    difference_list.append(round_difference)
    
print(difference_list)


# In[1]:


# Plotting a line graph for the top 5 Cz derivatives molecuels that satisfy the heterofission mechanism
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter

smiles_list_data = pd.read_csv('smiles_list.csv')
smiles_list_T1_values = [0.5524, 0.9177, 0.7682, 0.4988, 0.9452, 3.21]
smiles_list_S1_values = [2.2273, 2.5334, 2.2991, 2.2229, 2.1793, 4.02]
x_labels = ['1', '2', '3', '4', '5', '6']

fig, ax = plt.subplots(1, 1, figsize=(8, 4))

sns.scatterplot(x=np.arange(len(smiles_list_T1_values)), y=smiles_list_T1_values , ax=ax,marker="_",alpha=1,s=700,color="pink", label='T1')
sns.scatterplot(x=np.arange(len(smiles_list_S1_values)), y=smiles_list_S1_values, ax=ax,marker="_",alpha=1,s=700,color="purple", label='S1')

plt.ylabel('Energy (eV)', fontsize=16)
plt.xlabel('Molecule number', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)

ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
sns.set_style('white', {'axes.linewidth': 0.5})
ax.tick_params(bottom=True, left=True)

ax.legend(fontsize=14, loc='best')

# set the x-axis ticks and labels
ax.set_xticks(np.arange(len(smiles_list_T1_values)))
ax.set_xticklabels(x_labels)

# set the x-axis limit to include all the data points
ax.set_xlim(-0.5, len(smiles_list_T1_values)-0.5)

plt.plot([], linestyle='-')
plt.show()


# In[ ]:




