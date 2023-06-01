#!/usr/bin/env python
# coding: utf-8

# The data saved as the file 'dissertation data.csv' is imported and read using the pandas module 
import pandas as pd
import numpy as np

data = pd.read_csv('dissertation data.csv')
print(data)

# The HOMO and LUMO columns are stored in seperate variables 
x = data['HOMO']
print(x)
x2 = data['LUMO']
print(x2)

# The HOMO and LUMO data is plotted as a histogram using the seaborn and matplotlab modules
# They are first plotted with the frequency on the y-axis 
import seaborn as sns
import matplotlib.pyplot as plt

fig, ax=plt.subplots()
for a in [x]:
    sns.histplot(a, bins = range(-12, 2, 1), ax = ax, kde = True, legend = True, color = "#12A6D5").set(title='Distribution of HOMO and LUMO Energies in Dataset', xlabel = "Energy (eV)", ylabel = "Frequency")

for i in [x2]:
    sns.histplot(i, bins = range(-12, 2, 1), ax = ax, kde = True, legend = True, color = "#1256D5")
    
ax.legend(loc = 'upper right', labels=['HOMO', 'LUMO'])


# The histogram is normalised using by setting the stat to density so it can be compared to the carbazole derivatives
x = data['HOMO']
x2 = data['LUMO']

fig, ax = plt.subplots()
sns.histplot(x, bins=range(-12, 3), ax=ax, kde=True, legend=True, color="#12A6D5", stat="density")
sns.histplot(x2, bins=range(-12, 3), ax=ax, kde=True, legend=True, color="#1256D5", stat="density")

ax.set_xlim(-12, 2)
ax.set_ylim(0, 2.00)

ax.set(title='Distribution of HOMO and LUMO Energies in Dataset', xlabel="Energy (eV)", ylabel="Density")
ax.legend(loc='upper right', labels=['HOMO', 'LUMO'])

plt.show()


# The columns containing the energy levels of the three singlet states are imported and stored as variables.
s1 = data['E(S1)']
s2 = data['E(S2)']
s3 = data ['E(S3)']

print(s1)
print(s2)
print(s3)

# The variables containing the singlet energy data are plotted on a histogram and normalised 
fig, ax=plt.subplots()
for singlet1 in [s1]:
    sns.histplot(singlet1, bins = range(-1, 6, 1), ax = ax, kde = True, legend = True, color = "green", stat="density").set(title='Distribution of Singlet Energy Levels in Dataset', xlabel = "Energy (eV)", ylabel = "Denisty")

for singlet2 in [s2]:
    sns.histplot(singlet2, bins = range(-1, 6, 1), ax = ax, kde = True, legend = True, color = "red", stat="density")
    
for singlet3 in [s1]:
    sns.histplot(singlet3, bins = range(-1, 6, 1), ax = ax, kde = True, legend = True, color = "orange", stat="density")
   
ax.legend(loc = 'upper right', labels=['E(S1)', 'E(S2)', 'E(S3)'], title = 'Singlet Energy Level')



# The columns containing the oscillating strengths for each of the singlet states are imported and stored as variables.
f1 = data['f(S1)']
f2 = data['f(S2)']
f3 = data ['f(S3)']

print(f1)
print(f2)
print(f3)

# A historgram is plotted for the oscillator strength and normalised
fig, ax=plt.subplots()
for field1 in [f1]:
    sns.histplot(field1, bins = range(0, 6, 1), ax = ax, kde = True, legend = True, color = "blue", stat="density").set(title='Distribution of Oscillator Strength in Dataset', xlabel = "Oscillator Strength", ylabel = "Density")

for field2 in [f2]:
    sns.histplot(field2, bins = range(0, 6, 1), ax = ax, kde = True, legend = True, color = "pink", stat="density")
    
for field3 in [f1]:
    sns.histplot(field3, bins = range(0, 6, 1), ax = ax, kde = True, legend = True, color = "purple", stat="density")
    
ax.legend(loc = 'upper right', labels=['f(S1)', 'f(S2)', 'f(S3)'], title = 'Singlet Energy Level')

# The columns containing the energy levels of the three triplet states are imported and stored as variables.
t1 = data['E(T1)']
t2 = data['E(T2)']
t3 = data ['E(T3)']
print(t1)
print(t2)
print(t3)

# The variables containing the singlet energy data are plotted on a histogram and normalised.
fig, ax=plt.subplots()
for triplet1 in [t1]:
    sns.histplot(triplet1, bins = range(0, 6, 1), ax = ax, kde = True, legend = True, color = "#003CF6", stat = "density")

for triplet2 in [t2]:
    sns.histplot(triplet2, bins = range(0, 6, 1), ax = ax, kde = True, legend = True, color = "#CE1FF1", stat = "density")
    
for triplet3 in [t3]:
    sns.histplot(triplet3, bins = range(0, 6, 1), ax = ax, kde = True, legend = True, color = "#00FFAA", stat = "density")

ax.set_xlim(-1, 6)
ax.set_ylim(0, 3.4)

ax.set(title='Distribution of Triplet Energy Levels in Dataset', xlabel="Energy (eV)", ylabel="Density")
ax.legend(loc='upper right', labels=['HOMO', 'LUMO'])

plt.show()    
ax.legend(loc = 'upper right', labels=['E(T1)', 'E(T2)', 'E(T3)'], title = 'Triplet Energy Level')

# The column containing the value for the S1-T1 gap is imported and stored as a variable
energy_gap = data['S1-T1 gap']
print(energy_gap)

# The values for the S1-T1 gap are plotted on a histogram and normalised
fig, ax=plt.subplots()
sns.histplot(energy_gap, bins = range(0, 4, 1), ax = ax, kde = True, color = "#75FF94", stat = "density").set(title='Distribution of \u0394E(S1-T1) in Database', xlabel = "Energy (eV)", ylabel = "Density")
    
ax.set_ylim(0, 1.0)
ax.set_xlim(-0.15, 3.00)


# The column containing the value for the HOMO-LUMO gap is imported and stored as a variable
HOMO_LUMO = data['HOMO-LUMO Gap']
print(HOMO_LUMO)

# The values for the HOMO-LUMO gap are plotted on a histogram and normalised
fig, ax=plt.subplots()

sns.histplot(HOMO_LUMO, bins = range(0, 9, 1), ax = ax, kde = True, color = "orange", stat="density").set(title='Distribution of HOMO-LUMO gap in Dataset', xlabel = "Energy (eV)", ylabel = "Denisty")

# The values for the S1-T2 gap are plotted on a histogram and normalised
T2_gap = data['S1-T2 gap']
fig, ax=plt.subplots()
sns.histplot(T2_gap, bins = range(0, 4, 1), ax = ax, kde = True, color = "#9A6FF5", stat = "density").set(title='Distribution of \u0394E(S1-T2) in Database', xlabel = "Energy (eV)", ylabel = "Density")   
ax.set_ylim(0, 2.0)
ax.set_xlim(-0.15, 3.00)

# A list is created using the 'tolist()' function that converts the data imported from the T1 column from an array to a list
t1_list = data['E(T1)'].tolist()
print(t1_list)

# The list containing the T1 energies is sorted from largest to smallest 
t1_list.sort(reverse=True)
print(t1_list)

# The list containing the T1 energies is sorted from smallest to largest
t1_list.sort()
print(t1_list)


# RDKit will now be used to visualise the molecules with the largest and smallest T1 energies starting with the largest T1 energy
# The SMILES of the molecule is stored in an array called subArr
# The SMILES of the molecule is converted to a molecule object using the RDKit function MolFromSmiles
# The molecule object is then converted to a 2D image using the RDKit function Draw.MolToImage
# The image is then displayed using the IPython function display 
import rdkit
from rdkit import Chem

subArr = ['C[N+](C)(C)CCOS(=O)(=O)[O-]']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver

# Second largest T1 energy
subArr = ['OP(O)(=O)C[NH+](CP(O)(=O)[O-])C1CCCCC1']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Third largest T1 energy
subArr = ['CCCCOP(O)(=O)C[NH+](CC[NH+](CP(O)(=O)OCCCC)CP(=O)([O-])OCCCC)CP(=O)([O-])OCCCC']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Fourth largest T1 energy
subArr = ['CP(C)(=O)CCC[NH+](CP(O)(O)=O)CP(O)(=O)[O-]']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Fifth largest T1 energy
subArr = ['OC(CC[NH+]1CCCC1)(P(O)(O)=O)P(O)(=O)[O-]']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Sixth largest T1 energy 
from rdkit import Chem

subArr = ['CC(CNC(=O)NC1CC2CCC1(CC2)NC(=O)NCCC(=O)OCc1ccccc1)NC(=O)NC12CCC(CC1)CC2NC(=O)NCC(Cc1ccccc1)NC(=O)NC12CCC(CC1)CC2NC(=O)NCC(C)NC(=O)NC12CCC(CC1)CC2NC(=O)NCC(Cc1ccccc1)NC(=O)Nc1ccc(Br)cc1']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Seventh largest T1 energy
subArr = ['OP(O)(=O)C[NH+](CCCCCC[NH+](CP(O)(O)=O)CP(O)(=O)[O-])CP(O)(=O)[O-]']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Eighth largest T1 energy
subArr = ['CCOP(=O)(OCC)C(Nc1ccc(OC)cc1)C1=CN(N=C1c1ccc(OC)cc1)c1ccccc1']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Ninth largest T1 energy
subArr = ['O=C([O-])C1CCC[NH+]1Cc1ccccc1']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Tenth largest T1 energy
subArr = ['OCCN1CC[NH+](CC1)CCS(=O)(=O)[O-]']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Smallest T1 energy
subArr = ['CC(C)(C)N1[Si](C)(C)N(C(C)(C)C)[Si]21N=N[Si]1(N=N2)N(C(C)(C)C)[Si](C)(C)N1C(C)(C)C']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Second smallest T1 energy
subArr = ['BrC1=C(Br)C(=C(C#N)C#N)[C-]([n+]2ccccc2)C1=C(C#N)C#N']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Third smallest T1 energy
subArr = ['CC1=C(C(c2ccccc2)C(=C(N)O1)C#N)N1N=NC(=N1)N(=O)=O']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))
    
print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Fourth smallest T1 energy
subArr = ['CC(CCC(=C)C(C)C(O)=O)C1CCC2C3=C(C(=O)CC12C)C1(C)CCC(=O)C(C)C1CC3']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Fifth smallest T1 energy
subArr = ['CC1([N+](=C2C=CC(=CC2=[N+]1[O-])N(=O)=O)[O-])c1ccc(cc1)C(O)=O']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Sixth smallest T1 energy
subArr = ['O=N(=O)C1=CC2=[N+]([O-])C3(CCCC3)[N+](=C2C=C1)[O-]']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Seveneth smallest T1 energy
subArr = ['C[Si](C)(C)C1(CC1S(=O)(=O)c1ccccc1)[Si](=O)(=O)c1ccccc1']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Eighth smallest T1 energy
subArr = ['CC(C)(C)C1=CC(=CC=C2C=C(C(=O)C(=C2)C(C)(C)C)C(C)(C)C)C=C(C1=O)C(C)(C)C']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Ninth smallest T1 energy
subArr = ['O=N(=O)c1ccc(N=NC(c2ccccc2)=C(N=Nc2ccc(cc2N(=O)=O)N(=O)=O)c2ccccc2)c(c1)N(=O)=O']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver


# Tenth smallest T1 energy
subArr = ['Cc1ccc2SC3=C(Sc2c1)C(=NC#N)c1ccccc1C3=NC#N']
lists_subgroup = []
for j in subArr:
    lists_subgroup.append(Chem.MolToSmiles(Chem.MolFromSmiles(j)))

print('smiiiiy: ',lists_subgroup[0])
ver = Chem.MolFromSmiles(lists_subgroup[0])
ver



# The S1-T1 gap is stored in a new variable 
# A for loop is used to iterate over the list of S1-T1 gaps to check for values less than 0.2
# The values less than 0.2 are stored in a new list from smallest to largest
# The new list is printed
gap = data['S1-T1 gap']

small_gap = [item for item in gap if item < 0.2]
small_gap.sort()

print(small_gap)

# This programme uses linear regression to predict the LUMO energy of a molecule
# The modules used are numpy and sklearn
import numpy as np
from sklearn.linear_model import LinearRegression

# Each value in the HOMO and LUMO columns is casted as a float
HOMO_floats = [float(i) for i in x]

LUMO_floats = [float(i) for i in x2]


# The HOMO and LUMO lists are then converted to numpy arrays and reshaped to fit the model
HOMO = np.array(HOMO_floats).reshape(-1, 1)
LUMO = np.array(LUMO_floats).reshape(-1, 1)

# A linear regression model is created and fitted to the data
regressor = LinearRegression()
regressor.fit(HOMO, LUMO)

# A new HOMO energy is entered and the LUMO is predicted
new_HOMO = [[3.5]] 
predicted_LUMO = regressor.predict(new_HOMO)

# The predicted LUMO is printed
print("Predicted LUMO:", predicted_LUMO[0][0])


# The same logic as above is followed to predict the energy of S1 from the HOMO-LUMO gap
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

data = pd.read_csv('dissertation data.csv')
gap = data['HOMO-LUMO Gap']
S1 = data['E(S1)']

HOMO_LUMO_floats = [float (i) for i in gap]
S1_floats = [float (i) for i in S1]

HOMO_LUMO_gap_array = np.array(HOMO_LUMO_floats).reshape(-1, 1)
S1_array = np.array(S1_floats).reshape(-1, 1)

regressor = LinearRegression()
regressor.fit(HOMO_LUMO_gap_array, S1_array)


new_gap = [[3.2]]
predicted_S1 = regressor.predict(new_gap)

print("Predicted S1 energy: ", predicted_S1[0][0])


# The data is split into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(HOMO_LUMO_gap_array, S1_array, test_size=0.2, random_state=42)

# A linear regression model is created and fitted to the training data
regressor = LinearRegression()
regressor.fit(X_train, y_train)

# The performance of the model is evaluated on the testing data
score = regressor.score(X_test, y_test)
print("Model score on testing data: {:.2f}".format(score))

# A score of 0.35 was obtained, suggesting that prediciting the S1 from the HOMO-LUMO gap is not an accurate procedure to follow

