#!/usr/bin/env python
# coding: utf-8

# In[7]:


# This programme extracts the excited state T1 value from a Gaussian log file from a geormetry optimisation
# It then writes the T1 value to a csv file along with the file path of the log file

# The necessary modules are imported
import os, sys
import csv
import glob
import pandas as pd
import tempfile
import re
import xlrd


# A function is defined that takes the parent path of the directory containing the log files
# as well as the path of the output csv file as arguments
# The function then searches for the T1 value in the log file and writes the file path and T1 value to the csv file
# The file path contains the molecules SMILES so this could be used to identify which molecule the T1 value belongs to
# The function also checks if the csv file already exists and if it does it appends the new data to the end of the file
def findT1Gauss(parentPath, outputPath):
    searchParam = 'Excited State   1:'
    searchParam2 = '      Triplet-?Sym    '
    count = 0
    if os.path.exists(outputPath):
        outputFile = open(outputPath,'w').close()
    append_write = 'w'
       
    for root, dirs, filepaths in os.walk(parentPath):
        print(f"root: {root}")
        print(f"dirs: {dirs}")
        print(f"filepaths: {filepaths}")
        for filepath in filepaths:
            if filepath.startswith('SP_T1') and filepath.endswith(".log"):
                readFile = open(os.path.join(root, filepath), 'r')
               
                count += 1
                if count > 1:
                    if os.path.exists(outputPath):
                        append_write = 'a' # append if already exists
                    else:
                        append_write = 'w' # make a new file if not
                   
                outputFile = open(outputPath,append_write)
               
                for ES_Line in readFile:
                    if searchParam in ES_Line and searchParam2 in ES_Line:
                        ES_Line_Arr = ES_Line.split(' ')
                        T1_Val = ES_Line_Arr[15]
                        writer = csv.writer(outputFile)
                        if append_write == 'w':
                            writer.writerow([root, T1_Val])
                        else:
                            writer.writerow([root, T1_Val])


 # The parent path of the directory containing the log files is defined
 # The path of the output csv file is defined                           
parentPath = '/Users/artemiswebster/Downloads/carbazole_structures'
outputPath = '/Users/artemiswebster/Downloads/carbazole_structures/T1_output.csv'

# The function is called
findT1Gauss(parentPath, outputPath)


# In[13]:


# This programme follows the same logic as the programme above but extracts the S1 value from the log file
# As the singlet value is not always at the same excited sate, the programme searches for the line containing the
# singlet value and then extracts the value from the line. The value is then written to a csv file along with the
# file path. The programme also checks if the file already exists and if it does, appends the data to the file

import os
import csv

def findS1Gauss(parentPath, outputPath):
    searchParam2 = '      Singlet-?Sym    '
    count = 0
    if os.path.exists(outputPath):
        outputFile = open(outputPath,'w').close()
    append_write = 'w'
       
    for root, dirs, filepaths in os.walk(parentPath):
        print(f"root: {root}")
        print(f"dirs: {dirs}")
        print(f"filepaths: {filepaths}")
        for filepath in filepaths:
            if filepath.startswith('SP_T1') and filepath.endswith(".log"):
                readFile = open(os.path.join(root, filepath), 'r')
               
                count += 1
                if count > 1:
                    if os.path.exists(outputPath):
                        append_write = 'a' # append if already exists
                    else:
                        append_write = 'w' # make a new file if not
                   
                outputFile = open(outputPath,append_write)
                
                S1_found = False  # flag to track if S1 has been found in file
               
                for ES_Line in readFile:
                    if searchParam2 in ES_Line and not S1_found:
                        ES_Line_Arr = ES_Line.split(' ')
                        print(ES_Line_Arr)
                        S1_Val = ES_Line_Arr[14:16]
                        print(f"S1_Val: {S1_Val}")
                        writer = csv.writer(outputFile)
                        if append_write == 'w':
                            writer.writerow([root, S1_Val])
                        else:
                            writer.writerow([root, S1_Val])
                        S1_found = True  # set the flag to True after writing to file

# The parent path of the directory containing the log files is defined
# The path of the output csv file is defined  
parentPath = '/Users/artemiswebster/Downloads/carbazole_structures'
outputPath = '/Users/artemiswebster/Downloads/carbazole_structures/S1_output.csv'

# The function is called
findS1Gauss(parentPath, outputPath)


# In[ ]:




