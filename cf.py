#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue, 13/2/2024

@author: pappel

the latest incarnation CalcFormula in python
includes:
    1) general discrimination of the data based on a decision schema build on compositional criteria
    2) calculation of dataframes for the phases with their stoiciometric coefficients 
    3) estimation of the ferrous iron contents based on the Droop algorithm (includes back-calculation of oxides, completely according to Droop)
    4) plotting of some key diagrams for the phases

  ********** ToDo ************
1) Input data (wt-% from ems) can be stored in text files with file names of the sample (e.g. T56-1.txt). Always one file per sample.
When reading the files, a list of file names can be created, which is later used for the markers of the plots.




"""
import pandas as pd
import stoic as st
import os
import minplots as mp


# a dictionary that relates the minerals with their oxygen base
ox4min = {'gt':12, 'fsp':8, 'px':6, 'amph':23}

# simple python list for storing the dicts of formulas
ferrousform = list()


# example here: 3 datafiles - they will be concatenated to one. For this, create a list of dataframes ('frames')
# ToDo:   create a list of all datafiles in a data directory, read them all and create a large datafile from them

# Define the dir path
directory = './in/'
# get a list of all files in this directory
files = [f for f in os.listdir(directory) if f.endswith('txt')]

# Create an empty list to store DataFrames
dataframes = []

"""
 This section is only to remove the first 3 and the last 7 lines from ems files
and overwrite the original files with the modified. 

# Loop through the files and read each into a DataFrame
for file in files:
    file_path = os.path.join(directory, file)

    # Get the :total number of lines in the file
    with open(file_path, 'r') as f:
        total_lines = sum(1 for _ in f)
    
    # Calculate the number of rows to read
    nrows = total_lines - 7  # Total lines minus last 7 lines
    
    try:
        # Read the file while skipping the first 2 lines and ignoring the last 4 lines
        df = pd.read_csv(file_path, skiprows=3, nrows=nrows - 3)  # Subtract 2 more for skipped lines
        df.to_csv(file_path, index=False)

    except Exception as e:
        print(f"Error reading {file}: {e}")

"""

# Read and combine all files in the directory into one DataFrame using pd.concat. Field separator is tabstop \t
combined_df = pd.concat((pd.read_csv(os.path.join(directory, file), engine='python', sep='\s*\\t\s*') for file in files),  ignore_index=True)

# init some empty python lists for later use
gt_frmls = [] 
fsp_frmls = [] 
cpx_frmls = list()
glc_frmls = list()
amph_frmls = list()
amph_sites = []
amph_ferricf = []

# This is to replace nan fields ('Not a number') with zero
combined_df.fillna(0, inplace=True)
# ... and insert a new column on second place for a phase id
combined_df.insert(1, "min", "")


# these are the conditions to assigne each measurement a phase id, based on compositional ranges (for pyroxene not correct)

combined_df.loc[(35 <= combined_df.SiO2) & (41 >= combined_df.SiO2) & (combined_df.Al2O3 > 5), 'min'] = "gt" 
combined_df.loc[(42 <= combined_df.SiO2) & (70 >= combined_df.SiO2) & (combined_df.Al2O3 > 18) & ((combined_df.Na2O + combined_df.K2O + combined_df.CaO) > 10), 'min'] = "fsp" 
combined_df.loc[((combined_df.FeO + combined_df.MgO) > 35) & (58 > combined_df.SiO2) & (combined_df.SiO2 > 48) & (combined_df.Al2O3 < 14), 'min'] = "opx" 
combined_df.loc[(combined_df.SiO2 > 72) & (combined_df.Al2O3 > 8) & (combined_df.K2O <7), 'min'] = "unk"
combined_df.loc[(combined_df.Total > 98.6) & (combined_df.SiO2 < 43) & (combined_df.Al2O3 < 1) & (combined_df.CaO < 1) & ((combined_df.FeO + combined_df.MgO) > 58), 'min'] = "ol"
combined_df.loc[(combined_df.SiO2 > 72) & (combined_df.Al2O3 > 8) & (combined_df.K2O <7), 'min'] = "chl"
combined_df.loc[(combined_df.SiO2 < 52) & (combined_df.Al2O3 > 29) & (combined_df.K2O < 0.9) & (combined_df.CaO < 1), 'min'] = "crd"


combined_df.loc[(combined_df.SiO2 <= 52) & (combined_df.Al2O3 > 9) & (combined_df.K2O > 7), 'min'] = "ms"

combined_df.loc[((combined_df.CaO + combined_df.Na2O) > 13) & (combined_df.SiO2 > 46) & (combined_df.SiO2 < 61) & (combined_df.Al2O3 < 25), 'min'] = "cpx"
combined_df.loc[(combined_df.Total < 80), 'min'] = "fail"
combined_df.loc[(combined_df['SiO2'] < 60) & (combined_df['SiO2'] > 48) & (combined_df['Total'] < 99.6) & (combined_df.Al2O3 < 14) & (combined_df.Na2O > 5), 'min'] = "amph"
# combined_df.loc[(combined_df.SiO2 > 49) & (combined_df.Al2O3 < 14) & (combined_df.K2O < 2) & (combined_df.Na2O > 5), 'min'] = "glc"
combined_df.to_csv('./out/alldata.txt',sep='\t', index=False)

gt_df = combined_df[combined_df['min'] == "gt"]         # extract a separate dataframe containing only garnet
fsp_df = combined_df[combined_df['min'] == "fsp"]       # ... the same here for fsp etc...
cpx_df = combined_df[combined_df['min'] == "cpx"]
amph_df = combined_df[combined_df['min'] == "amph"]
glc_df = combined_df[combined_df['min'] == "glc"]

if len(gt_df) > 0:      # check if there are data available
    for index, row in gt_df.iterrows():
            frml = st.minform(row.to_dict(), 12)
            frml['sample'] = str(row['Comment'])[:5]
            gt_frmls.append(frml)
    mp.gt_tern(pd.DataFrame.from_dict(gt_frmls), True)

if len(fsp_df) > 0:
    for index, row in fsp_df.iterrows():
        frml = st.minform(row.to_dict(), 8)
        fsp_frmls.append(frml)
    mp.fsp_tern(pd.DataFrame.from_dict(fsp_frmls), True)

    """
if len(glc_df) > 0:
    for index, row, in glc_df.iterrows():
        frml = st.minform(row.to_dict(), 23)
        frml_ferric = st.amph_ferric_minform(row.to_dict())
        glc_frmls.append(frml)
    """
if len(cpx_df) > 0:
    for index, row in cpx_df.iterrows():
        frml = st.minform(row.to_dict(), 6)
        cpx_frmls.append(frml)
    mp.px_tern(pd.DataFrame.from_dict(cpx_frmls), True)


if len(amph_df) > 0:
    for index, row in amph_df.iterrows():
        frml = st.minform(row.to_dict(), 23)
        frml_ferric, site_frml = st.amph_ferric_minform(row.to_dict())
        frml_ferric['sample'] = str(row['Comment'])[:5]
        frml['sample'] = str(row['Comment'])[:5]
        amph_frmls.append(frml)
        amph_ferricf.append(frml_ferric)
        amph_sites.append(site_frml)

    frmls_df = pd.DataFrame(amph_frmls)
    af_df = pd.DataFrame(amph_ferricf)
    as_df = pd.DataFrame(amph_sites)
    frmls_df.to_csv('./out/frmls.txt', sep= "\t", index=False)
    af_df.to_csv('./out/ferric.txt', sep= "\t", index=False)
    as_df.to_csv('./out/sites.txt', sep = "\t", index=False)
    mp.sodicAmphs(af_df, True)

        

