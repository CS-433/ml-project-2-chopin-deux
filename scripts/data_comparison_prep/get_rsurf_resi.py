"""
Created on 25.11.21
This script returns a csv of all residues used by Andreas in RosettaSurf. It reads his slurm commands to get these.
'rosetta_res.csv' is complete and contains 200 of each of the 20 AA's. Also columns for pdb_id, aa type and position.
It also returns a csv of the residues the we managed to get a surface for, using dMaSIF. 99.6% were retrieved.

Main output is done: '../data/rosetta_res.csv'
@author: maxjansen
"""
from pathlib import Path
import re
import glob
from itertools import chain
import pandas as pd

# Get RosettaSurf comparison residues and the pdb they came from
data_dir= Path('../../data/')
df = pd.DataFrame(columns=['pdb_id', 'res_id'])
for filepath in data_dir.glob("rosetta_slurm/*.slurm"):
    aa_type = str(filepath).split(".slurm")[0][-3:]
    aa_series = pd.Series(aa_type, name = 'restype')
    aa_series = aa_series.repeat(200)
    aa_series = aa_series.reset_index()
    aa_series = aa_series.drop(['index'], axis=1)
    f = open(filepath, "r")
    for line in f.readlines():
        if line.startswith("pdb_id"):
            line = line.split("pdb_id=(")[1]
            line = line.split(")")[0]
            line = re.sub('"', '', line)
            line = list(line.split(" "))
            pdb_series = pd.Series(line, name="pdb_id")
        elif line.startswith("mut_res"):
            line = line.split("mut_res=(")[1]
            line = line.split(")")[0]
            line = re.sub('"', '', line)
            line = list(line.split(" "))
            resi_series = pd.Series(line, name="res_id")

    # Save data from different residue types, resid, pdb. All in one pandas dataframe.
    res_df = pd.concat([pdb_series, resi_series, aa_series], axis=1)
    df = pd.concat([df, res_df])
df = df.reset_index()

# Touch up the pandas df to make it more useful
def splitAN(series):
    resi = re.split('(\d+)', series)[1]
    chain = re.split('(\d+)', series)[2]
    return resi, chain
split_resi = df['res_id'].apply(splitAN)
df['resi'], df['chain'] = split_resi.str
df = df.drop(['res_id'], axis=1)
df = df.drop(['index'], axis=1)

# Get dMaSIF surface prediction data, make a set of all pdb_id's for which you have a surface
surf_list = []
surf_dir= Path('../../data/rosetta_surfaces/')
for filepath in surf_dir.glob("*predcoords.npy"):
     print(filepath)
     pdb_id = str(filepath).split("rosetta_surfaces/")[1]
     pdb_id = str(pdb_id).split("_")[0]
     surf_list.append(pdb_id)

surf_set = sorted(set(surf_list))
print(surf_set)

# Make a set of all pdb_id's occurring in the DataFrame. DataFrame is directly based on the residues Andreas used for
# RosettaSurf.
rosetta_set = sorted(set(list(df['pdb_id'])))
print(rosetta_set)
# 1495 unique pdbs in Rosetta dataset
# 319 pdbs that occur in Rosetta dataset AND dmasif_surface set, used to keep those


# This is what we will use for dataset generation in Max's AA prediction
final_df = df
# Keep the pdb_id column, but copy a split where you only keep the first 4 characters before the underscore
final_df['pdb'] = final_df.pdb_id.str.split("_").str[0]
comparison_df = final_df

# Use this new 'pdb' column to only keep rows that occur in the surface dataset you got from dMaSIF
print("Do we have 20 x 200=4000 residues in df?: ", final_df['pdb'].sum())

print("rosettasurf residue counts: \n ", final_df['restype'].value_counts())
final_df.to_csv(data_dir / 'rosetta_res.csv')


# Use this new 'pdb' column to only keep rows that occur in the surface dataset you got from dMaSIF
print("How many of 4000 rosettasurf intersect with dmasif?: ", comparison_df['pdb'].isin(surf_set).sum())
comparison_df = comparison_df[comparison_df['pdb'].isin(surf_set)]

print("dmasif overlap rosettasurf residue counts: \n ", comparison_df['restype'].value_counts())
comparison_df.to_csv(data_dir / 'dmas_rosetta_res.csv')


