"""
Created on 03.12.21

@author: maxjansen
"""
### This script compares the pdb_id's I put in a list to replace the dMaSIF test set to the actual surfaces I got back.
### It seems like there is a discrepancy, and I wonder why dMaSIF does not process some pdbs's.

from pathlib import Path
import re
import glob
from itertools import chain
import pandas as pd

# Get all files in the surfaces directory and chop off the chain, description and filetype.
# Get RosettaSurf comparison residues and the pdb they came from
data_dir= Path('../data/surfaces')
surface_pdb_list = []
for filepath in data_dir.glob("*predcoords.npy"):
    surf_pdb_id = str(filepath).split("_")[0]
    surf_pdb_id = surf_pdb_id.split("surfaces/")[1]
    surface_pdb_list.append(surf_pdb_id)
# Turn this list into a set
print("Set of processed pdb surfaces: " , len(set(surface_pdb_list)))

# Get the list that you originally used to run dMaSIF and turn this into a set.
# Get RosettaSurf comparison residues and the pdb they came from
data_dir= '../data/'
rosetta_pdb_list = data_dir + "rosetta_surf_select.txt"
ros_f = open(rosetta_pdb_list, "r")
pdb_og_list = []
for line in ros_f.readlines():
    pdb_id = line.split("_")[0]
    pdb_og_list.append(pdb_id)

print("Set of original rosettasurf surfaces: ", len(set(pdb_og_list)))
print("Difference between the two sets: ", len(set(pdb_og_list)) - len(set(surface_pdb_list)))
print("Intersection of rosettasurf and processed surfaces: ", len(set(pdb_og_list).intersection(set(surface_pdb_list))))
print("List of rosettasurf pdbs that were not processed to surfaces: ",
      list(sorted(set(pdb_og_list).difference(set(surface_pdb_list)))))

no_surface_ls = list(sorted(set(pdb_og_list).difference(set(surface_pdb_list))))

textfile = open("../data/no_pred_surface.txt", "w")
for element in no_surface_ls:
    textfile.write(element + "\n")
textfile.close()

# for filepath in data_dir.glob("*.slurm"):
#     aa_type = str(filepath).split(".slurm")[0][-3:]
#     aa_series = pd.Series(aa_type, name = 'restype')
#     aa_series = aa_series.repeat(200)
#     aa_series = aa_series.reset_index()
#     aa_series = aa_series.drop(['index'], axis=1)
#     f = open(filepath, "r")
#     for line in f.readlines():
