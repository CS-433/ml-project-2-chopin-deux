"""
Created on 03.12.21
This script allows you to only get the relevant pdbs from the cluster. Does this based on the surface data,
@author: maxjansen
"""

from pathlib import Path
import re
import glob
from itertools import chain
import pandas as pd

# Get dMaSIF surface prediction data, make a set of all pdb_id's for which you have a surface
surf_list = []
surf_dir= Path('../data/Rosetta_surfaces')
for filepath in surf_dir.glob("*predcoords.npy"):
     pdb_id = str(filepath).split("_pred")[0]
     pdb_id = str(pdb_id).split("surfaces/")[1]
     surf_list.append(pdb_id)

surf_set = list(sorted(set(surf_list)))

with open('../data/rosetta_pdb_names.txt', 'w') as f:
     for item in surf_set:
          item += ".pdb"
          f.write("%s\n" % item)

print("Length of surf list:", len(surf_list))
print("Length of surf set:", len(surf_set))
print("Number of pdbs:", len(surf_set)/2)
