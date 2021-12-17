"""
Created on 26.11.21
This script creates a test and training set to compare to RosettaSurf (bound) as a benchmark.
The test set is what we will perform predictions on, the ~200 * 20 residues from 1495 pdb's.
All other residues from those pdb's will serve as a training set.
@author: maxjansen
"""

# Get the specific amino acids that Andreas used for the RosettaSurf comparison
# Import packages
from pathlib import Path
import re
import glob
from itertools import chain
import pandas as pd
import numpy as np
from Bio.PDB import *
from get_aa_helper import *
from get_aa_loader import *
from tqdm import tqdm
from scipy.spatial import distance
from zander import *

# Load the directories containing general data, surface-npy-files, and pdbs.
data_dir, surf_dir, pdb_dir = load_directories()
# Load DataFrame of residues and pdbs
df, pdb_chain_list = prepare_df(data_dir)

# Prepare a new df that will contain only verified rows, useful for checking npy later
try:
    clean_df = pd.read_csv(data_dir / "cleaned_df.csv")
    print("Loaded previously cleaned df successfully!")
except:
    print("No cleaned df to load, create it from scratch.")
    clean_df = df_cleaner(df, pdb_chain_list, pdb_dir)
    clean_df.to_csv(data_dir / "cleaned_df.csv")


# Create empty list of restypes, only append to this if restype in df matches pdb later on.
target_dict = {'Target type': []}
feat_dict = dict.fromkeys(['feat_data'])



#surf_dict = load_surface_np(opp_chain, surf_dir)

#####
#target_dict['Target type'].append(res2num[restypes[j].upper()])
# Get C-beta coordinate with function. Uses pseudo coordinate for Glycine
#print(structure[0][chain][k]["CB"])

# Get opposing chain and load its surface
# opp_chain = surf_finder(pdb_fpath.stem, pdb_dir)
######

print(clean_df.columns)







    # except:
    #     print('No pdb found')
    #     counter_no += 1

### Important!!! Save test and training set here.


    # Load all residues from a pdb once
    #res_loader(pdb_filepath)

# Get the surface for the opposing chain



# Load the pdb you need.
# Get relevant surface
# Within a given pdb, make a function to get specific residue based on ID

#
##########################
# def aa_selector(pdb_path, chain, resi, list):
#     """Loads pdb and returns selected residue based on """
#     parser = PDBParser(QUIET=True)
#     structure = parser.get_structure("structure", pdb_path)
#     if list == True:
#
#         for i in resi:
#             resi = int(resi)
#             residue = structure[0][chain][resi]
#     return residue
##############################




# for i in df.index:
#     # Get pdb and chain for iteration
#     pdb_stem = df['pdb_id'][i].split('_')[0]
#     pdb_chain = df['chain'][i]
#     resi = df['resi'][i]
#     restype = df['restype'][i]
#
#     # Look for pdb_files containing your selected chain in directory
#     pdb_name = pdb_dir.glob(pdb_stem + '*')
#     for j in pdb_name:
#         chains = str(j).split('_')[1]
#         pdb_chain = str(pdb_chain)
#         # Only keep files if pdb chain matches chain from df row(i)
#         if (pdb_chain in chains) == True:
#             aa = aa_selector(j, pdb_chain, resi)
#             print("res chain " + str(j))
#         else:
#             print("surf chain" + str(j))