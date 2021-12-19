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
    df = pd.read_csv(data_dir / "cleaned_df.csv")
    pdb_chain_list = list(df.groupby(['pdb', 'chain']).count().index)
    print("Loaded previously cleaned df successfully!")
except:
    print("No cleaned df to load, create it from scratch.")
    df, pdb_chain_list = df_cleaner(df, pdb_chain_list, pdb_dir)
    df.to_csv(data_dir / "cleaned_df.csv")


# Create empty list of restypes, only append to this if restype in df matches pdb later on.
target_dict = {'target_type': []}
feat_dict = {'feat_data': []}

print("target dict", target_dict)
print("feat dict", feat_dict)
error_count = 0

# Return usable training and test set. Load pdbs in ascending order by pdb, chain, residue.
for i in tqdm(pdb_chain_list):
    # Iterate through unique pdb-chain combinations in DataFrame.
    pdb_id = i[0]
    chain = i[1]
    # Find indexes of residues in dataframe for pdb-chain combo.Then list all these residue positions and types.
    res_idxs = (df['pdb'] == pdb_id) & (df['chain'] == chain)
    resids, restypes = list(df[res_idxs]["resi"]), list(df[res_idxs]["restype"])

    # try:
    # Get pdb-filepaths in pdb directory that match df "pdb_id + chain" combinations with wildcard.
    # IMPORTANT WARNING: Do this way because df 'chains' column only contain single chain_id's,
    # whereas some pdb-files have more than two chains in the filename. Hence the asterisks.
    pdb_fpath = next(pdb_dir.glob(str("*" + (pdb_id) + "*" + str(chain) + "*")))
    # Load the entire pdb_chain
    structure = struct_loader(pdb_fpath)
    # Get opposing chain and load its surface
    opp_chain = surf_finder(pdb_fpath.stem, pdb_dir)
    surf_dict = load_surface_np(opp_chain, surf_dir)
    print("pdb_fpath", pdb_fpath)
    print(resids, restypes)
    print("opp_chain", opp_chain)

    for j,k in enumerate(resids):
        # Check whether residue in df matches residue in structure at position!
        res_obj = structure[0][chain][k]
        if restypes[j].upper() == res_obj.get_resname():
            target_dict['target_type'].append(res2num[restypes[j].upper()])
            # Get Cb coordinate, then nearest points and features.
            near_points = get_n_points(res_obj, surf_dict, pdb_id, chain)
            print("near dists", near_points[2])
            # Append to feat_data
            #feat_dict['feat_data']
        else:
            error_count += 1
            print("Warning! pdb-chain-res mismatch in df vs. pdb-file. Fix this!")
            print(pdb_fpath.stem, k, restypes[j], structure[0][chain][k].get_resname())

    # except:
    #     error_count += 1
    #     print(i)
    #     print("pdb not found")

    # HERE get all other amino acids at interface

print("counter_no", error_count)
print(target_dict)

print(df.columns)



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