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
    df = df_cleaner(df, pdb_dir, surf_dir)
    pdb_chain_list = list(df.groupby(['pdb', 'chain']).count().index)
    df.to_csv(data_dir / "cleaned_df.csv")


# Create empty list of restypes, only append to this if restype in df matches pdb later on.
#IMPORTANT
test_dict = {'target_type': [], 'feats': []}
train_dict = {'target_type': [], 'feats': []}

error_count = 0
print(len(pdb_chain_list))
# Return usable training and test set. Load pdbs in ascending order by pdb, chain, residue.
for i in tqdm(pdb_chain_list):
    # Iterate through unique pdb-chain combinations in DataFrame.
    pdb_id = i[0]
    chain = i[1]
    # Find indexes of residues in dataframe for pdb-chain combo.Then list all these residue positions and types.
    res_idxs = (df['pdb'] == pdb_id) & (df['chain'] == chain)
    resids, restypes = list(df[res_idxs]["resi"]), list(df[res_idxs]["restype"])
    print("test resids", resids)

    try:
        # Get pdb-filepaths in pdb directory that match df "pdb_id + chain" combinations with wildcard.
        # IMPORTANT WARNING: Do this way because df 'chains' column only contain single chain_id's,
        # whereas some pdb-files have more than two chains in the filename. Hence the asterisks.
        pdb_fpath = next(pdb_dir.glob(str("*" + (pdb_id) + "*" + str(chain) + "*")))
        # Load the entire pdb_chain
        structure = struct_loader(pdb_fpath)
        # Get opposing chain and load its surface
        opp_chain = surf_finder(pdb_fpath.stem, pdb_dir)
        surf_dict = load_surface_np(opp_chain, surf_dir)
        for j,k in enumerate(resids):
            # Check whether residue in df matches residue in structure at position!
            res_obj = structure[0][chain][k]
            if restypes[j].upper() == res_obj.get_resname():
                test_dict['target_type'].append(res2num[restypes[j].upper()])
                # Get Cb coordinate, then nearest points and features.
                near_points = get_n_points(res_obj, surf_dict, pdb_id, chain)
                # Append to feat_data
                test_dict['feats'].append(near_points)
            # Warn and count mismatch between df and pdb-file.
            else:
                error_count += 1
                print("Warning! pdb-chain-res mismatch in df vs. pdb-file. Fix this!")
                print(pdb_fpath.stem, k, restypes[j], structure[0][chain][k].get_resname())

        # Make training set, take all residues from loaded structure (at interface!) that are not in test
        inter_coords = surface_interface(surf_dict)
        tr_residues = structure.get_residues()
        tr_res_ls = []
        for resi in tr_residues:
            if resi.get_full_id()[3][1] not in resids:
                tr_res_ls.append(resi)
            else:
                pass
        close_res = find_inter_res(tr_res_ls, inter_coords)
        #Get N nearest points to residue CB and save
        for res_obj in close_res:
            train_dict['target_type'].append(res2num[res_obj.get_resname()])
            train_dict['feats'].append(get_n_points(res_obj, surf_dict, pdb_id, chain))
    # Keep going if anything fails, but warn where it happened.
    except:
        error_count += 1
        print(i)
        print("pdb not found")

print("error counter_no", error_count)


# Save training and test sets.
# Convert the 2 lists in `test_dict` to np arrays
test_dict['target_type'] = np.stack(test_dict['target_type'])
test_dict['feats'] = np.stack(test_dict['feats'])
# Convert target types to one-hot encoded entries
te_types_array = np.zeros((len(test_dict['target_type']), len(res2num)))
for i, t in enumerate(test_dict['target_type']):
        te_types_array[i, t] = 1.0

print("length test check: ", len(test_dict['target_type']), len(test_dict['feats']), len(te_types_array))

# Convert the 2 lists in `train_dict` to np arrays
train_dict['target_type'] = np.stack(train_dict['target_type'])
train_dict['feats'] = np.stack(train_dict['feats'])
# Convert target types to one-hot encoded entries
tr_types_array = np.zeros((len(train_dict['target_type']), len(res2num)))
for i, t in enumerate(train_dict['target_type']):
        tr_types_array[i, t] = 1.0

print("length train check: ", len(train_dict['target_type']), len(train_dict['feats']), len(tr_types_array))
print(df.columns)



### Important!!! Save test and training set here.
np.save("../../data/dataset/test_feats.npy", test_dict['feats'])
np.save("../../data/dataset/test_target.npy", te_types_array)

np.save("../../data/dataset/train_feats.npy", train_dict['feats'])
np.save("../../data/dataset/train_target.npy", tr_types_array)
