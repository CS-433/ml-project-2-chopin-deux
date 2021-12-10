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
from tqdm import tqdm
from scipy.spatial import distance
from zander import *




# Get the pandas df containing all the pdb's you need.
data_dir = Path('../../data/')
surf_dir = Path('../../data/rosetta_surfaces')
pdb_dir = Path('../../data/rosetta_pdbs/')
df = pd.read_csv(data_dir / "rosetta_res.csv")
df = df.drop(['Unnamed: 0'], axis=1)


# Sort dataframe by pdb and chain to prepare for efficient iteration
df = df.sort_values(by=['pdb', 'chain'])
df = df.reset_index()
df = df.drop(['index'], axis=1)
pd.set_option('display.max_rows', df.shape[0])

# Iterate through dataframe. First get unique combinations of pdb_id and chains that residues belong to.
pdb_chain_list = list(df.groupby(['pdb', 'chain']).count().index)

# Count for how many pdb's in df you got the right residues
counter = 0
counter_no = 0


for i in tqdm(pdb_chain_list[2:3]):
    # Get pdb-filepaths in pdb directory that match df "pdb_id + chain" combinations.
    # Do this because df 'chains' column only contain single chain_id's,
    # whereas some pdb-files have multiple chains in the filename. Hence the asterisks.
    pdb_id = i[0]
    chain = i[1]
    res_idxs = (df['pdb'] == pdb_id) & (df['chain'] == chain)
    resids = list(df[res_idxs]["resi"])
    restypes = list(df[res_idxs]["restype"])

    try:
        pdb_fp = pdb_dir.glob(str("*" + (pdb_id) + "*" + str(chain) + "*"))
        # Get surface of opposite chain
        surf_fp = list(surf_dir.glob(str("*" + pdb_id + "*" )))
        # (pdb_id + chain) not in

        print(surf_fp[0])
        counter += 1
        structure = struct_loader(next(pdb_fp))
        for j in resids:
            print(structure[0][chain][j])


        # # Load one pdb_filepath at a time, returns generator of residues
        # resi_gen = res_loader(next(pdb_fp))
        # print(list(resi_gen))

    # From all residues
    except:
        print('No pdb found')
        counter_no += 1
print("counter ", counter)
print("counter_no", counter_no)
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


# def load_surface_np(pdb_id_chain):
#     surf_path = Path('/Users/maxjansen/EPFL/AA_at_Interface/dataset/surfaces')
#     coord_array = np.load(surf_path / str(pdb_id_chain +
#                                           "_predcoords.npy"))
#     feat_array = np.load(surf_path / str(pdb_id_chain +
#                                          "_predfeatures.npy"))
#     return {"xyz": coord_array, "feats": feat_array}

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