"""
Created on 10.12.21
Load data for the script that prepares dataset by getting specific residues based on dataframe.
@author: maxjansen
"""

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

def load_directories():
    """Loads the necessary directories to get residues for dataset preparation. Directories are:
    General data, surfaces and pdbs."""
    # Get the pandas df containing all the pdb's you need.
    data_dir = Path('../../data/')
    surf_dir = Path('../../data/surfaces')
    pdb_dir = Path('../../data/rosetta_pdbs/')
    return data_dir, surf_dir, pdb_dir

def prepare_df(data_dir):
    """Prepares the DataFrame containing residues, pdb's they belong to, chains they belong to, etc. For
    further use."""
    df = pd.read_csv(data_dir / "rosetta_res.csv")
    df = df.drop(['Unnamed: 0'], axis=1)

    # Sort dataframe by pdb and chain to prepare for efficient iteration
    df = df.sort_values(by=['pdb', 'chain', 'resi'])
    df = df.reset_index()
    df = df.drop(['index'], axis=1)

    pd.set_option('display.max_rows', df.shape[0])

    # Iterate through dataframe. First get unique combinations of pdb_id and chains that residues belong to.
    pdb_chain_list = list(df.groupby(['pdb', 'chain']).count().index)
    return df, pdb_chain_list

def df_cleaner(df, pdb_dir, surf_dir):
    """Remove rows of the df so you only keep pdbs where you can retrieve every single residue."""
    # Keep count of how many pdb's in df got the right residues
    counter = 0
    counter_no = 0
    clean_df = df
    pdb_chain_list = list(df.groupby(['pdb', 'chain']).count().index)
    # First iteration is just for checking whether everything matches
    for i in tqdm(pdb_chain_list):
        # Iterate through unique pdb-chain combinations in DataFrame.
        pdb_id = i[0]
        chain = i[1]

        # Find indexes of residues in dataframe for pdb-chain combo. Get residue positions and types.
        res_idxs = (df['pdb'] == pdb_id) & (df['chain'] == chain)
        resids = list(df[res_idxs]["resi"])
        restypes = list(df[res_idxs]["restype"])

        try:
            # Get pdb-filepaths in pdb directory that match df "pdb_id + chain" combinations with wildcard.
            # IMPORTANT WARNING: Do this because df 'chains' column only contain single chain_id's,
            # whereas some pdb-files have multiple chains in the filename. Hence the asterisks.
            pdb_fpath = next(pdb_dir.glob(str("*" + (pdb_id) + "*" + str(chain) + "*")))
            # Load the entire pdb_chain
            structure = struct_loader(pdb_fpath)
            opp_chain = surf_finder(pdb_fpath.stem, pdb_dir)
            surf_dict = load_surface_np(opp_chain, surf_dir)

            for j,k in enumerate(resids):
                # Check whether residue in df matches residue in structure at position!
                if restypes[j].upper() == structure[0][chain][k].get_resname():
                    counter += 1
                elif restypes[j].upper() != structure[0][chain][k].get_resname():
                    counter_no += 1
                    # After check, remove this one from df to ensure df indexes correspond to target npy-files
                    print("Warning, mismatch in df vs. pdb-file")
                    print(pdb_fpath.stem, k, restypes[j], structure[0][chain][k].get_resname())
                    clean_df.drop(res_idxs[res_idxs].index.to_list(), inplace=True)
        except:
            counter_no += 1
            clean_df.drop(res_idxs[res_idxs].index.to_list(), inplace=True)

            print("pdb not found")
    print("counter ", counter)
    print("counter_no", counter_no)
    return clean_df