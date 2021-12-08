"""
Created on 26.11.21

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


# Get the pandas df containing all the pdb's you need.
data_dir = Path('../data/')
surf_dir = Path('../data/surfaces')
pdb_dir = Path('../data/pdbs/')
df = pd.read_csv(data_dir / "rosetta_res.csv")


#iterate through dataframe
df = df.sort_values(by=['pdb_id'])
df = df.reset_index()
df = df.drop(['index'], axis=1)

def aa_selector(pdb_path, chain, resi):
    """Loads pdb and returns selected residue"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)
    resi = int(resi)
    residue = structure[0][chain][resi]
    return residue

def load_surface_np(pdb_id_chain):
    surf_path = Path('/Users/maxjansen/EPFL/AA_at_Interface/dataset/surfaces')
    coord_array = np.load(surf_path / str(pdb_id_chain +
                                          "_predcoords.npy"))
    feat_array = np.load(surf_path / str(pdb_id_chain +
                                         "_predfeatures.npy"))
    return {"xyz": coord_array, "feats": feat_array}

for i in df.index:
    # Get pdb and chain for iteration
    pdb_stem = df['pdb_id'][i].split('_')[0]
    pdb_chain = df['chain'][i]
    resi = df['resi'][i]
    restype = df['restype'][i]

    # Look for pdb_files containing your selected chain in directory
    pdb_name = pdb_dir.glob(pdb_stem + '*')
    for j in pdb_name:
        chains = str(j).split('_')[1]
        pdb_chain = str(pdb_chain)
        # Only keep files if pdb chain matches chain from df row(i)
        if (pdb_chain in chains) == True:
            aa = aa_selector(j, pdb_chain, resi)
            print("res chain " + str(j))
        else:
            print("surf chain" + str(j))





    # Get relevant surface

# Load pdb
#def pdb_load:

# Function to extract AA from pdb

# Load the pdb you need.

# Within a given pdb, make a function to get specific residue based on ID

#