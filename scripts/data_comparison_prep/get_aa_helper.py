"""
Created on 09.12.21
Helper to get the specific amino acids that Andreas used for the RosettaSurf comparison
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

res2num = {"ALA": 0, "ASX": 1, "CYS": 2, "ASP": 3, "GLU": 4, "PHE": 5, "GLY": 6, "HIS": 7,
           "ILE": 8, "LYS": 9, "LEU": 10, "MET": 11, "ASN": 12, "PRO": 13, "GLN": 14,
           "ARG": 15, "SER": 16, "THR": 17, "SEC": 18, "VAL": 19, "TRP": 20, "XAA": 21,
           "TYR": 22, "GLX": 23}

def res_loader(fname):
    """Loads pdb and returns residues in a useful format."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", fname)
    residues = structure.get_residues()
    return residues

def struct_loader(fname):
    """Loads pdb in filepath, returns structure"""
    parser = PDBParser(QUIET=True)
    return parser.get_structure("structure", fname)

def load_surface_np(pdb_id_chain, surf_path):
    """Returns a dictionary of cartesian coordinates and features (both .npy-files) for a pdb-chain combination.
    Provide pdb_id_chain and directory containing surfaces"""
    coord_array = np.load(surf_path / str(pdb_id_chain +
                                          "_predcoords.npy"))
    feat_array = np.load(surf_path / str(pdb_id_chain +
                                         "_predfeatures.npy"))
    return {"xyz": coord_array, "feats": feat_array}

def surf_finder(pdb_chain, pdb_dir):
    """Finds surfaces for opposing chain of a given pdb. Provide the pdb_id and the chain
     of the pdb-file you would like to have the opposing surface for. Next, provide directory containing surfaces."""
    pdb_id = pdb_chain.split('_')[0]
    all_stem = pdb_dir.glob(pdb_id + "*.pdb")
    for i in all_stem:
        if i.stem != pdb_chain:
            surf_pdb_chain = i.stem
    return surf_pdb_chain

# def get_res_points()
