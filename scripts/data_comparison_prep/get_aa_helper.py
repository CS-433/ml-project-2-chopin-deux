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
from zander import *
from scipy.spatial import distance

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
    """Returns a dictionary of cartesian coordinates and features (both .npy-files)
     for a pdb-chain combination. Provide pdb_id_chain and directory containing surfaces"""
    coord_array = np.load(surf_path / str(pdb_id_chain +
                                          "_predcoords.npy"))
    feat_array = np.load(surf_path / str(pdb_id_chain +
                                         "_predfeatures.npy"))
    return {"xyz": coord_array, "feats": feat_array}

def surf_finder(pdb_chain, pdb_dir):
    """Finds surfaces for opposing chain of a given pdb. Provide the pdb_id and the chain
     of the pdb-file you would like to have the opposing surface for. Next, provide
     directory containing surfaces."""
    pdb_id = pdb_chain.split('_')[0]
    all_stem = pdb_dir.glob(pdb_id + "*.pdb")
    for i in all_stem:
        if i.stem != pdb_chain:
            surf_pdb_chain = i.stem
    return surf_pdb_chain

def get_ncac(residue):
    """"Takes residue class after loading structure via biopython. Returns N,
    Ca and C coordinates in numpy array"""
    try:
        ncac = np.stack([residue["N"].get_coord(), residue["CA"].get_coord(),
                         residue["C"].get_coord()])
    except KeyError:
        print("Missing backbone coordinate in: ", residue.get_resname(),
              residue.get_full_id()[3][1])
        ncac = None
    return ncac

def vectors_angle(v1, v2):
    """ Gives angle between two vectors. Gets unit vectors first.
    Max angle= 180 degrees (see cosine). Consider larger angles
    for side-chains that are flipped in."""
    v1_u= v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    return np.arccos(np.dot(v1_u,v2_u))

def extend_cb(residue):
    """If a residue has no C beta, extend to make one up, based on known angles."""
    N, CA, C = get_ncac(residue)
    return extend(C, N, CA, 1.522, 1.927, -2.143)

def get_cb_coord(residue, pdb_id, chain):
    """Returns C-beta coord of a residue object. Take pdb_id and chain in case
    you can not get Cb of non-Gly residue."""
    restype = residue.get_resname()
    if restype != 'GLY':
        try:
            cb_coord = residue['CB'].get_coord()
        except:
            # Might still not work, even if  not Gly. Find a way to get CB in this edge case.
            print("Not Gly and still no cb, extend")
            print(residue.get_full_id()[3][1], residue.get_resname())
            print("pdb_id and chain", pdb_id, chain)
            cb_coord = extend_cb(residue)
    elif restype == 'GLY':
        # extend coordinates of CB
        cb_coord = extend_cb(residue)
    else:
        print("Residue is not glycine, still no C-beta.")
    ca_coord = residue['CA'].get_coord()
    vector = cb_coord - ca_coord
    return cb_coord, vector

def dist_selector(c_coord, surf_dict, n_points):
    """From pdb-chain surface, selects closest N points to a given c_coord.
    Returns the coordinates, distances, indices and features of these points."""
    point_coords = surf_dict['xyz']
    point_feats = surf_dict['feats']
    # Make an array of all distances of all surface points to the c-boord. Gives index after sorting.
    dist_array = distance.cdist([c_coord], point_coords,
                                'euclidean')[0]
    dist_index = np.argsort(dist_array)
    # Get distances of selected points to c-coord and their coordinates
    selected_coords = point_coords[dist_index][0:n_points]
    selected_dists = dist_array[dist_index][0:n_points]
    # Get features of selected points, based on index
    selected_feats = point_feats[dist_index][0:n_points]
    return selected_coords, selected_dists, dist_index, selected_feats

def get_n_points(res_obj, surf_dict, pdb_id, chain):
    """Returns features of the nearest surface points from an opposing protein, for a given residue.
    Also adds features such as distance and angles."""
    # res_obj is a biopython class. Retrieved with: structure[0][chain][k].
    # get coordinates of CB atom.
    n_points = 50
    cb_coord, c_vector = get_cb_coord(res_obj, pdb_id, chain)
    selected_coords, selected_dists, dist_index, selected_feats = dist_selector(
        cb_coord, surf_dict, n_points)
    # Get carbon vector angles. First make an array of vectors between carbon beta and selected points.
    point_c_vectors = selected_coords - cb_coord
    angles = np.array([vectors_angle(c_vector, i) for i in point_c_vectors])
    # Get backbone atoms for theta angles (similar to backbone nitrogen)
    ncac = get_ncac(res_obj)
    N, CA, C = ncac[0], ncac[1], ncac[2]
    thetas = np.array([to_dih(N, CA, cb_coord, i) for i in selected_coords])
    final_feats = np.column_stack((selected_feats, selected_dists, angles, thetas))
    # =============================================================================
    # # UNCOMMENT this if you want to save point patch per AA, for visualisation in pymol
    # seq = res_obj.get_full_id()[3][1]
    # AA_type = res_obj.get_resname()
    # point_dir = Path("../../data/aa_points")
    # np.save(point_dir / (pdb_id  + chain + AA_type + str(seq) + "_predcoords.npy"), selected_coords)
    # np.save(point_dir / (pdb_id + chain + AA_type + str(seq) + "_predfeatures.npy"), selected_feats)
    # =============================================================================
    return final_feats

def surface_interface(surf_dict):
    """Takes dict of surface points coordinates and features and only returns
    features with interface feature == True"""
    interface_idx = np.where(surf_dict['feats'][:, 33] == 1)[0]
    interface_coords = surf_dict['xyz'][interface_idx]
    return interface_coords

def get_ca_coord(resi):
    """Gets C-alpha coordinates for residue"""
    return resi['CA'].get_coord()

def find_inter_res(res_ls, inter_coords):
    """Returns residues from a protein structure (chain) that are close to opposing chain.
    Takes list of residues on one side and dict of interface points on the other chain as input. """
    dist_thres = 4
    ca_arr = np.stack(list(map(get_ca_coord, res_ls)))
    # Compute distances of all protein carbon alphas to interface points. Only keep residues within
    # dist threshold. Use this index for residue selection.
    dist_bool = np.any(distance.cdist(ca_arr, inter_coords, 'euclidean') < dist_thres, axis = 1)
    close_res = np.array(res_ls)[dist_bool]
    return close_res