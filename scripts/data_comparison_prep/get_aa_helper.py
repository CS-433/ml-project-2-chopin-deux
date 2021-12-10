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