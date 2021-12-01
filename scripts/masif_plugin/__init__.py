# Pablo Gainza Cirauqui 2016 LPDI IBI STI EPFL
# This pymol plugin for Masif just enables the load ply functions.

from pymol import cmd
from .loadPLY import *
from .loadDOTS import *
from .simple_mesh import *
from .next_seed import *
from .dmasif_pymol import *
import sys

def __init_plugin__(app):
    cmd.extend('loadply', load_ply)
    cmd.extend('loaddots', load_dots)
    cmd.extend('simple_mesh', simple_mesh)
    cmd.extend('loadgiface', load_giface)
    cmd.extend('nextseed', next_seed)
    cmd.extend('next', next_)
    #cmd.extend("loadvtk", load_vtk)
    #cmd.extend("loadpred", load_pred)
