# Pablo Gainza Cirauqui 2016 LPDI IBI STI EPFL
# This pymol plugin for Masif just enables the load ply functions.

#sys.path.append('C:\\Users\\antho\\OneDrive\\Documents\\Pymol\\masif_pymol_plugin\\')

from pymol import cmd
from .loadPLY import *
from .loadDOTS import *
from .next_seed import *
#from .dmasif_pymol import *
import sys


cmd.extend('loadply', load_ply)
cmd.extend('loaddots', load_dots)
cmd.extend('loadgiface', load_giface)
cmd.extend('n', n_)
cmd.extend('next_', next_)
cmd.extend('next', next_)
cmd.extend('next_seed', next_seed)

#cmd.extend("loadvtk", load_vtk)
#cmd.extend("loadpred", load_pred)
