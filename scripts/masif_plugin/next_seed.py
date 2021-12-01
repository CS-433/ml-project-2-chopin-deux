# Julius Upmeier zu Belzen 2019 LPDI IBI STI EPFL
# This pymol plugin for Masif conveniently loads the pdb and ply files


"""
CAUTION: currently requires to be run from a directory in which the loadPLY and loadDOTS plugins are present

This plugin displays the results of a masif search ordered by their score


ids are generally the names of the respecitive directories

Usage:
======
cd to the directory with the masif outputs of a specific protein (child directories would be 'site_0, site_1, ...')

initialize the plugin for the protein and a specific site:
PyMOL> next <id of the protein>, <site>

bring up the next hit:
PyMOL> next

move backwards one hit:
PyMOL> back

move to a specific hit:
PyMOL> back <id of the hit>

save a note in the directory corresponding to a hit (this is displayed, when the hit is loaded again):
PyMOL> note 'note'

display all notes:
PyMOL> readnotes

naive selection of the interface based on a distance cutoff (default 3 A)
PyMOL> getiface <cutoff>

print residue numbers of a selection
PyMOL> printsele <selection>

reset:
PyMOL> reset_next
"""
import pymol
from pymol import cmd
from pymol import stored
import os
import sys
import shutil
sys.path.append('/Users/maxjansen/anaconda3/share/pymol/data/startup/masif_pymol_plugin/')

from .loadPLY import *
from .loadDOTS import *

class next_obj:
    def __init__(self, name, site, ranking_criteria=3, target_focus='pb', binder_focus='None', note_only=False, hide=None):
        print('Initialising for presenting the next proteins: {}, {}'.format(name, site))
        self.name = name
        self.site = site
        self.target_focus = target_focus
        self.binder_focus = binder_focus
        self.note_only = note_only

        self.target_pdb = name + '.pdb'
        self.target_ply = name + '.ply'

        self.already_seen = set()

        # show target:
        cmd.load(self.target_pdb)
        cmd.set('cartoon_side_chain_helper', 'on', name)
        cmd.show('sticks', name)
        #cmd.hide('sticks', name + ' and backbone')

        cmd.hide('(h.) and ' + name)
        if hide:
            cmd.select('hide', hide)
            cmd.hide('sticks', 'hide')
            cmd.hide('cartoon', 'hide')

        from .loadPLY import load_ply
        load_ply(self.target_ply)

        # hide everything except focus
        cmd.disable('*_' + self.target_ply)
        #cmd.disable(name)
        cmd.enable(self.target_focus + '_' + name)

        # show target dots:
        self.target_vert = os.path.join(site, 'target.vert')
        from .loadDOTS import load_dots
        load_dots(os.path.join(site, 'target.vert'), name='target')
        cmd.set_name('vert_' + site + '_target.vert', name + '.vert')
        cmd.disable(name + '.vert')

        cmd.group('target', ' and '.join([name, name + '.ply', name + '.vert']))

        # prioritize
        binders = [d for d in os.listdir(site) if os.path.isdir(os.path.join(site, d))]
        print('Found {} distinct potential binders'.format(len(binders)))

        self.data = {}
        for binder in binders:
            binder_vars = set([f.split('.')[0] for f in os.listdir(os.path.join(site, binder)) if os.path.isfile(os.path.join(site, binder, f)) \
                and 'patch' not in f])
            print(binder_vars)
            for binder_var in binder_vars:
                with open(os.path.join(site, binder, binder_var + '.score'), 'r') as ifile:
                    line = ifile.readline().split(',')
                    score = line[int(ranking_criteria)].split(':')[1]
                    clashes_ca = line[3].split(':')[1]
                    clashes = line[4].split(':')[1]
                    score = float(score)
                    desc_dist = float(line[5].split(':')[1])
                    print(line)
                self.data[binder_var] = {'score': score, 'desc_dist_score':desc_dist, 'binder': binder, 'clashes_ca': clashes_ca, 'clashes': clashes}
        self.prios = sorted(self.data.keys(), key=lambda x: -self.data[x]['score'])
        self.prio_pos = -1


    def next_(self, pdb_id=None):
        if pdb_id is not None:
            for ix, key in enumerate(self.prios):
                if key.startswith(pdb_id):
                    self.prio_pos = ix
                    binder_var = key
                    break
            if binder_var is None:
                return
        else:
            self.prio_pos += 1
            binder_var = self.prios[self.prio_pos]

        pdb_id = binder_var.split('_')[0]
        chain = binder_var.split('_')[1]

        site = self.site
        binder = self.data[binder_var]['binder']

        print('Loading {} with score {}, descriptor_distance score: {}, clashes_ca: {},  and clashes (heavy_atoms) {}'.format(binder_var, self.data[binder_var]['score'],
            self.data[binder_var]['desc_dist_score'], self.data[binder_var]['clashes_ca'], self.data[binder_var]['clashes']))
        p = os.path.join(site, binder, binder_var)

        # display notes
        if os.path.isfile(p + '.txt'):
            with open(p + '.txt', 'r') as ifile:
                print('Note:')
                print(ifile.read())
                print('\n')
        else:
            if self.note_only:
                print('No note.')
                self.prio_pos += 1
                return

        if '{}_{}'.format(pdb_id, chain) in self.already_seen:
            return
        else:
            self.already_seen.add('{}_{}'.format(pdb_id, chain))

        cmd.load(p + '.pdb')
        #cmd.set('cartoon_side_chain_helper', 'on', binder_var)
        #cmd.set('cartoon_transparency', '0.5', binder_var)
        cmd.show('sticks', binder_var)
        #cmd.hide('sticks', binder_var + ' and backbone')
        cmd.hide('(h.) and ' + binder_var)

        loadPLY.load_ply(p + '.ply')

        # hide everything except focus
        cmd.disable('*_' + binder_var + '.ply')
        cmd.enable(self.binder_focus + '_*_' + binder_var + '.ply')

        # show target dots:
        #loadDOTS.load_dots(p + '.vert', color='yellow', name=binder_var)
        loadPLY.load_ply(p + '_patch.ply')
        cmd.delete('mesh_' + site + '_' + binder + '_' + binder_var + '_patch.ply')
        cmd.set_name('vert_' + site + '_' + binder + '_' + binder_var + '_patch.ply', binder_var + '.vert')
        cmd.center(binder_var + '.vert')

        # fix grouping
        cmd.delete(site + '_' + binder + '_' + binder_var + '.ply')
        cmd.delete(site + '_' + binder + '_' + binder_var + '_patch.ply')
        cmd.group(binder_var + '.ply', '*_' + binder_var + '.ply' + ' and ' + binder_var + ' and ' + binder_var + '.vert')


    def save(self, binder_var=None, out_parent="/home/gainza/lpdi_fs/protein_designs/masif_seeds/by_seed/"):

        if not binder_var:
            binder_var = self.prios[self.prio_pos]

        out_dir = os.path.join(out_parent, self.name, self.site, binder_var)

        pdb_id = binder_var.split('_')[0]
        chain = binder_var.split('_')[1]

        site = self.site
        binder = self.data[binder_var]['binder']

        p = os.path.join(site, binder, binder_var)
        pdbfile = p + '.pdb'
        plyfile = p + '.ply'
        patchfile = p + '_patch.ply'
        scorefile = p +'.score'

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        shutil.copy(pdbfile, out_dir)
        shutil.copy(plyfile, out_dir)
        shutil.copy(patchfile, out_dir)
        shutil.copy(scorefile, out_dir)

        shutil.copy(self.target_pdb, out_dir)
        shutil.copy(self.target_ply, out_dir)
        shutil.copy(self.target_ply, out_dir)
        shutil.copy(self.target_vert, out_dir)

        site = self.site
        binder = self.data[binder_var]['binder']
        return binder_var, self.name+'_'+self.site+'_'+binder_var, out_dir

    def back(self, back_to=None):
        if back_to:
            try:
                back_to = int(back_to)
                self.prio_pos = back_to

            except:
                self.prio_pos = self.prios.index(back_to)
        else:
            self.prio_pos -= 1

    def note(self, note, binder_var=None):
        if not binder_var:
            binder_var = self.prios[self.prio_pos]
        # write file to dir
        p = os.path.join(self.site, self.data[binder_var]['binder'], binder_var + '.txt')
        print('Wrote to {}'.format(p))
        with open(p, 'a') as ofile:
            ofile.write(note)

    def readnotes(self):
        for binder_var in self.prios:
            p = os.path.join(self.site, self.data[binder_var]['binder'], binder_var + '.txt')
            if os.path.isfile(p):
                with open(p, 'r') as ifile:
                    print('Reading {} with score {}\n{}\n'.format(p, self.data[binder_var]['score'], ifile.read()))


    def getiface(self, cutoff):
        # make selection
        binder_var = self.prios[self.prio_pos]

        cmd.select('iface_{}'.format(binder_var), binder_var + ' within ' + cutoff + ' of ' + '(' + self.name + ' and not hide)')




# wrap everything
def next_(name=None, site=None, ranking_criteria=2, target_focus='pb', binder_focus='None', note_only=False, hide=None, remove=None):

    try:
        if name:
            if stored.next_obj.name:
                print('Name mismatch: {} != {}'.format(name, stored.next_obj.name))
                raise AttributeError
        if site:
            if stored.next_obj.site != site:
                print('Site mismatch: {} != {}'.format(site, stored.next_obj.site))
                raise AttributeError
#        stored.next_obj.next()

    except AttributeError:
        if not site:
            site = 'site_0'
        stored.next_obj = next_obj(name, site, ranking_criteria, target_focus, binder_focus, note_only, hide)

def n_(pdb_id=None):
    stored.next_obj.next(pdb_id)

def back(back_to=None):
    stored.next_obj.back(back_to)


def note(note, binder_var=None):
    stored.next_obj.note(note, binder_var=binder_var)


def readnotes():
    stored.next_obj.readnotes()


def getiface(cutoff='3'):
    stored.next_obj.getiface(cutoff)


def printsele(sele='sele'):
    stored.resis=[]
    cmd.iterate(sele, 'stored.resis.append(resi)')
    resis = sorted(set(stored.resis), key=lambda x: int(x))
    print(','.join(resis))

def design(sele='sele'):
    stored.resis=[]
    cmd.iterate(sele, 'stored.resis.append((chain, resi))')
#    resis = #sorted(set(stored.resis), key=lambda x: int(x))
    resis = set(stored.resis)#sorted(set(stored.resis), key=lambda x: int(x))
    outstring = "NATAA\n"
    outstring += "USE_INPUT_SC\n"
    outstring += "start\n"

    for tup in resis:
        outstring += "{} {} PIKAA ADEFGHIKLMNOQRSTVWY\n".format(tup[1], tup[0])

    outfile = open('resfile', 'w')
    outfile.write(outstring)
    outfile.close()
    print(outstring)


def reset():
    stored.next_obj=None

def save_seed(count=1, binder_var=None, out_parent="/home/gainza/lpdi_fs/protein_designs/masif_seeds/by_seed/"):
    count = int(count)
    if count == 1:
        hotspot_res=['hot']
        seed=['seed']
    else:
        hotspot_res=['hot1', 'hot2']
        seed=['seed1', 'seed2']
    binder_name, name, output_dir = stored.next_obj.save(binder_var=binder_var, out_parent=out_parent)
    template_dir = "/home/gainza/lpdi_fs/protein_designs/masif_seeds/template_scripts/"
    # Copy the target to a new object
    cmd.copy('context', stored.next_obj.name )
    # Rename the chain
    cmd.alter('context', 'chain=\'A\'')
    out_graft_dir = os.path.join(output_dir, 'graft')
    if not os.path.exists(out_graft_dir):
        os.makedirs(out_graft_dir)
    # Save the target to the target directory as 'context.pdb'
    cmd.save(os.path.join(out_graft_dir, 'context.pdb'), stored.next_obj.name)

    # Save the seed to the target directory as 'seed.pdb'
    if count == 1:
        cmd.create('myseed', '(seed)')
    else:
        cmd.create('myseed', '(seed1) or (seed2)')
        cmd.set('pdb_use_ter_records', 'on')

    cmd.alter('myseed', 'chain=\'B\'')

    cmd.save(os.path.join(out_graft_dir, 'seed.pdb'), 'myseed')

    # Iterate through seed, and figure out the index of each hotspot.
    if count == 1:
        stored.vals = []
        cmd.iterate('(seed)', 'stored.vals.append(resi)')
        resi_asint = [int(x) for x in stored.vals]
        minresi = min(resi_asint)

        cmd.select('myhotsel', '(hot) and name CA')
        stored.vals = []
        cmd.iterate('(myhotsel)', 'stored.vals.append(resi)')
        hot = [int(x)-minresi+1 for x in stored.vals]

        hotstring = ''
        for ix in hot[0:-1]:
            hotstring+='{}:'.format(ix)
        hotstring += '{}'.format(hot[-1])
        print(hotstring)
    # Iterate through seed, and figure out the index of each hotspot.
    else:
        stored.vals = []
        hotstring = ''
        for ix in range(len(seed)):
            cmd.iterate('(seed{})'.format(ix+1), 'stored.vals.append(resi)')
            resi_asint = [int(x) for x in stored.vals]
            minresi = min(resi_asint)

            cmd.select('myhotsel', '(hot{}) and name CA'.format(ix+1))
            stored.vals = []
            cmd.iterate('(myhotsel)', 'stored.vals.append(resi)')
            hot = [int(x)-minresi+1 for x in stored.vals]

            for ix in hot[0:-1]:
                hotstring+='{}:'.format(ix)
            hotstring += '{}'.format(hot[-1])
            if ix == 0:
                hotstring += ', '
            print(hotstring)


    # Open the motif graft script and modify the hotspots
    target_file = open(os.path.join(template_dir, 'grafting.xml'), 'r')
    out_target_file = open(os.path.join(out_graft_dir, 'grafting.xml'), 'w')
    all_lines = target_file.readlines()
    for line in all_lines:
        if 'HOTSPOT_RES' in line:
            new_line = line.replace('HOTSPOT_RES', hotstring)
        else:
            new_line = line
        out_target_file.write(new_line)
    out_target_file.close()
    target_file.close()

    # Copy the script to run the design.
    shutil.copy(os.path.join(template_dir, 'run_motif_graft.slurm'), out_graft_dir)

    out_exelogs_err = os.path.join(out_graft_dir, 'exelogs/err')
    out_exelogs_out = os.path.join(out_graft_dir, 'exelogs/out')
    if not os.path.exists(out_exelogs_err):
        os.makedirs(out_exelogs_err)
        os.makedirs(out_exelogs_out)

def sele_exists(sele):
    return sele in cmd.get_names("selections")

#def sele_exists(sele):
#    sess = cmd.get_session()
##    for i in sess["names"]:
#        if type(i) is ListType:
#            if sele==i[0]:
#                return 1
#    return 0


def save_full(polar='polar', hphob='hphob', allaa='allaa',  binder_var=None, out_parent="/home/gainza/lpdi_fs/protein_designs/masif_seeds/by_scaffold/"):
    stored.next_obj.save(binder_var=binder_var, out_parent=out_parent)
    binder_name, name, output_dir = stored.next_obj.save(binder_var=binder_var, out_parent=out_parent)
    template_dir = "/home/gainza/lpdi_fs/protein_designs/masif_seeds/template_scripts/"

    # Copy the target to a new object
    cmd.copy('context', stored.next_obj.name )
    # Rename the chain
    cmd.alter('context', 'chain=\'A\'')

    # Create  binder
    cmd.copy('binder', '{}'.format(binder_name))
    cmd.alter('binder', 'chain=\'B\'')

    # Select the target and binder.
    cmd.select('complex_sel', 'binder or context')

    out_des_dir = os.path.join(output_dir, 'design')

    if not os.path.exists(out_des_dir):
        os.makedirs(out_des_dir)

    # Save it to new directory
    cmd.save(os.path.join(out_des_dir, name+'.pdb'), 'complex_sel')

    outstring = "NATAA\n"
    outstring += "USE_INPUT_SC\n"
    outstring += "start\n"

    # Create resfile from objects in polar, hphob, and allaa
    if sele_exists(allaa):
        stored.resis=[]
        cmd.iterate(allaa, 'stored.resis.append(resi)')
        resis = sorted(set(stored.resis), key=lambda x: int(x))
        for res in resis:
            outstring += "{} B PIKAA ADEFGHIKLMNQRSTVWY\n".format(res)

    if sele_exists(polar):
        stored.resis=[]
        cmd.iterate(polar, 'stored.resis.append(resi)')
        resis = sorted(set(stored.resis), key=lambda x: int(x))
        for res in resis:
            outstring += "{} B PIKAA DEHKNQRST\n".format(res)

    if sele_exists(hphob):
        stored.resis=[]
        cmd.iterate(hphob, 'stored.resis.append(resi)')
        resis = sorted(set(stored.resis), key=lambda x: int(x))
        for res in resis:
            outstring += "{} B PIKAA AFLIMVWYT\n".format(res)

    outfile = open(os.path.join(out_des_dir,'resfile'), 'w')
    outfile.write(outstring)
    outfile.close()

    # Copy design rosettascript and script
    shutil.copy(os.path.join(template_dir, 'design.sh'), out_des_dir)
    shutil.copy(os.path.join(template_dir, 'design.xml'), out_des_dir)
