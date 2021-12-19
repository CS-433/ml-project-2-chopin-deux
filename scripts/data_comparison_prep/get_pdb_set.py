"""
Created on 24.11.21
A simple script to compare the pdb's used in the RosettaSurf comparison to the training and test set of dMaSIF
@author: maxjansen
"""
from pathlib import Path
import re
import glob
from itertools import chain

data_dir= Path('../../data/')
pdb_list = []
for filepath in data_dir.glob("*.slurm"):
    f = open(filepath, "r")
    for line in f.readlines():
        if line.startswith("pdb_id"):
            line = line.split("pdb_id=(")[1]
            line = line.split(")")[0]
            line = re.sub('"', '', line)
            pdb_list.append(list(line.split(" ")))



unnested_pdb = list(chain.from_iterable(pdb_list))
unique_pdb = (list(set(list(unnested_pdb))))


############
test_list = []
test_f = open(data_dir / "testing_ppi.txt", "r")
for line in test_f.readlines():
    line = line.split("\n")[0]
    test_list.append(line)


#######
unique_pdb_set = set(unique_pdb)
dmasif_test_set = set(test_list)



training_list = []
training_f = open(data_dir / "training_ppi.txt", "r")
for line in training_f.readlines():
    line = line.split("\n")[0]
    training_list.append(line)
dmasif_train_set = set(training_list)

# All pdb's used by RosettaSurf in one list
rosetta_surf_list = list(sorted(unique_pdb_set))

# Save pdb/names of protein complexes used by RosettaSurf that intersect with dMaSIF training set (use for surfaces!)
textfile = open(data_dir / "rosetta_surf_select.txt", "w")
for element in rosetta_surf_list:
    textfile.write(element + "\n")
textfile.close()

# Check how long this is
#print(len(list(sorted(dmasif_train_set.intersection(unique_pdb_set)))))

first_chains = []
second_chains = []

# for i in rosetta_train_list:
#     pdb_id = i.split('_')[0]
#     first_chain = i.split('_')[1]
#     first_chains.append(pdb_id + '_' + first_chain + '.pdb')
#     second_chain = i.split('_')[2]
#     second_chains.append(pdb_id + '_' + second_chain + '.pdb')
#
# all_pdbs = first_chains + second_chains
# Save all pdb_ids used by RosettaSurf to get pdb's from dMaSIF preprocessing.
# textfile = open(data_dir / "pdb_chain_select.txt", "w")
# for element in all_pdbs:
#     textfile.write(element + "\n")
# textfile.close()