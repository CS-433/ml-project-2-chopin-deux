# ml-project-2-chopin-deux

This README describes the main scripts and data-files present in this project folder.

The relevance of each directory and file will be apparent when the two main goals of this project are kept in mind:
1. Comparing the classification accuracy of our protein surface-based machine learning algorithm to commonly used methods, mainly Rosetta.

2. Improving the accuracy of our algorithm by dataset augmentation. The key concept here is sampling of multiple patches per residue at a protein-protein interface.

## Scripts for main goal 1. Comparison
To perform a fair comparison, identical datasets are required.
To do this, one needs to know which residues were used for the single amino acid retrieval in the RosettaSurf publication. 

**Step 1:** In `data/rosetta_slurm` there are slurm files for every amino acid. `scripts/data_comparison_prep/get_rsurf_resi.py` returns a csv of all residues used by Andreas in RosettaSurf. It reads his slurm commands to get these. `data/rosetta_res.csv` is complete and contains 200 of each of the 20 AA's. Also columns for pdb_id, aa type and position.

**Step 2:** In `scripts/data_comparison_prep/get_specific_aa.py`, the test and training set is created to compare the performance of Chopin software to the RosettaSurf. The test set is what we will perform predictions on, the ~200 * 20 residues from 1495 pdb's. All other residues from those pdb's will serve as a training set. The pdbs are loaded in ascending order by pdb, chain, residue.

### Generating and retrieving surfaces using dMaSIF

**Step 2:** `/scripts/data_comparison_prep/find_pdbs.py` allows you to retrieve the relevant pdbs from the cluster after maSIF has generated the relevant surface data files. 

## Scripts for main goal 2. Data set augmentation

The `data_augmentation_prep` folder contains "jack_multi_patch_1.py" script 
##consider changing the name to something more descriptive
which generates multiple new patches of N surface points closest to the center of mass, starting from the fifth closest.

----
## Additional scripts for checks

- `/scripts/data_comparison_prep/get_pdb_set.py` A simple script to compare the pdb's used in the RosettaSurf comparison to the training and test set of dMaSIF.

- After generating surfaces with dMaSIF, we checked how many were succesfully retrieved, see `/scripts/data_comparison_prep/pdb_id_check.py`. It saves a list of pdb_ids for which no surfaces were generated. 

- The `/scripts/masif_plugin/` folder contains scripts to visualise protein surfaces generated using the MaSIF framework in PyMol. Most of  these scripts are from the original MaSIF publication (Gainza et al., 2020). The `/scripts/masif_plugin/dmasif_pymol.py` script is written by Arne Schneuing and allows for visualization of our surfaces. 


