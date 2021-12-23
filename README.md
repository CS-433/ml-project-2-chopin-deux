# ml-project-2-chopin-deux

This README describes the main scripts and data-files present in this project folder.

The relevance of each directory and file will be apparent when the main goal of this project are kept in mind:
*Comparing the classification accuracy of our protein surface-based machine learning algorithm to commonly used methods, mainly Rosetta.*

### Compare the model with Rosetta

**Step 1:** Clone the git repository to your local device.

**Step 2:** Open `/scripts/Test_set_and_comparison.ipynb` and run the notebook.

**Step 3:** Find the comparison at the end of the notebook.

## Used scripts for main goal Comparison
To perform a fair comparison, identical datasets are required.
To do this, one needs to know which residues were used for the single amino acid retrieval in the RosettaSurf publication. 

- In `data/rosetta_slurm` there are slurm files for every amino acid. `scripts/data_comparison_prep/get_rsurf_resi.py` returns a csv of all residues used by Andreas in RosettaSurf. It reads his slurm commands to get these. `data/rosetta_res.csv` is complete and contains 200 of each of the 20 AA's. Also columns for pdb_id, aa type and position.

- In `scripts/data_comparison_prep/get_specific_aa.py`, the test and training set is created to compare the performance of Chopin software to the RosettaSurf. The test set is what we will perform predictions on, the ~200 * 20 residues from 1495 pdb's. All other residues from those pdb's will serve as a training set. The pdbs are loaded in ascending order by pdb, chain, and then residue.

- The `/scripts/data_comparison_prep/zander.py` adds new features to the data set, including geometric angles between the residues.

- `/scripts/RosettaSurf_benchmark.ipynb` is the script from the original RosettaSurf github which contains the data for the Rosetta comparison.

*Generating and retrieving surfaces using dMaSIF*

- `/scripts/data_comparison_prep/find_pdbs.py` allows you to retrieve the relevant pdbs from the cluster after maSIF has generated the relevant surface data files. 

----
## Additional scripts for checks

- `/scripts/data_comparison_prep/get_pdb_set.py` A simple script to compare the pdb's used in the RosettaSurf comparison to the training and test set of dMaSIF.

- After generating surfaces with dMaSIF, we checked how many were succesfully retrieved, see `/scripts/data_comparison_prep/pdb_id_check.py`. It saves a list of pdb_ids for which no surfaces were generated. 

- The `/scripts/masif_plugin/` folder contains scripts to visualise protein surfaces generated using the MaSIF framework in PyMol. Most of  these scripts are from the original MaSIF publication (Gainza et al., 2020). The `/scripts/masif_plugin/dmasif_pymol.py` script is written by Arne Schneuing and allows for visualization of our surfaces. 

## Helpers and other auxiliary files

- The `scripts/data_comparison_prep/get_aa_helper.py` contains helper functions that are then used in `scripts/data_comparison_prep/get_specific_aa.py` to obtain training and testing dataset for comparison with RosettaSurf.
- The `scripts/data_comparison_prep/get_aa_loader.py` is used to load, sort, and pre-clean the dataframe containing Pdbs used in Rosetta and the data points containing surfaces' features by removing the unnecessary columns.


------
### Procedure to train the model

**Step 1:** Get the training data from the google drive (https://drive.google.com/drive/folders/16iuOe2GBY9-flGjP7FfbjYjIabSq2Uy1?usp=sharing)
Place all the directories in the "data" directory. 

**Step 2:** Open `/scripts/train.ipynb` on a cluster and use the tensorboard to monitor the progress.
For instructions on how to open a notebook and the cluster and a tensorboard on a cluster see the instructions below:

## Tensorboard
```
module load gcc
module load python
source amino_env/python_env/bin/activate
``` 
Get hostname (it is a number) with:
```
hostname -I | awk '{print $1}'
```
Fill it in here:
```
tensorboard --logdir=runs --host hostname --port $SERVER_PORT
```
On you local device:
```ssh -L 5555:hostname:8888 username@izar.epfl.ch```


## Jupyter notebook
```
jupyter notebook --no-browser --port=1234 --ip=$(hostname -i)
ssh -NL 1234:hostname:1234 user@izar.epfl.ch
```
Copy the url that will appear in remote terminal and paste it in your browser.
