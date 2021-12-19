# ml-project-2-chopin-deux

This README describes the main scripts and data-files present in this project folder.

The relevance of each directory and file will be apparent when the two main goals of this project are kept in mind:
1. Comparing the classification accuracy of our protein surface-based machine learning algorithm to commonly used methods, mainly Rosetta.

2. Improving the accuracy of our algorithm by dataset augmentation. The key concept here is sampling of multiple patches per residue at a protein-protein interface.

## Scripts for main goal 1. Comparison
To perform a fair comparison, identical datasets are required.
To do this, one needs to know which residues were used for the single amino acid retrieval in the RosettaSurf publication. 

*Step 1:* In `data/rosetta_slurm`
 
