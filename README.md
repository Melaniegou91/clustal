# clustal
Short project whose aim is to reprogram the clustal software

This script allows you to display a phylogenetic tree in order to make a multiple alignment.
Unfortunately, it does not allow multiple alignment, it only displays the tree
## How to download my conda environment


To download the conda environment from the environment.yml file you can use the following command:


```{}
conda env create -f environment.yml
```
## First step 

To start, you need to create a fasta file containing the different sequences you want to compare.
These sequences must be protein sequences.
To retrieve the sequences to compare, you can go to the uniprot website
URL : https://www.uniprot.org/

## second step
The script uses the following parameters : 

- --file_name The name of your file containing the sequences
- --seq1 Indicate the first sequence you want to compare for alignment: 0 for the first sequence, 1 for the second sequence...
- --seq2 Indicate the second sequence you want to compare for alignment: 0 for the first sequence, 1 for the second sequence...
- --file_matrix_blosum Indicate the path leading to your file containing the blosum62 sequence

### Example

The sequences.fasta file is present in the github directory
You can run the script with the following command : 
```{R}
 python src/GOU_clustal_project.py --file_name data/sequences.fasta --seq1 0 --seq2 1 --file_matrix_blosum data/blosum62.txt
```

## Results

The script will open a window containing the phylogenetic tree.
If you execute the command line shown above as an example you should obtain the following tree

![front-page](img/phylogenetic_tree.png)


The sequence ids represent your sequences in the order they are in the file.
You must add +1 to the id to find your sequence.
That is to say that id 0 corresponds to sequence 1 in the file, id 1 to the second sequence and so on.

The script returns the following alignment :

![front-page](img/alignement.png)
