This folder contains:

1) A script called 'task3_map_pdb_germlines.py'. This code assigns germline gene names to PDBs by comparing the output list from Task 1 to PyIgClassify db. It then matches on germline gene name to IMGT/3D dataset to select only those that interact with antigen.

2) A text file called 'SqLite.commands' containing the PyIgClassify database search terms.

3) A folder called 'fasta' with input files of germline gene sequences.

4) An input file called 'PyIgClassify_data.csv.'

5) An input file called 'IMGT_3D_res1to3Ang.csv.'

6) An output file called 'matches.csv.' These are matches of germline gene names to PDB code.

7) An output file called 'Final_list.csv.' This is the final list of 42 used for downstream analysis.
