import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
import csv


### Match allele name in genes within fasta directory with allele names and associated PDBs in PyIgClassify

fasta_dir = './fasta/'

pyg_df = pd.read_csv('PyIgClassify_data.csv', index_col = 0, usecols=['PDB','frame_ig']) ## specify PDB as index
#print(pyg_df)#642

list_alleles = []
for gene_name in os.listdir(fasta_dir):
    for seq_record in SeqIO.parse(fasta_dir + gene_name, 'fasta'):
        list_alleles.append(seq_record.id)
#print(len(list_alleles))   ## 147 alleles in 25 gene_groups
matches_df = pyg_df[pyg_df['frame_ig'].isin(list_alleles)]
matches_df = matches_df.reset_index().drop_duplicates().set_index(['PDB']) ## remove duplicate PDB
print (matches_df) #211
matches_df.to_csv('matches.csv', index = True)



      
## match germline genes in PygClassify to IMGT/3D, i.e those that interact with antigen

data_dir = '../Data_downloads/'
csv_name = 'IMGT_3D_Res1to3Ang.csv'
#imgt3D_list = []

matches = pd.read_csv('matches.csv')
imgt3D_df = pd.read_csv(f'{data_dir}{csv_name}')
imgt3D_df ['IMGT entry ID'] = imgt3D_df['IMGT entry ID'].str.upper()
#print(imgt3D_df) #392

#comb_df = pd.merge(matches, imgt3D_df, left_on = 'PDB', right_on = 'IMGT entry ID')
comb_df = matches.merge(imgt3D_df, left_on = 'PDB', right_on = 'IMGT entry ID', left_index=False, right_index=False)
comb_df = comb_df[['PDB','frame_ig', 'IMGT receptor description','Ligand(s)','Resolution']]
print(comb_df) #42
#comb_df.to_csv('Final_list.csv', index = True)
#matches.info() ##check datatype 
#imgt3D_df.info()
 
