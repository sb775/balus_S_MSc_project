import pandas as pd

import Bio
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser



import pathlib
from pathlib import Path # a good way to move out of Windows

main_folder = Path('C:/Users/b/Desktop/Suchu/MSc_projects/')  #to move out of Windows


### Create mutation lists to feed into mCSM-PPI2 webserver. Loop through each allele. If aa is not a match to germline, add to list.

#Open up py1g_classify.  PyIg_classify needed to find the original_chain for each pdb
classify_fname = 'PyIgClassify_and_final_list_auto.csv'
classify_file = main_folder / 'Code' / classify_fname
#print(classify_file)
classify_df = pd.read_csv(classify_file)

#Open up final_list to get PDB and frame_ig codes
list_fname = 'Final_list.csv'
list_file = main_folder / 'Code' / list_fname
list_df = pd.read_csv(list_file)
#print(list_df)
for index, row in list_df.iterrows():
    frame_ig = str(row['frame_ig']).upper() ##Force upper
    pdb = str(row['PDB']).upper()
    print(frame_ig + ' ' + pdb)

    if True:  #pdb in ['6PCU']:

        #Get rows from classify_df where the pdb code matches
        select_df = classify_df[classify_df['PDB'] == pdb]
        #Check for multiple matches
        match_count = len(select_df.index) #if more than one will need to check for the right one later
        if match_count > 1:
            print ( 'ERROR - multiple records')
        original_chain = select_df.iloc[0]['Original_Chain'] #just choose the first one for now
        #print(original_chain)

        #Open the file with the germline sequence in it.  It should be in fasta format.
        frame_ig_new = frame_ig.replace('-', '_')
        frame_ig_new = frame_ig_new.replace('*', '_')
        frame_ig_new = frame_ig_new.lower()
        germline_fname = frame_ig_new + '_pdb_seq.txt'


        germline_file = main_folder / 'Data_downloads' / 'pdbseq&germline'/ germline_fname
        #print (germline_file)
        r_dict = SeqIO.to_dict(SeqIO.parse(germline_file, 'fasta'))
        seq_germline = r_dict[frame_ig].seq
         
        #Open the file with the pdb sequence in it (i.e sequence in pdb format).
        pdb_fname = pdb.lower() + '.pdb'
        pdb_file = main_folder / 'Data_downloads' / 'pdbs'/ pdb_fname
        
        #Read in three-letter to one-letter file mapping file
        key_dict = {}
        with open('ThreeLetterAACode.csv', 'r') as f:
            for line in f:
                (key, val) = line.rstrip().split(',')
                key_dict[key] = val

        #Now read from pdb file 
        pdb_d = {}  ##{position,1-letter code)
        high_pos = 0 #to store the highest position read in the sequence from the residues in the Atom records
        low_pos = 99999999 #to store the lowest position read in the sequence from the residues in the Atom records
        parser = PDBParser(PERMISSIVE = True, QUIET = True)
        ignore_ls = ['MPD', 'HOH', 'EDO', ' CL', 'GOL', 'LDA', 'SO4']  #ignore these when returned.  Not amino acids.

        s = parser.get_structure(pdb, pdb_file)
        model = s[0]
        chain = model[original_chain]
        count=0
        
        for res in chain:
            count+=1
            if count <=225:  #germline sequences are about 100 long.  Go to 225 here because 6PCU start of germline sequence is at position 119 in the PDB file.
                resname = res.get_resname()
                if resname not in ignore_ls:
                    position = str(res.id[1]) + res.id[2].strip()
                    pdb_d[position] = key_dict[res.get_resname()]
                    high_pos = max(high_pos, res.id[1]) ## Find highest position
                    low_pos = min(low_pos, res.id[1]) #Find lowest position

        #add in the missing residues (identified by REMARK 465 in the pdb file)
        for missing_res in s.header['missing_residues']:
            if missing_res['chain'] == original_chain:
                pdb_d[str(missing_res['ssseq'])] =  key_dict[missing_res['res_name']]

        position_d = {} #stores the position of the residue + insertion code against the integer position in the sequence
        pdb_seq = ''
        count = 0
        if pdb in ['2CMR', '2XRA', '5BV7', '3U2S']:  #for these four the germline starts with records that are missing from the PDB file
            low_pos = 1
        elif pdb == '6PCU':
            low_pos = 119  #Sequence in PDB file starts later - position 119 corresponds to the start of the germline sequence.    
        for i in range(low_pos,high_pos+1):
            num  = str(i) ## pdb_dict key is a string (a lot of positions have a,b,cs)
            if num in pdb_d:  #construct the sequence from the pdb using Xs for deletions (to help line up with the germline sequence)
                pdb_seq += pdb_d[num]  
            else:
                pdb_seq += 'X'  #for deletions e,g 3NPS position 40 (no entries in pdb file).
            count += 1
            position_d[count] = num
            for letter in 'ABCDEFGH':
                alt_num = num + letter
                if alt_num in pdb_d:
                    pdb_seq += pdb_d[alt_num]
                    count += 1
                    position_d[count] = alt_num ## Adds actual position to a dictionary (i.e of a,b etc)
        
                
        print(seq_germline)
        print(pdb_seq[:102])
                       
            
        #Now compare the mutant and wild-type strings
       
        i = 0
        result_list = []
        for residue in seq_germline:
            pdb_s = pdb_seq[i]
            if residue != pdb_s:
                #result_list.append(original_chain.upper() + ' ' + pdb_s + str(i + 1) + residue)
                result_list.append(original_chain.upper() + ' ' + pdb_s + position_d[i + 1] + residue)
            i +=1

        #Now save in text file
        muts_fname = pdb.lower() + '_mut_list.txt'
        muts_file = main_folder / 'mscm_ppi2_mutlists_auto' / muts_fname

        if False:  #Decide whether to write new results or check them against existing stored results 
            with open(muts_file, 'w') as text_file:
                for res in result_list:
                    #write results to file
                    text_file.write(res + '\n')
        else:
            with open(muts_file, 'r') as text_file:
                for res in result_list:
                    #check results against existing file
                    if res != text_file.readline().strip():
                        print(res)
                    

                
                


