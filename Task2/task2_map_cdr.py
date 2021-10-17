import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
#from IPython.display import display, HTML
import pandas as pd
import re
import os
from tabulate import tabulate




dict_cdr1 = {}
dict_cdr2 = {}
fasta_dir = './fasta/'
for fasta_name in os.listdir(fasta_dir):
    print(fasta_name)  
    for seq_record in SeqIO.parse(fasta_dir + fasta_name, 'fasta'):
        s = str(seq_record.seq).replace('.','')
        p = re.compile(r'^\w{18,24}?C\w{3}(\w{10,12}?)W[AIV]')
        i = 0 #set up counter to see if ever more than one match
        for m in p.finditer(s):
            i += 1
            if i >= 2:
                print (i, ' CDR1 matches found for ', seq_record.id)
            dict_cdr1[seq_record.id] = m # store the match object (if more than one match overwrites with last found)
            #print(seq_record.id,m.group(1)) 
            #print(m.group(1)) # 1 for CDR, 0 for whole regex match
            #print(m.span(1)) # span position

        if i >= 1:
            cdr1_match = dict_cdr1[seq_record.id]
            cdr1_end= cdr1_match.end(1) #retrieve the end position of CDR1
            #print(cdr1_end)
            start_pos = str(cdr1_end + 9) #9 after the end of CDR1 (to allow for five LEWIG AAs) CDR2 starts on 15th AA following
            regex_str = r'^\w{' + start_pos + '}LE\w{3}(\w{16,19}?)[KR][LIVFTA][TSIA]' #? mark is minimal matching
            p2 = re.compile(regex_str) 
            i = 0 #set up counter to see if ever more than one match
            for m2 in p2.finditer(s):
                i += 1
                #print(i)
                if i >= 2:
                    print (i, ' CDR2 matches found for ', seq_record.id)
                dict_cdr2[seq_record.id] = m2
                #print(seq_record.id,m2.group(1)) 
                #print(m2.group(1)) # 1 for CDR, 0 for whole regex match
                #print(m2.end(1) + 30) # span position + 30 to start CDR3
                #print(len(s)) ## not long enough for CDR3. CDR3 needs V,J recomb.




dict_comb={}  #Combine two dictionaries to one with important data only
for key in dict_cdr1:
    CDR1Start =dict_cdr1[key].start(1)
    CDR1End=dict_cdr1[key].end(1)
    CDR1Sequence=dict_cdr1[key].group(1)
    if key in dict_cdr2:
        CDR2Start =dict_cdr2[key].start(1)
        CDR2End=dict_cdr2[key].end(1)
        CDR2Sequence=dict_cdr2[key].group(1)
    else:
        CDR2Start =0
        CDR2End=0
        CDR2Sequence=''
    dict_comb[key] = [CDR1Start, CDR1End, CDR1Sequence, CDR2Start, CDR2End, CDR2Sequence]


cdrs_df =  pd.DataFrame.from_dict(dict_comb, orient = 'index', columns = ['CDR1Start', 'CDR1End', 'CDR1Sequence', 'CDR2Start', 'CDR2End', 'CDR2Sequence'])
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
#print(cdrs_df)
cdrs_df.to_csv('cdrs.csv', index = True)
    
    
    
             
   

    

