import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
#import numpy as np
import pandas as pd
from Bio.Align import AlignInfo
from IPython.display import display, HTML




### Option   - Straight into dictionaries

### A.S: 20/5/21 - leave IMGT for now. May need to come back if not enough data in VDJ.
### I have removed 9 entries that were problematic to tanslate (see problematictranslation_imgt.fasta)


## 2.1 - imgt dataset into dictionaries. Data file name 'imgtref_200521.fasta' (gapped data).

## Prints 499 lines. Dirty data. Stop codons in the middle.

"""
add_prefix = '_imgt'

dict_imgtref = {}
for seq_ref in SeqIO.parse('imgtref_200521.fasta', 'fasta'):
    temp_imgt_id = seq_ref.id.split("|")[1]
    #print(temp_imgt_id)
    tmp_imgt_shortid = temp_imgt_id.split("*")[0] + add_prefix ##split before *
    #print(tmp_imgt_shortid)
    imgt_upper = seq_ref.seq.upper()
    #print(imgt_upper)
    triplet = imgt_upper[:len(imgt_upper)-len(imgt_upper) % 3]
    #print(triplet)
    triplet_temp = triplet.translate(gap='.')
    #print(triplet_temp)
    dict_imgtref[temp_imgt_id] = (triplet_temp,tmp_imgt_shortid) ## tuple value field
print(dict_imgtref)


## Put dict_imgtref into a dataframe

#imgt_df = pd.DataFrame.from_dict(dict_imgtref, orient = 'index')
#print(imgt_df.head(4))
"""


## 2.2 - vdjbase dataset into dictionaries. Data input file name 'VDJbase_dataset130521.fasta' (gapped data).

## Prints 290 lines 

dict_vdjbase = {}

for seq_record in SeqIO.parse('VDJbase_dataset130521.fasta', 'fasta'):
    temp_record_id = seq_record.id.split("|")[0]
    #print(temp_record_id)
    temp_short_id = temp_record_id.split("*")[0]
    #print (temp_short_id)
    seq_upper = seq_record.seq.upper()## Change case
    #print(seq_upper)
    aa_triplet = seq_upper[:len(seq_upper)-len(seq_upper) % 3]
    #print(aa_triplet)
    aa_temp = aa_triplet.translate(gap='.') 
    #print(len(aa_temp)) ## column size for freq table
    dict_vdjbase[temp_record_id]= (aa_temp,temp_short_id) 
#print(dict_vdjbase)

vdjbase_df = pd.DataFrame.from_dict(dict_vdjbase, orient = 'index', columns=['Sequence','GeneName'])
#print(vdjbase_df)
#print(vdjbase_df.describe())

gene_groups = vdjbase_df.groupby('GeneName').count()
pd.set_option('display.max_rows', None, 'display.max_columns', None)
print(gene_groups)
mult_alleles = gene_groups[gene_groups['Sequence']>3] ## Select genes with x number of alleles
print(mult_alleles)


count = 0
dict_alignstats = {}
for row in mult_alleles.iterrows():
    gene = row[0]
    #print(gene)
    f_out = './fasta/' + gene + '_out.fasta'
    seq_selected = []
    for key in dict_vdjbase:
        if dict_vdjbase[key][1] == gene:
            #print(dict_vdjbase[key][0]) ## aa sequences for all genes in mult_alleles
            record = SeqRecord(dict_vdjbase[key][0],id = key, description = '')
            #print(record)
            seq_selected.append(record)

    #SeqIO.write(seq_selected, f_out, 'fasta')
    
    maxlen = max(len(rcd.seq) for rcd in seq_selected)
    for rcd in seq_selected:
        if len(rcd.seq) != maxlen:
            padded_seq = str(record.seq).ljust(maxlen, '.')
            rcd.seq = Seq(padded_seq)
    #print(count)
    #print(gene)
    align = MultipleSeqAlignment(seq_selected)
    summary_align = AlignInfo.SummaryInfo(align)
    consensus = summary_align.dumb_consensus(ambiguous='.')
    ls = []
    ls.append(len(seq_selected))
    allele_count = 0
    ls_diffs = []
    ignoreambiguous = True #if true will not count the dots in alleles as differences
    for allele in align:
        char_count = 0
        diff_count = 0
        
        #print(consensus)
        #print(allele)
        #print(align)
        for aa in consensus:
            if aa != allele[char_count]:
                if ignoreambiguous and allele[char_count] != '.':
                    diff_count +=1
            char_count += 1
        ls_diffs.append(diff_count)
        #print(diff_count)
        allele_count +=1
    change_count = len([x for x in ls_diffs if x > 0])
    if change_count > 0:
        polymorph_ave = sum(ls_diffs) / change_count
    polymorph_max = max(ls_diffs)
    dict_alignstats[gene] = [allele_count, change_count,polymorph_ave,polymorph_max]
    observed_frequencies = align.substitutions
    observed_frequencies = observed_frequencies.select("AEIKQPRTVW")
    #print(gene)
    #print(observed_frequencies)
#print(dict_alignstats)
summary_df = pd.DataFrame.from_dict(dict_alignstats, orient = 'index', columns=['NumAlleles','NumWithDiffsFromConsensus', 'AveNumPolymorhisms', 'MaxNumPolymorphisms'])
#print(summary_df)
#summary_df.to_csv('summary_table.csv', index = True)





    



           


