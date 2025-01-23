# Get median read counts of each gene by tissue
## input(s):
## - Metadata: HEAD/dat/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
## - GeneCounts: HEAD/dat/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct

import sys
sys.path.insert(0, '../../src/') # make src available for loading scripts
from task_tracker import startTask, endTask, TaskTimes
from data_peek import peekTable
import numpy as np
import pandas as pd
import re


tsk = "load data"

startTask(tsk)
Metadata = pd.read_csv('../../dat/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t')
Metadata = Metadata[['SAMPID','SMTSD']] # Extract sample id and tissue id
GeneCounts = pd.read_csv('../../dat/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct', skiprows=2, sep='\t') #skiprows to avoid metdata
GeneCounts['Name'] = GeneCounts['Name'].str.split('.').str[0]
peekTable(GeneCounts)
endTask(tsk)

# Generate tissue sample membership {'Adipose-Subcut': (Sample1,2,5), 'Heart': (3,4), ...}

tsk = "generate sample-tissue membership"

startTask(tsk)
tissue_sampling_mapping = {} 
unique_tissues = tuple(Metadata.SMTSD.unique())
for tis in unique_tissues:
    tissue_sampling_mapping[tis] = Metadata[Metadata['SMTSD'] == tis]['SAMPID']
endTask(tsk)

# Calculate median gene count for every tissue

tsk = "calculate median gene count per tissue"

startTask(tsk)
MedianCount_df = GeneCounts[['Name','Description']]
for tis in unique_tissues:
    tis_samples = tissue_sampling_mapping[tis].to_list()
    tis_median_gene_count = GeneCounts[GeneCounts.columns.intersection(tis_samples)].median(axis=1).round()
    MedianCount_df[tis] = tis_median_gene_count
endTask(tsk)

peekTable(MedianCount_df)
tsk = "save median gene count gtex"

startTask(tsk)
MedianCount_df.to_csv('../../dat/MedianGeneCount_GTEx_v8.tsv', sep='\t', index=False)
endTask(tsk)

TaskTimes()
