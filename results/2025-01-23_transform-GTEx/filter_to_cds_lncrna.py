# Filter MedianGeneCounts data to protein coding regions + lncRNAs
## input(s):
## - GeneCounts: HEAD/dat/MedianGeneCount_GTEx_v8.tsv
## - Annotation file: HEAD/dat/Homo_sapiens.GRCh38.113.chr.gtf

import sys
sys.path.insert(0, '../../src/') # make src available for loading scripts
from task_tracker import startTask, endTask, TaskTimes
from data_peek import peekTable
import numpy as np
import pandas as pd
import re

r
