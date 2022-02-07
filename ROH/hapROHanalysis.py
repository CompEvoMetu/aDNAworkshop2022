# hapROHanalysis.py:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os as os
import sys as sys
import multiprocessing as mp


from hapsburg.PackagesSupport.hapsburg_run import hapsb_ind
iids = ['Ash128','cay007','cay013','GD13a','Bar8','Bar31','CCH294','CCH311','cth739','cth747']
for iid in iids:
    hapsb_ind(iid=iid, chs=range(1,22), processes=6,
              path_targets='/mnt/NEOGENE1/toTransfer/workshopfiles/roh/cay_hapROH',
              h5_path1000g='/mnt/NEOGENE1/toTransfer/workshopfiles/roh/chr',
              meta_path_ref='/mnt/NEOGENE1/toTransfer/workshopfiles/roh/meta_df_all.csv',
              folder_out='/mnt/NEOGENE1/home/user', prefix_out='',
              e_model="haploid", p_model="Eigenstrat", n_ref=2504,
              random_allele=True, readcounts=False,
              delete=False, logfile=True, combine=True)
# post-processing
from hapsburg.PackagesSupport.pp_individual_roh_csvs import pp_individual_roh
# hapROHanalysis.py file end
