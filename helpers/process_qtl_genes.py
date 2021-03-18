import pandas as pd
import numpy as np
import re
import itertools
from scipy.stats import norm

##################################################
## Helper script to make gene:qtl:normscore tables
##################################################
##################################################

def find_qtl_genes(row,gff=gff):
    row = row[1]
    start = row['start']
    end = row['end']
    chrom = row['chromosome']
    gff_genes = gff.loc[(gff['start'] >= start) & 
                        (gff['end'] <= end) & 
                        (gff['chromosome']==chrom)]['gene'].tolist()
    return [(row['trait_name'],g) for g in gff_genes]

def apply_norm_transform(l,min_score=0.5,max_score=1):
    n = len(l)+1
    ranks=range(1,n)
    cuts = [float(c)/float(n) for c in ranks]
    ea_scores = [norm.ppf(r,0,1) for r in cuts]
    if len(ea_scores) % 2 == 0:
        midpt = int(len(ea_scores)/2)
        reflect_pt = midpt+1
        half_ea_scores = ea_scores[:midpt]
    else:
        midpt = int(len(ea_scores)/2 + 1)
        reflect_pt = midpt - 1
        half_ea_scores = ea_scores[:midpt]
    ea_np = np.array(half_ea_scores)
    a = list(np.interp(ea_np, (ea_np.min(), ea_np.max()), (min_score, max_score)))
    b = a[:reflect_pt][::-1]
    return a+b