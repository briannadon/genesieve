import pandas as pd
import numpy as np
import re
import itertools
from scipy.stats import norm

##################################################
## Helper script to make gene:qtl:normscore tables
##################################################
## The inputs should be of the form:
## trait: "trait_desc,chrom,start,end"
## and a gff of "chrom,start,end,gene"
##################################################

def find_qtl_genes(row,gff=gff):
    # takes in a row of the trait_df file, and uses the gff file
    # to find the genes that belong in that qtl. 
    # returns those as a table row.
    row = row[1]
    start = row['start']
    end = row['end']
    chrom = row['chromosome']
    gff_genes = gff.loc[(gff['start'] >= start) & 
                        (gff['end'] <= end) & 
                        (gff['chromosome']==chrom)]['gene'].tolist()
    return [(row['trait_name'],g) for g in gff_genes]

#the trait_df must be in the format discussed in the header of this file
def process_trait_df(trait_df):
    df = pd.DataFrame(list(itertools.chain(*map(find_qtl_genes,
                                                trait_df.iterrows()))),
                     columns=['trait','gene'])
    return df

def apply_norm_transform(l,min_score=0.5,max_score=1):
    ## This method will take in a LIST of genes/items,
    ## and return a list of norm scores of length len(list)
    ## to be applied manually to the final table.
    ## You will probably have to write another method to do that.
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

