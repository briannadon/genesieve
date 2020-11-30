import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from multiprocessing import Pool
from itertools import product
from sys import getsizeof

def read_quants(species):
    path = 'data/quants/'
    fullpath = path + species + "_quants.csv"
    quants = pd.read_csv("soy_quants.csv",index_col=0)
    
def get_genes(quant_table):
    return quant_table.index.tolist()

def find_all_coexps(gene,allgenes,quants):
    try:
        othergenes = allgenes
        othergenes.remove(gene)
    except ValueError:
        print(gene)
    g1 = np.array(quants.loc[gene])
    coexp_dict = {}
    for g in othergenes:
        g2 = np.array(quants.loc[g])
        coexp = spearmanr(g1,g2)
        coexp_dict[g] = coexp
    return coexp_dict

def get_all_coexpressions(target_gene,genelist,quants,threads=2):
    with Pool(processes=threads) as p:
        args = ((tg,genelist,quants) for tg in target_genes)
        gcoexpdict = p.starmap(find_all_coexps,args)
    return gcoexpdict
    