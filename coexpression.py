import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from multiprocessing import Pool
from itertools import product
from sys import getsizeof

#unused
def find_all_coexps(gene,allgenes,quants):
    try:
        othergenes = allgenes
        othergenes.remove(gene)
    except ValueError:
        print(gene + ' is not in the list of genes')
    g1 = np.array(quants.loc[gene])
    coexp_dict = {}
    for g in othergenes:
        g2 = np.array(quants.loc[g])
        coexp = spearmanr(g1,g2)
        coexp_dict[g] = coexp
    return coexp_dict

def get_sql_coexp(conf_file,species,gene1,gene2):
    record_result = []
    db = mysql.connector.connect(option_files = conf_file, use_pure = True)
    cursor = db.cursor()
    table = f"{species}_coexp" # at/maize/rice/soy_coexp
    query = f"SELECT coexp FROM {table} WHERE gene1 = {gene1} AND gene2 = {gene2}"
    if gene1 == gene2:
        return 1.0
    cursor.execute(query, (gene1, gene2))
    result = cursor.fetchall()
    if len(result) > 0:
        result = [float(x) for x in result]
        result = sum(result)/len(result)
        return result
    #If gene1,gene2 didn't work, reverse it
    else:
        cursor.execute(query, (gene2,gene1))
        result = cursor.fetchall()
        if len(result) > 0:
            result = [float(x) for x in result]
            result = sum(result)/len(result)
            return result
        else:
            print("\t\tboth queries failed to return results")
    cursor.close()
    db.close()

def 
