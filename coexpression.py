import pandas as pd
import numpy as np
from sys import getsizeof
from itertools import product
import mysql.connector

def get_sql_coexp(cursor,species,gene1,gene2,coexp_min=0.6):
    record_result = []
    table = f"{species}_coexp" # at/maize/rice/soy_coexp
    #YOU CANNOT PARAMETERIZE TABLE OR COLUMN NAMES! Be careful with string formatting below.
    query = f"SELECT coexp FROM {table} WHERE gene1 = %s AND gene2 = %s AND coexp > %s"
    if gene1 == gene2:
        return 1.0
    cursor.execute(query, (gene1,gene2,coexp_min))
    result = cursor.fetchall()
    if len(result) > 0:
        print(result)
        result = [float(x) for x in result[0]]
        result = sum(result)/len(result)
        return result
    #If gene1,gene2 didn't work, reverse it
    else:
        cursor.execute(query, (gene2,gene1,coexp_min))
        result = cursor.fetchall()
        if len(result) > 0:
            print(result)
            result = [float(x) for x in result[0]]
            result = sum(result)/len(result)
            return result
        else:
            return None


def get_all_coexps(gene_l_1, gene_l_2, conf_file,species,coexp_min=0.6):
    print(f"gene_l_1 is {len(gene_l_1)} long and gene_l_2 is {len(gene_l_2)} long")
    df_list = []
    db = mysql.connector.connect(option_files = conf_file, use_pure = True)
    cursor = db.cursor()
    all_combinations = list(product(gene_l_1,gene_l_2))
    for g1,g2 in all_combinations:
        coexp = get_sql_coexp(cursor,species,g1,g2,coexp_min)
        if coexp != None:
            df_list.append([g1,g2,coexp,'db gene','db gene','coexpression'])
    df = pd.DataFrame(df_list,columns=['gene1','gene2','coexpression','type1','type2','connection'])
    cursor.close()
    db.close()
    return df
