import pandas as pd
import numpy as np
from sys import getsizeof
from itertools import permutations

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
    print(f"the result of {gene1},{gene2} was coexp={result}")


def get_all_coexps(gene_l_1, gene_l_2, conf_file,species):
    df_list = []
    all_combinations = [list(zip(each_permutation, gene_l_2)) 
            for each_permutation in permutations(gene_l_1, len(gene_l_2))]
    print(all_combinations)
    for g1,g2 in all_combinations:
        coexp = get_sql_coexp(conf_file,species,g1,g2)
        df_list.append([g1,g2,coexp,'db gene','db gene','coexpression'])
    df = pd.DataFrame(df_list,columns=['gene1','gene2','coexpression','type1','type2','connection'])
    return df
