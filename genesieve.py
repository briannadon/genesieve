import subprocess
import time
from sys import argv
import pandas as pd
import os
import pathlib
import re

import annotate
import blast
import sanitize
import scoring
import phenotype
import coexpression

if __name__=="__main__":
    script, in_fasta, in_pheno = argv
    timestamp = str(time.time())
    timestamp = f"genesieve_{timestamp}"
    
    cwd = os.getcwd()
    
    in_fasta=f"{cwd}/{in_fasta}"
    in_pheno=f"{cwd}/{in_pheno}"
    

    augustus_protein_script = f"{cwd}/genesieve/helpers/get-fasta.sh"
    coexp_min = 0.6
    pheno_sim_min = 0.7
        
    db_conf = "/home/bnadon/mysql.test/db.conf" 
    gene_db = f'{cwd}/testdata/Oryza_sativa.faa'

    qtl_db = f"{cwd}/testdata/rice_qtl_genes.csv"
    coexp_db = f"{cwd}/testdata/ricecoexp_sample_1000_LOC_Os05g04990.1.csv"
    #blast_file = 'rice_height_test/scripts/Ph5-1_vs_rice.blast'
    pheno_model = f'{cwd}/testdata/model/dbow_wiki_2.model'
    
    pathlib.Path(f"{cwd}/{timestamp}").mkdir(parents=True, exist_ok=True)

    os.chdir(f'{cwd}/{timestamp}')

    #run annotation first
    augout = timestamp+".augustus"
    aquery = annotate.augustus_query(in_fasta,"rice",augout)
    a = subprocess.Popen(aquery)
    
    #While augustus is running, load up the phenotype stuff
   
    # Read in every single QTL, then get a list of all the traits
    qtl_table = pd.read_csv(qtl_db)
    trait_list = list(set(list(qtl_table.trait)))

    # get the distance of your input pheno to every pheno in the database 
    pheno_table = phenotype.get_pheno_results(in_pheno,pheno_model,
                                              trait_list,pheno_sim_min)


    #Here, we should filter the phenotype table. It's already filtered,
    # but I want to pare down results for now.
    pheno_table_filtered = pheno_table.nlargest(int(.1 * len(pheno_table)),'similarity') 
    
    #wait till annotation is done
    a.wait()


    aa_file = annotate.get_proteins(timestamp,get_fasta_script=augustus_protein_script)
 
    #Get the BLAST running in the background as "b"
    blast_output = timestamp + ".blast"
    bquery = blast.blast_query(aa_file,gene_db,blast_output)
    b = subprocess.Popen(bquery)
    
    #load in the coexp data
    #We will NOT be filtering it right now
    #TK: filter the coexp values
    #coexp_table = pd.read_csv(coexp_db)
    ##TK: Delete this line when we have properly reformatted the data
    #coexp_table['gene1'] = [re.sub('\.\d','',x) for x in list(coexp_table['gene1'])]
    #coexp_table['gene2'] = [re.sub('\.\d','',x) for x in list(coexp_table['gene2'])]
    
    #if the BLAST is still running, wait till it finishes
    b.wait()

    blast_table = blast.process_blast(blast_output)
    blast_table = blast_table.assign(type1='input gene',type2='db gene')
    ##TK: delete this line when we have properly formatted the data
    blast_table['subject'] = [re.sub('\.\d','',x) for x in list(blast_table['subject'])]
    blast_hits = blast_table.subject.tolist()
    
    
    #########################
    #Create the Results Table
    #########################
    
    #select the QTL that are associated with the phenotypes
    good_pheno_matches = pheno_table_filtered['db pheno'].tolist()
    selected_qtl = qtl_table.loc[
            qtl_table['trait'].isin(good_pheno_matches)]
    selected_qtl = selected_qtl.assign(
            type1='db pheno',type2='db gene',connection='qtl')

    print(f"the length of selected_qtl is {len(selected_qtl)}.")
    
    #select the coexpressed genes with the BLAST hits
    ##TK: refactor when full database is available to filter for blast
    ##hits being in either column gene1 or gene2

    #Getting coexpression from the sql db (db.conf set above)
    #tk: implement all species, change "rice" to 'species'
    
    qtl_connected_genes = selected_qtl['gene'].tolist()
    selected_coexp = coexpression.get_all_coexps(blast_hits,qtl_connected_genes,db_conf,'rice',coexp_min=coexp_min) 
    print(f"selected_coexp has {len(selected_coexp)} rows")

    #selected_coexp = coexp_table.loc[
    #    coexp_table['gene2'].isin(blast_hits)]
    #selected_coexp = selected_coexp.assign(type1='db gene',type2='db gene',connection='coexpression')

    #the second depth of qtl: coexpressed genes with QTL connections
    #check if any of the selected coexpressed genes have phenotype connections

    #Fix this - isin doesn't work this way.
    if not selected_coexp.empty and not selected_qtl.empty:
        selected_qtl_genes = set(selected_qtl['gene'].tolist())
        coexp_g1s = set(selected_coexp['gene1'].tolist())
        coexp_g2s = set(selected_coexp['gene2'].tolist())
        full_coexps = coexp_g1s | coexp_g2s
        coexp_qtl_hits = selected_qtl.loc[selected_qtl.isin(full_coexps)['gene']]



    #Cleaning up headers, aligning them to be renamed
    header = ['item1','item2','weight','source1','source2','connection']
    pheno_table = pheno_table_filtered[['input pheno','db pheno','similarity','type1','type2','connection']]
    selected_qtl = selected_qtl[['trait','gene','norm_score','type1','type2','connection']]
    blast_table = blast_table[['query','subject','pid','type1','type2','connection']]
    try:
        s_coexp_empty = False
        selected_coexp = selected_coexp[['gene1','gene2','coexpression','type1','type2','connection']]
        coexp_qtl_hits = coexp_qtl_hits[['trait','gene','norm_score','type1','type2','connection']]
    except (TypeError,NameError):
        s_coexp_empty = True
    
    
    #rename the columns
    pheno_table.columns=selected_qtl.columns=blast_table.columns=header
    if not s_coexp_empty:
        selected_coexp.columns=coexp_qtl_hits.columns=header
        tables = [pheno_table,blast_table,selected_qtl,selected_coexp,coexp_qtl_hits]
    else:
        tables = [pheno_table,blast_table,selected_qtl]


    #The results table
    if not s_coexp_empty:
        results_table = pd.concat([pheno_table,blast_table,selected_qtl,selected_coexp,coexp_qtl_hits],
                             ignore_index=True)
    else:
        results_table = pd.concat([pheno_table,blast_table,selected_qtl],
                                             ignore_index=True)
    results_table = results_table.drop(results_table[(results_table['connection'] == 'text') & 
                                                     (results_table['weight'] < pheno_sim_min)].index)
    num = results_table._get_numeric_data()
    num[num < 0] = 0
    
    results_table.to_csv(f"{timestamp}_resultstable.csv",index=False)
