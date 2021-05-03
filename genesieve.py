import subprocess
import time
from sys import argv
import annotate
import blast
import sanitize
import scoring
import phenotype
import coexpression
import pandas as pd
import os
import pathlib

if __name__=="__main__":
    script, in_fasta, in_pheno = argv
    timestamp = str(time.time())
    timestamp = f"genesieve_{timestamp}"
    
    cwd = os.getcwd()
    
    in_fasta=f"{cwd}/{in_fasta}"
    in_pheno=f"{cwd}/{in_pheno}"

    augustus_protein_script = f"{cwd}/genesieve/helpers/get-fasta.sh"
    coexp_min = 0.65
    pheno_sim_min = 0.7
        
    
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
    
    qtl_table = pd.read_csv(qtl_db)
    trait_list = list(set(list(qtl_table.trait)))

    #sanitize pheno input:
    
    pheno_table = phenotype.get_pheno_results(in_pheno,pheno_model,
                                              trait_list,pheno_sim_min)
    
    
    #check if annotation is done
    while not a.poll():
        time.sleep(5)
    
    aa_file = annotate.get_proteins(timestamp,get_fasta_script=augustus_protein_script)
 
    #Get the BLAST running in the background as "b"
    blast_output = timestamp + ".blast"
    subprocess.run(blast.blast_query(aa_file,
                                      'data/testdata/Orzya_sativa.faa',
                                      blast_output
                                     ))
    
    #load in the coexp data
    #We will NOT be filtering it right now
    #TK: filter the coexp values
    coexp_table = pd.read_csv(coexp_db)
    
    ##TK: Delete this line when we have properly reformatted the data
    coexp_table['gene1'] = [re.sub('\.\d','',x) for x in list(coexp_table['gene1'])]
    coexp_table['gene2'] = [re.sub('\.\d','',x) for x in list(coexp_table['gene2'])]
    
    #if the BLAST is still running, wait till it finishes
    while not b.poll():
        time.sleep(5)

    blast_table = blast.process_blast(blast_output)
    blast_table = blast_table.assign(type1='input gene',type2='db gene')
    
    ##TK: delete this line when we have properly formatted the data
    blast_table['subject'] = [re.sub('\.\d','',x) for x in list(blast_table['subject'])]
    
    blast_hits = blast_table.subject.tolist()
    
    
    #########################
    #Create the Results Table
    #########################
    
    #select the QTL that are associated with the BLAST hits
    selected_qtl = qtl_table.loc[
        qtl_table['gene'].isin(blast_hits)]
    selected_qtl = selected_qtl.assign(type1='db pheno',type2='db gene',connection='qtl')
    
    #select the coexpressed genes with the BLAST hits
    ##TK: refactor when full database is available to filter for blast
    ##hits being in either column gene1 or gene2
    selected_coexp = coexp_table.loc[
        coexp_table['gene2'].isin(blast_hits)]
    selected_coexp = selected_coexp.assign(type1='db gene',type2='db gene',connection='coexpression')

    #the second depth of qtl: coexpressed genes with QTL connections
    #check if any of the selected coexpressed genes have phenotype connections
    coexp_qtl_hits = qtl_table.loc[
        (qtl_table['gene'].isin(selected_coexp['gene1'])) | 
        (qtl_table['gene'].isin(selected_coexp['gene2']))
        ]
    coexp_qtl_hits = coexp_qtl_hits.assign(type1='db pheno',type2='db gene',connection='qtl')


    #Cleaning up headers, aligning them to be renamed
    pheno_table = pheno_table[['input pheno','db pheno','similarity','type1','type2','connection']]
    selected_qtl = selected_qtl[['trait','gene','norm_score','type1','type2','connection']]
    selected_coexp = selected_coexp[['gene1','gene2','coexp','type1','type2','connection']]
    coexp_qtl_hits = coexp_qtl_hits[['trait','gene','norm_score','type1','type2','connection']]
    blast_table = blast_table[['query','subject','pid','type1','type2','connection']]
    
    header = ['item1','item2','weight','source1','source2','connection']
    
    #rename the columns
    pheno_table.columns=header
    selected_qtl.columns=header
    selected_coexp.columns=header
    coexp_qtl_hits.columns=header
    blast_table.columns=header
    
    #The results table
    results_table = pd.concat([pheno_table,blast_table,selected_qtl,selected_coexp,coexp_qtl_hits],
                             ignore_index=True)
    results_table = results_table.drop(results_table[(results_table['connection'] == 'text') & 
                                                     (results_table['weight'] < pheno_sim_min)].index)
    num = results_table._get_numeric_data()
    num[num < 0] = 0
    
    
