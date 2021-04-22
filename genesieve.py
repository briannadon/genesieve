import subprocess
import time
from sys import argv
import annotate
import blast
import sanitize
import scoring
import phenotype
import coexpression


def pheno_input(phenofile='tests/testdata/rice_qtl_genes.csv'):
    with open(phenofile) as f:
        phenotype = f.read()
    return phenotype


if __name__=="__main__":
    script, in_fasta, in_pheno = argv
    timestamp = str(time.time())
    
    augustus_protein_script = "helpers/get-fasta.sh"
    coexp_min = 0.65
    pheno_sim_min = 0.5
    header = ['item1','item2','weight','type1','type2','connection']
    
    
    qtl_db = "testdata/rice_qtl_genes.csv"
    coexp_db = "testdata/ricecoexp_sample_1000_LOC_Os05g04990.1.csv"
    blast_file = 'rice_height_test/scripts/Ph5-1_vs_rice.blast'
    pheno_model = 'testdata/model/dbow_wiki_2.model'
    
    #run annotation first
    augout = timestamp+".augustus"
    aquery = annotate.augustus_query(augout,in_fasta,species="rice")
    a = subprocess.Popen(aquery)
    
    #While augustus is running, load up the phenotype stuff
    
    qtl_table = pd.read_csv(qtl_db)
    trait_list = list(set(list(qtl_table.trait)))

    #sanitize pheno input:
    pheno_text = pheno_input(in_pheno)
    pheno_text = phenotype.sanitize_text(pheno_text)

    pheno_dists = phenotype.get_distances(pheno_model,trait_list,pheno_text)
    pheno_dist_d =[(pheno_text,p.description,p.distance,'text') for p in pheno_dists]
    pheno_table = pd.DataFrame(pheno_dist_d,columns=['input pheno','db pheno','similarity','connection'])
    pheno_table = pheno_table.assign(type1='input pheno',type2='db pheno')
    
    #check if annotation is done
    while not a.poll():
        time.sleep(5)
    
    aa_file = annotate.get_proteins(timestamp,bash_script="get-fasta.sh")
 
    #Get the BLAST running in the background as "b"
    blast_output = timestamp + ".blast"
    subprocess.run(blast.blast_query(aa_file,
                                      'data/testdata/Orzya_sativa.faa',
                                      blast_output
                                     ))
    
    #load in the coexp data
    coexp_table = pd.read_csv(coexp_db)
    coexp_table['gene1'] = [re.sub('\.\d','',x) for x in list(coexp_table['gene1'])]
    coexp_table['gene2'] = [re.sub('\.\d','',x) for x in list(coexp_table['gene2'])]
    
    #if the BLAST is still running, wait till it finishes
    while not b.poll():
        time.sleep(5)

    blast_table = blast.process_blast(blast_file)
    blast_table = blast_table.assign(type1='input_gene',type2='gene')
    blast_table['subject'] = [re.sub('\.\d','',x) for x in list(blast_table['subject'])]
    blast_hits = blast_table.subject.tolist()
    
    
    ###
    #Create the Results Table
    ###
    
    #select the QTL that are associated with the BLAST hits
    selected_qtl = qtl_table.loc[
        qtl_table['gene'].isin(blast_hits)]
    selected_qtl = selected_qtl.assign(type1='db pheno',type2='db gene',connection='qtl')
    
    #select the coexpressed genes with the BLAST hits
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


    #Cleaning up headers
    pheno_table = pheno_table[['input pheno','db pheno','similarity','type1','type2','connection']]
    selected_qtl = selected_qtl[['trait','gene','norm_score','type1','type2','connection']]
    selected_coexp = selected_coexp[['gene1','gene2','coexp','type1','type2','connection']]
    coexp_qtl_hits = coexp_qtl_hits[['trait','gene','norm_score','type1','type2','connection']]
    blast_table = blast_table[['query','subject','pid','type1','type2','connection']]
    pheno_table.columns=header
    selected_qtl.columns=header
    selected_coexp.columns=header
    coexp_qtl_hits.columns=header
    blast_table.columns=header
    
    #The results table
    results_table = pd.concat([pheno_table,blast_table,selected_qtl,selected_coexp,coexp_qtl_hits],
                             ignore_index=True)
