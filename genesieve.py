import subprocess
import time

from sys import argv

import blast
import sanitize
import scoring
import phenotype
import annotate

def pheno_input(phenofile='tests/testdata/rice_qtl_genes.csv'):
    with open(phenofile) as f:
        phenotype = f.read()
    return phenotype

'''def merge_qtl_blast(qtl,blast):
    merged_df = pd.merge(left=qtl,right=blast,left_on="gene",right_on="gene")
    merged_df = merged_df[['query','gene','species','pid','trait','norm_score']]
    return merged_df'''
    
augustus_protein_script = "/software/7/apps/augustus/3.3.3/scripts/getAnnoFasta.pl"
coexp_min = 0.65

if __name__=="__main__":
    script, in_fasta, in_pheno = argv
    timestamp = str(time.time())
    dbheader = ['item1','item2','weight','connection']
    
    #run annotation first
    augout = timestamp+".augustus"
    aquery = annotate.augustus_query(augout,in_fasta,species="rice")
    a = subprocess.Popen(aquery)
    
    #While augustus is running, load up the phenotype stuff
    pheno_model = 'data/models/dbow_wiki_2.model'
    
    #Load all our QTL data (genes:phenotypes previously known)
    #TODO change this to SQL
    pheno_list = "SELECT trait FROM QTL_TABLE".tolist()
    
    # sanitize pheno input:
    pheno_text = pheno_input(in_pheno)
    pheno_text = phenotype.sanitize_text(pheno_text)
    
    pheno_dists = phenotype.get_distances(pheno_model,trait_file,pheno_text)
    distlist =((pheno_text,p.description,p.distance,'text') for p in pheno_dists)
    
    #delete the memory-hungry model
    del pheno_model
    
    pheno_table = pd.DataFrame(distlist,columns=dbheader)
    
    while not a.poll():
        time.sleep(5)
    
    aa_file = annotate.get_proteins(timestamp,perl_script)
 
    #Get the BLAST running in the background as "b"
    blast_output = timestamp + ".blast"
    subprocess.run(blast.blast_query(aa_file,
                                      'data/testdata/Orzya_sativa.faa',
                                      blast_output
                                     ))
    
    #if the BLAST is still running, wait till it finishes
    #while not b.poll():
    #    time.sleep(5)

    blast_table = blast.process_blast(blast_output)    
    blast_hits = blast_table.subject.tolist()
    
    #filtered_qtl_table  = qtl_table.loc[qtl_table.subect.isin(blast_hits)]
    
    qtl_blast_query = "SELECT trait, gene, normscore FROM QTL_TABLE WHERE gene IN ({}) and ;".format(','.join(blast_hits))
    coexp_query = "SELECT gene1, gene2, coexp FROM COEXP_TABLE WHERE gene1 IN ({}) AND coexp > {};".format(','.join(blast_hits),coexp_min)
    qtl_coexp_query = "SELECT gene1, gene2, 