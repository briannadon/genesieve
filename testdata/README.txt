# Explanation of files

height_proteins.fa:
This is a multi-FASTA file, the result of an AUGUSTUS (ab initio, i.e. sequence-only, annotation) of the user's input sequence.
It is the "query" genes the user has given us. They give a fasta file with the sequence they mapped, we find genes on it. That's what these are.

height_filtered.blast:
This is the result of BLASTing the heigh_proteins.fa against all proteins (genes) in the rice (oryza sativa) genome.
The columns are [query gene,database gene,% identity (weight), connection].

Oryza_sativa.faa*:
these are just all the proteins in the rice genome from Phytozome. nothing special.

ph_dists.csv:
This is a table of the gensim (https://radimrehurek.com/gensim/) neural-network-derived distance of the input trait "plant height"
to every trait in the database. Each line is a database trait, followed by the cosine distance (0 to 1) of that trait to "plant height".

rice_coexp_dist.csv:
This is the giant coexpression distance matrix. Each row is a gene, each column is a gene. 
**THE NAMES NEED TO BE CHANGED!** they are currently e.g. "LOC_Os11g19770.1" but they SHOULD be "LOC_Os11g19770" -- Remove the trailing ".digit"

rice_height_pheno_source.txt:
this just explains where i got this data from. ignore.

rice_height_qtl.txt:
This is the text description of the phenotype. You can probably ignore it.

rice_qtl_genes.csv:
This is the table that associates previously-known rice genes to traits. 
the "norm_score" column is the "weight" of the connection between that gene and that trait.
the "species" column is probably unnecessary.


# explanation of process

The idea is this: every "hit" consists of two items from the database or user query, a weight between them, and the type of connection.
The weight is always 0 to 1. The connection can be any of ['homology','qtl','text','coexpression']
examples:
[user_gene,database_gene,% identity,"homology"] -- this is a BLAST result so it's "homology" and the "weight" is % identity.
[database_gene,database_trait,norm_score,"qtl"] -- this is a gene previously known to be associated with a trait in our database, so it's "qtl".

# Pseudocode:

for each gene in "height_proteins.fa":
    get all the BLAST hits for that gene from "height_filtered.blast"
    for each BLAST hit that matches a gene in "rice_qtl_genes.csv" (i.e. "qtl_genes"):
        return a value [qtl_gene,trait,norm_score]
        then, note each trait that matched and:
            look up user_trait <-> database_trait, find the NLP model's text distance between them (distances already computed)
            return as [user_trait,database_trait_hit,text_distance] from ph_dists.csv
        then, find the coexpression of each BLAST hit gene with ALL rice genes from rice_coexp_dist.csv:
            filter for genes that have >0.6 coexpression with this subject gene
            report all these connections as [subject_gene,coexpressed_gene,coexpression]
    return a table with all the rows resulting as [item1,item2,weight,connection_type]
        
        