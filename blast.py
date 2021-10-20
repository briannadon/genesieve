import subprocess
from time import time
import pandas as pd
from sys import argv
import shlex


def blast_query(query,db,out,evalue='1E-10'):
    #argstr = f"blastp -outfmt 6 -db {db} -query {query} -out {out} -evalue {evalue}"
    argstr = (f'blastp -outfmt "6 qseqid sseqid pident ' + 
    'length qcovhsp mismatch gapopen qstart qend sstart send evalue bitscore" ' + 
                f"-db {db} -query {query} -out {out} -evalue {evalue}")
    arg = shlex.split(argstr)
    return arg
    
    
def run_blast(blast_query):
    subprocess.Popen(blast_query, shell=True)

def process_blast(blast_file,hom_cutoff=.65):
    columns = ['query',
              'subject',
              'pid',
              'length',
              'qcovhsp',
              'mismatch',
              'gapopen',
              'qstart',
              'qend',
              'sstart',
              'send',
              'evalue',
              'bitscore']
    with open(blast_file) as b:
        blast_dict = {}
        for line in b:
            line = line.strip().split('\t')
            line_dict = dict(zip(columns,line))
            if (line_dict['query'],line_dict['subject']) not in blast_dict:
                b_score = (float(line_dict['pid']) / 100) * (float(line_dict['qcovhsp']) / 100)
                blast_dict[(line_dict['query'],line_dict['subject'])] = b_score
            else:
                b_score = blast_dict[(line_dict['query'],line_dict['subject'])]
                new_score = (float(line_dict['pid']) / 100) * (float(line_dict['qcovhsp']) / 100)
                blast_dict[(line_dict['query'],line_dict['subject'])] = b_score + new_score
    df_list = []
    for k,v in blast_dict.items():
        df_list.append([k[0],k[1],v])
    
    blast = pd.DataFrame.from_records(df_list,columns = ['query','subject','hom_score'])
    blast = blast.loc[blast['hom_score'] > hom_cutoff]        
    blast['connection'] = 'homology'
    
    return blast

if __name__ == "__main__":
    script, query, db, blastfile = argv
    q = blast_query(query, db, blastfile)    
    run_blast(q)
    process_blast(blastfile)
    
