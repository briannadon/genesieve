import subprocess
from time import time
import pandas as pd

def blast_query(query,db,out,evalue='1E-10'):
    argstr = f"blastp -outfmt 6 -db {db} -query {query} -out {out} -evalue {evalue}"
    arg = argstr.split(" ")
    return arg
    print(arg)
    
def run_blast(blast_query):
    subprocess.Popen(blast_query)

def process_blast(blast_file,pid_cutoff=65.0):
    columns = ['query',
              'subject',
              'pid',
              'length',
              'mismatch',
              'gapopen',
              'qstart',
              'qend',
              'sstart',
              'send',
              'evalue',
              'bitscore']
    blast = pd.read_table(blast_file,sep='\t',names=columns)
    blast = blast.loc[blast.pid > pid_cutoff]
    blast = blast[['query','subject','pid']]
    blast.pid = blast.pid/100
    blast['connection'] = 'homology'
    return blast
    
