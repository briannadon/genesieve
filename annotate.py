import subprocess

def augustus_query(fasta,species,outfile,partial=True):
    if partial:
        model = "partial"
    else:
        model = "complete"
    argstr = "augustus --outfile={} --species={} --genemodel={} {}".format(
        outfile,species,partial,fasta)
    return argstr.split(" ")

def run_augustus(query,output):
    subprocess.run(query,stdout=output)


def get_proteins(basename,get_fasta_script='helpers/get-fasta.sh'):
    args = ['bash',get_fasta_script,'-p',basename + ".gff"]
    protfile = basename + ".aa"
    subprocess.run(args)
    return protfile
        
