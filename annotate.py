import subprocess

def augustus_query(fasta,species,outfile,partial=False):
    if partial:
        model = "partial"
    else:
        model = "complete"
    argstr = (f"augustus --outfile={outfile} "
                f"--species={species} --genemodel={model} {fasta}")
    return argstr.split(" ")

def run_augustus(query,output):
    subprocess.run(query,stdout=output)


def get_proteins(basename,get_fasta_script='helpers/get-fasta.sh'):
    args = ['bash',get_fasta_script,'-p',basename + ".faa",basename+".augustus"]
    print("args: ",' '.join(args))
    protfile = basename + ".faa"
    subprocess.run(args)
    return protfile
        
