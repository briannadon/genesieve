#!/bin/bash
PROG="${0##*/}"

help() {
    echo -e "$PROG [options] <Agustus gff>\n"
    echo "Options:"
    echo "-c <coding file>    Output file for CDS fasta"
    echo "-p <coding file>    Output file for protein fasta"
    exit 0
}

while getopts "c:p:h" arg ; do
    case "$arg" in
        c)  cfile="$OPTARG" ;;
        p)  pfile="$OPTARG" ;;
        *)  help ;;
    esac
done

shift $((OPTIND-1))

# pre-process the file to "join" all coding and protein sequences
# output will be
# >gene1
# coding <coding-sequence>
# protein <protein-sequence>
# coding <coding-sequence>
# protein <protein-sequence>
# ....
# >gene2
# ...
awk '/^# start gene g/ {
    print ">" $NF
    f = ""
}

/^# (coding|protein) sequence/ {
    f = $NF                     # store the start of sequence
    t = $2                      # t is either "coding" or "protein"
    if ($NF ~ /\]$/) {          # for one line sequences, ] is on the same line, so print it out
        print t, f; f = ""
    }
}

length(f) && $1 == "#" && NF == 2 {    # multiline sequences will continue appending to f
    f = f $NF                          # first field is '#' second is sequence,
    if ($NF ~ /\]$/) {                 # once we get to ']' thats the end of sequence, so print it out
        print t, f; f = ""
    }
} ' $1  | tr -d '[]' > preprocessed

# after that Pick the gene line and the first two coding and protein lines
# and write then out
awk -v pfile="$pfile" -v cfile="$cfile" 'BEGIN {
    pfile = !length(pfile) ? "protein.fasta" : pfile;
    printf "" > pfile
    cfile = !length(cfile) ? "coding.fasta" : cfile;
    printf "" > cfile
}

NF == 1 && /^>/ {
    gene = $1
    c = p = 0
    total++
}

$1 == "coding" && c == 0 && NF == 2 {
    print gene "\n" $NF >> cfile
    c++; codings++;
}

$1 == "protein" && p == 0 && NF == 2 {
    print  gene "\n" $NF >> pfile
    p++; proteins++;
}

$1 == "protein" && (c == 0 || p == 0) {
    if (c == 0) {
        missing_codings++;
    } else {
        missing_proteins++;
    }
    print "WARNING: Gene", substr(gene, 2), "has", c, "coding sequences and", p, "protein sequences"
}

END {
    print "Statistics:\nTotal Genes:", total, "\nTotal Coding Sequences:", codings, "\nTotal Protein Sequences:", proteins, "\nEmpty Coding Sequences:", missing_codings, "\nEmpty Protein Sequences:", missing_proteins
}' preprocessed

echo "Total genes in $1: $(grep -c '^# start gene' $1)" 
