#!/bin/bash

# Ryan Groussman
# Armbrust Lab, University of Washington
# July 2016

# usage: mmetsp-jgi_blastp.sh ref_handle
# Runs in a folder containing a fasta file, uses BLASTP to search JGI and MMETSP databases.
# Example: if you have mygene.fasta, call script with mmetsp-jgi_blastp.sh mygene
# Input sequence must be in protein space!

# check for the reference fasta file, quits if not found
if [ ! -f $1.fasta ]
  then
    echo "$1.fasta does not exist"
    exit
fi

# create blastp_out/ if it doesn't exit
if [ ! -d blastp_out/ ]
 then
  mkdir blastp_out/
  echo "Created output directory: blastp_out/"
fi

# BLASTP fasta against JGI gene models
# CHANGE_ME to match your own path for blastdb
echo "Running BLASTP of $1.fasta against 172,441 JGI gene models..."
blastp -db /Users/rgroussman/data/JGI/blastdb/JGI_filtered_aa -query $1.fasta -out blastp_out/$1.vsJGI.tab -outfmt "6 qseqid sseqid evalue stitle" -evalue 1e-05
echo "Done: blastp_out/$1.vsJGI.tab"

# BLASTP fasta against MMETSP
# CHANGE_ME to match your own path for blastdb
echo "Running BLASTP of $1.fasta against 16,228,928 MMETSP gene models..."
blastp -db /Users/rgroussman/data/MMETSP/blastdb/MMETSP_aa -query $1.fasta -out blastp_out/$1.vsMMETSP.tab -outfmt "6 qseqid sseqid evalue stitle" -evalue 1e-05
echo "Done: blastp_out/$1.vsMMETSP.tab"

# concatenate results
echo "Concatenating results to blastp_out/$1.vsMMETSP-JGI.tab"
cat blastp_out/$1.vsJGI.tab blastp_out/$1.vsMMETSP.tab > blastp_out/$1.vsMMETSP-JGI.tab
