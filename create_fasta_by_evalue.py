#!/usr/bin/env python

# Ryan Groussman
# Armbrust Lab, University of Washington
# July 2016

# usage create_fasta_by_evalue.py infile.tab outfile.fasta e-value
# e-value in format: 1e-10

# load packages
from os import sys
from Bio import SeqIO

infile_path = sys.argv[1]
outfile_path = sys.argv[2]
lowest_eval = float(sys.argv[3])

# hard-coded links: CHANGE_ME to match your own paths
jgi_fasta = "/Users/rgroussman/data/JGI/JGI_microeuks.aa.fasta"
mmetsp_fasta = "/Users/rgroussman/data/MMETSP/CAM_P_0001000.pep.fa"

def determine_type(seq_header):
    """Determines from a fasta header whether it is from the JGI or MMETSP collections.
    JGI example: jgi|Thaps3|23115|estExt_fgenesh1_pg.C_chr_60346
    MMETSP example: CAMPEP_0172341316
    Returns 'jgi' or 'mmetsp'"""

    if seq_header.startswith("jgi"):
        return "jgi"
    elif seq_header.startswith("CAMPEP"):
        return "mmetsp"

def add_seq_id_to_list(seq_id):
    """Will add a sequence id to lists for either MMETSP or JGI databases"""

    db_type = determine_type(seq_id)

    if db_type == "jgi":
        jgi_getlist.add(seq_id)
    elif db_type == "mmetsp":
        mmetsp_getlist.add(seq_id)

def extract_jgi_seqs(jgi_getlist):
    """Extracts sequences from the linked JGI fasta using the list"""

    for seq_record in SeqIO.parse(jgi_fasta, "fasta"):
        if seq_record.id in jgi_getlist:
            SeqIO.write(seq_record, output_fasta, "fasta")
            jgi_getlist.remove(seq_record.id)

def extract_mmetsp_seqs(mmetsp_getlist):
    """Extracts sequences from the linked MMETSP fasta using the list"""

    for seq_record in SeqIO.parse(mmetsp_fasta, "fasta"):
        if seq_record.id in mmetsp_getlist:
            SeqIO.write(seq_record, output_fasta, "fasta")
            mmetsp_getlist.remove(seq_record.id)

# initialize blank sets
jgi_getlist = set([])
mmetsp_getlist = set([])

# load input file and output file
input_tab = open(infile_path, 'r')
output_fasta = open(outfile_path, 'w')

# go through the input tab, adding seq id if it passes eval threshold
for line in input_tab:
    line_elts = line.split("\t")
    seq_id = line_elts[1]
    evalue = float(line_elts[2])
    if evalue <= lowest_eval:
        add_seq_id_to_list(seq_id)

# now run through the large FASTA files and extract the sequences
print "Extracting", len(jgi_getlist), "sequences from", jgi_fasta
extract_jgi_seqs(jgi_getlist)
print "Sequences unclaimed:", len(jgi_getlist)
print
print "Extracting", len(mmetsp_getlist), "sequences from", mmetsp_fasta
extract_mmetsp_seqs(mmetsp_getlist)
print "Sequences unclaimed:", len(mmetsp_getlist)
