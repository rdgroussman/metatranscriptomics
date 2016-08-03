#!/usr/bin/env python

# Ryan Groussman
# Armbrust Lab, University of Washington
# July 2016

# usage blastp_result.py blastp_output

from os import sys
import matplotlib.pyplot as plt

blastp_output = sys.argv[1]

# hard-coded links to taxa lists:
## CHANGE_ME to match your own file path
jgi_taxa_file = "/Users/rgroussman/Dropbox/Armbrust/bioinfo/scripts/taxa_lists/jgi_taxa.list"
mmetsp_taxa_file = "/Users/rgroussman/Dropbox/Armbrust/bioinfo/scripts/taxa_lists/mmetsp_taxa.list"
# campep2mmetsp_file = "/Users/rgroussman/data/MMETSP/CAMPEP_to_MMETSP.csv" # not used


def process_blast_output(blastp_output):
    """This function will load the given argument for BLASTP output, and then retrieve the e-values
    for each hit attributed to a taxon. Returns the taxa_dict to include the e-value matches for each taxon.
    Note that this is currently configured for just a single-sequence BLASTP query"""

    # load the blastp output
    input_tab = open(blastp_output)

    # build dictionary and list of taxa
    taxa_dict, jgi_taxa, mmetsp_taxa = build_taxa_dict()

    # go through each line of the tabular blastp output
    # at this time, we are only interested in the 1st and 10th (python counting) fields
    for line in input_tab:
        # parse the BLASTP output from JGI and MMETSP header styles
        line_elts = line.split("\t")
        seq_header = line_elts[3]
        db_type = determine_type(seq_header) # Determine whether hit is from JGI or MMETSP
        taxon_id = get_taxon_id(seq_header, db_type)
        evalue = line_elts[2]

        # record results to taxa_dict
        if taxon_id in taxa_dict.keys():
            # print taxon_id, evalue
            taxa_dict[taxon_id].append(evalue)
        else:
            print "WARNING! Not in taxa_dict:", taxon_id

    return taxa_dict

def get_taxon_id(seq_header, db_type):
    """From the seq_header and db_type, return the unique taxon id to which it belongs"""

    if db_type == "jgi":
        return seq_header.split("|")[1]
    elif db_type == "mmetsp":
        return parse_mmetsp_header(seq_header)

def parse_mmetsp_header(seq_header):
    """Returns the taxon id from the NCGR fasta headers.
    Note that /NCGR_PEP_ID= field is ubiquitous and used for this purpose"""

    header_elts = seq_header.split("/")
    ncgr_elts = header_elts[1]
    mmetsp_id, pep_id = get_ncgr_elts(ncgr_elts)
    ncbi_id = header_elts[2][9:]
    # organism = header_elts[3][10:] # not used in this script
    return mmetsp_id

def get_ncgr_elts(ncgr_elts):
    """Example:
    'NCGR_PEP_ID=MMETSP0168-20121206|4096_1 '
    """
    ncgr_elts = ncgr_elts.split("|")
    mmetsp_id = ncgr_elts[0].split("-")
    mmetsp_id = mmetsp_id[0][12:]
    pep_id = ncgr_elts[1].strip()

    return mmetsp_id, pep_id

### We don't have to use this slow-arse function anymore :) ###
# def link_CAMPEP_to_MMETSP_id():
#     """Loads the CSV linking each CAMPEP id to the MMETSP id. Returns as a
#     dictionary with CAMPEP ids as keys and MMETSP id, and seq id as values in a tuple"""
#
#     # CHANGE_ME to match your own file path
#     CAMPEP_to_MMETSP_path = "/Users/rgroussman/data/MMETSP/CAMPEP_to_MMETSP.csv"
#     CAMPEP_to_MMETSP_file = open(CAMPEP_to_MMETSP_path, 'r')
#
#     global CAMPEP_to_MMETSP_dict
#     CAMPEP_to_MMETSP_dict = {}
#     for line in CAMPEP_to_MMETSP_file:
#         line_elts = line.split(",")
#         CAMPEP_to_MMETSP_dict[line_elts[0]] = (line_elts[1], line_elts[2])
#     return CAMPEP_to_MMETSP_dict

def determine_type(seq_header):
    """Determines from a fasta header whether it is from the JGI or MMETSP collections.
    JGI example: jgi|Thaps3|23115|estExt_fgenesh1_pg.C_chr_60346
    MMETSP example: CAMPEP_0172341316
    Returns 'jgi' or 'mmetsp'"""

    if seq_header.startswith("jgi"):
        return "jgi"
    elif seq_header.startswith("CAMPEP"):
        return "mmetsp"

def build_taxa_dict():
    """From a list containing taxonomic identifies (JGI or MMETSP handles),
    create a dictionary with the taxa handles as keys and an empty list as values.
    Return the dictionary."""

    # Load list of taxon identifiers
    jgi_taxa_list = open(jgi_taxa_file, 'r')
    mmetsp_taxa_list = open(mmetsp_taxa_file, 'r')

    # Initialize blank dictionary
    global taxa_dict
    taxa_dict = {}

    # Add JGI IDs to a specific list and the general dictionary
    global jgi_taxa
    jgi_taxa = []
    for taxon in jgi_taxa_list:
        taxon = taxon.strip()
        taxa_dict[taxon] = []
        jgi_taxa.append(taxon)

    # Add MMETSP IDs to a list and the general dictionary
    global mmetsp_taxa
    mmetsp_taxa = []
    for taxon in mmetsp_taxa_list:
        taxon = taxon.strip()
        taxa_dict[taxon] = []
        mmetsp_taxa.append(taxon)

    return taxa_dict, jgi_taxa, mmetsp_taxa

def calculate_percentage_by_evalue():
    """For a series of e-values, calculate the percentage of taxa that have
    at least one hit at that e-value or smaller. Return the number.
    """

    eval_list = []
    for i in range(5, 101, 5):
        eval_list.append(10 ** -i)

    # count total number of taxa
    total_taxa = len(taxa_dict.keys())
    print
    print "Total number of transcriptome samples:", total_taxa
    print

    eval_pcts = []
    print "### Proportion of taxa with at least one BLASTP hit ###"
    for i in eval_list:
        eval_matches = count_matches_at_eval(i)
        print "At e-value " + str(i) + ": " + "{0:.2f}".format(eval_matches/total_taxa)
        eval_pcts.append(eval_matches/total_taxa)
    return eval_list, eval_pcts

def count_matches_at_eval(i):
    """Counts the number of taxa with >=1 match at given eval, returns integer"""

    count = 0.0
    for taxon in taxa_dict.keys():
        counted = False
        if len(taxa_dict[taxon]) > 0:
            for evalue in taxa_dict[taxon]:
                if float(evalue) <= i and counted == False:
                    count += 1
                    counted = True

    return count

def plot_pct_by_evalue(eval_list, eval_pcts):
    """Generates a plot of pct hits versus evalue criteria"""

    x_vals = [str(i)[2:] for i in eval_list]
    plt.bar(range(len(eval_pcts)), eval_pcts) #, align='center')
    plt.xticks(range(len(eval_pcts)), x_vals, size='small')
    plt.xlabel("E-value (exponent)")
    plt.ylabel("Proportion of samples with at least one hit")
    plt.title("Frequency of BLASTP hits by e-value")
    png_filename = blastp_output[:-3] + "eval_pct.png"
    plt.savefig(png_filename)

    plt.clf()
    print "Generated figure for presence ratio by e-value:", png_filename
    print

def count_gene_copies():
    """Determines the frequency of gene copies identified by BLASTP search
    """

    counts_dict = {} # Dictionary for storing count frequencies

    for taxon in taxa_dict:
        n = len(taxa_dict[taxon])
        if n in counts_dict.keys():
            counts_dict[n] += 1
        elif n not in counts_dict.keys():
            counts_dict[n] = 1

    print "### Frequency of copy number ###"
    for i in counts_dict:
        print "Taxa with", i, "copies:", counts_dict[i]

    plot_copy_freq(counts_dict)

def plot_copy_freq(counts_dict):
    """Given a dictionary with frequency of copy numbers, will output a
    PNG file with a histrogram"""

    x = counts_dict.keys()
    y = [counts_dict[i] for i in counts_dict.keys()]
    plt.bar(x, y, align='center')
    plt.xlabel("Number of sequences found")
    plt.ylabel("Frequency")
    plt.title("Number of copies found per taxon at e-value 1e-05")
    png_filename = blastp_output[:-3] + "copy_freq.png"
    plt.savefig(png_filename)
    plt.clf()
    print "Generated figure for copy frequency:", png_filename
    print

def main():
    """Main script"""

    taxa_dict = process_blast_output(blastp_output)

    eval_list, eval_pcts = calculate_percentage_by_evalue()
    plot_pct_by_evalue(eval_list, eval_pcts)
    count_gene_copies()

if __name__ == "__main__":
    main()



## Pseudo-code to be implemented
# 	calculate percentage of taxa represented by at least one hit at intervals of e-value (or other metric):
# 		for each gradation of e-value, count percentage of taxa with at least one hit
# 			gradations: (e-val = 5, 10, 15 etc > > to >100)
# add percentage to list for gradation

## Pseudo-code
# Also graph the number of gene copies per organism - histogram?
