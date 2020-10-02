#!/usr/bin/python

''' Python script for running merging paired sequences from Genewiz 

'''
import argparse
from Bio.Seq import Seq
from Bio import pairwise2
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import os.path


__author__ = 'Joshua E. Goldford'

def get_args():
    '''This function parses and return arguments passed in'''
    # Assign description to the help doc
    parser = argparse.ArgumentParser(
        description='Script merges forward and reverse reads together from GeneWiz data')
    # Add arguments
    parser.add_argument(
        '-i', '--seq_header', type=str, help='header for fasta', required=True)
    parser.add_argument(
        '-f', '--forward', type=str, help='forward read file', required=True)
    parser.add_argument(
        '-r', '--reverse', type=str, help='reverse read file',required=True)
    parser.add_argument(
        '-o', '--output', type=str, help='output file', required=False, default='merged.fasta')
    # Array for all arguments passed to script
    args = parser.parse_args()
    # Assign args to variables
    header = args.seq_header
    forward = args.forward
    reverse = args.reverse
    output = args.output
    # Return all variable values
    return header,forward, reverse, output


def merge_seqs(seqF,seqR):
	alignments = pairwise2.align.localms(seqF,seqR,1,0,-10,-1);
	s1 = alignments[0][0];
	s2 = alignments[0][1];
	s1 = list(s1);
	s2 = list(s2);
	s_merged = s1;
	for idx in range(0,len(s1)):
		#s_merged[idx] = s1[idx];
		if s1[idx] == '-':
			s_merged[idx] = s2[idx];
		if s2[idx] == '-':
			s_merged[idx] = s1[idx];
		if s1[idx] == 'N' and s2[idx] in ['A','C','G','T']:
			s_merged[idx] = s2[idx];
		if s2[idx] == 'N' and s1[idx] in ['A','C','G','T']:
			s_merged[idx] = s1[idx];
	return ''.join(s_merged)


header,forward,reverse,output = get_args()

# get forward read
r1 = SeqIO.parse(forward, "fasta");
r1 = list(r1);
r1 = r1[0];

# get reverse read
r2 = SeqIO.parse(reverse, "fasta");
r2 = list(r2);
r2 = r2[0];

# merge sequences
s = merge_seqs(r1.seq,r2.seq.reverse_complement());

# make new seq io record
record = SeqRecord(Seq(s,IUPAC.IUPACAmbiguousDNA()),id=header,name="",description="");

handle = open(output,"w");
fasta_out = FastaIO.FastaWriter(handle, wrap=None);
SeqIO.write(record,handle,"fasta")
handle.close()


