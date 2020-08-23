#!/usr/bin/python

''' Python script for running merging paired sequences from Genewiz 
    Input ab1 sequence format
    Output merged fasta files with all the merged paired sequences
    Used the pairing fasta sequence python code from Joshua E. Goldford
'''
import os, glob, sys
from Bio import SeqIO

# go to the directory containing the raw sequence
#os.chdir("./Users/sylvie/Desktop/")
os.chdir("/Users/sylvie/Desktop/SangerSeq_Genewiz")

os.chdir("Data/")
# for record in SeqIO.parse("1-16S-rRNA-F.ab1","abi"):
#     print("%s %s" %(record.id,record.seq))
#     print(record)

# convert *.ab1 files to *.fasta
for file in glob.glob("*.ab1"):
    #print(file)
    name = os.path.splitext(os.path.basename(file))[0]
    #print(name)
    var1 = ".ab1"
    var2 = ".fasta"
    count = SeqIO.convert("%s%s" % (name,var1), "abi", "%s%s" % (name,var2), "fasta")
    print("Converted %i records" % count)

# merged fasta file and extract file names to save
for file in glob.glob("*F.fasta"):
    name = os.path.splitext(os.path.basename(file))[0]
    save_name = name.split("-F")[0]
    print(save_name)
    reverse_file = name.split("F")[0]+'R.fasta'
    sys.argv = ['../../mergeGenewizSeqs.py', '-i', save_name, '-f', file, '-r', reverse_file, '-o', save_name+"-merged.fasta"]
#  execfile('../../mergeGenewizSeqs.py')
    exec(open('/Users/sylvie/Desktop/SangerSeq_Genewiz/mergeGenewizSeqs.py').read())

# Terminal
# cat *-merged.fasta > merged.fasta
# submit to RDP classifier together: http://rdp.cme.msu.edu/classifier/hierarchy.jsp;jsessionid=00A6C6E30C2828C56401B137A8527881.radiant
