#!/usr/bin/python
__author__ = 'Jean Vila(jeanccvila@gmail.com)'
__version__ = '1.1.0'

''' Python script for running merging paired sequences from Genewiz 
    Input ab1 sequence format
    Output merged fasta files with all the merged paired sequences
    Used the pairing fasta sequence python code from Joshua E. Goldford, that was modified by Sylvei E. Estrela.
    JV added trimming algorithm (mott modified trimming algorithm) to remove poor quality ends of reads)
    Also updated merge seqs to use higher quality region of sequencing if there is a BP conflict..
    to run simply call python Process_GENEWIZ.py directory_name containing raw reads. NB make sure there is only one copy of each sequence ID. 
    
    20220913 Chang-Yu Chang modified the directory 
    On mac, install the EMBOSS needle following the instruction on http://emboss.open-bio.org/html/use/ch02s07.html
    direct to the directory of EMBOSS and run the following code in command line to install EMBOSS
     
    $ cd DIRECTORY_TO_EMBOSS
    $ ./configure --without-x
    $ make
    
    for the pairwise alignment algorithm `NeedleCommandline`
    
'''
import os
import glob
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO, pairwise2, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib import pyplot as plt
from Bio.Align.Applications import ClustalwCommandline
from Bio.Emboss.Applications import NeedleCommandline


def abi_trim(seq_record,cutoff=0.05):
    """
    Trims the sequence using Richard Mott's modified trimming algorithm (PRIVATE). adapted from Biopython.SeqIO
    Arguments:
        - seq_record - SeqRecord object to be trimmed.
    Trimmed bases are determined from their segment score, which is a
    cumulative sum of each base's score. Base scores are calculated from
    their quality values.
    More about the trimming algorithm:
    http://www.phrap.org/phredphrap/phred.html
    http://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Quality_trimming.html
    """
    start = False  # flag for starting position of trimmed sequence
    segment = 20  # minimum sequence length
    trim_start = 0  # init start index

    if len(seq_record) <= segment:
        return seq_record
    else:
        # calculate base score
        score_list = [
            cutoff - (10 ** (qual / -10.0))
            for qual in seq_record.letter_annotations["phred_quality"]
        ]

        # calculate cumulative score
        # if cumulative value < 0, set it to 0
        # first value is set to 0, because of the assumption that
        # the first base will always be trimmed out
        cummul_score = [0]
        for i in range(1, len(score_list)):
            score = cummul_score[-1] + score_list[i]
            if score < 0:
                cummul_score.append(0)
            else:
                cummul_score.append(score)
                if not start:
                    # trim_start = value when cumulative score is first > 0
                    trim_start = i
                    start = True

        # trim_finish = index of highest cumulative score,
        # marking the end of sequence segment with highest cumulative score
        trim_finish = cummul_score.index(max(cummul_score))

        return seq_record[trim_start:trim_finish]
def calculate_mismatches(seq1,seq2,t=20):
    # seq1 = [x for x in merged_sequence_list if x.id =='225'][0]
    # seq2 = [x for x in merged_sequence_list if x.id =='317'][0]
    # SeqIO.write(merged_sequence_list[0],'../Temp/alpha.faa','fasta')
    # SeqIO.write(merged_sequence_list[2],'../Temp/beta.faa','fasta')
    SeqIO.write(seq1, folder_data + 'temp/14-needle/alpha.faa','fasta')
    SeqIO.write(seq2, folder_data + 'temp/14-needle/beta.faa','fasta')
    needle_cline = NeedleCommandline(asequence = folder_data + 'temp/14-needle/alpha.faa', bsequence = folder_data + 'temp/14-needle/beta.faa', gapopen = 10, gapextend = 0.5, outfile = folder_data + 'temp/14-needle/needle.txt')
    stdout, stderr = needle_cline()
    align = AlignIO.read(folder_data + 'temp/14-needle/needle.txt', 'emboss')
    seq1_aligned = align[0].seq
    seq2_aligned =align[1].seq
    nucleotides = ['A','C','T','G']
    mismatches = [x for x in range(0,len(seq1_aligned)) if seq1_aligned[x] !=seq2_aligned[x] and seq1_aligned[x] in nucleotides and seq2_aligned[x] in nucleotides] 
    mismatches = [x for x  in mismatches if x>t and x< (len(seq1_aligned) - t)] 
    return len(mismatches)

folder_data = str(sys.argv[1])
#folder_data = '/Users/cychang/Dropbox/lab/emergent-coexistence/data/'
#directory = '/Users/cychang/Dropbox/lab/emergent-coexistence/data/raw/Sanger/sanger_seq_16S_communities/'
os.chdir(folder_data + 'raw/sanger/sanger_seq_16S_communities/')

if False:
    merged_sequence_list =[]
    for file in glob.glob("*F.ab1"):
        name = os.path.splitext(os.path.basename(file))[0][:-1]
        var1_f = "F.ab1"
        var1_r = "R.ab1"
        var2 = ".fasta"
        var3 = ".png"
        seqF = abi_trim(SeqIO.read("%s%s" % (name,var1_f),'abi'))
        seqR = abi_trim(SeqIO.read("%s%s" % (name,var1_r),'abi')).reverse_complement()
        alignments = pairwise2.align.localms(seqF.seq,seqR.seq,1,0,-10,-1);
        s1 = alignments[0][0];
        s2 = alignments[0][1];
        s1_phred = seqF.letter_annotations["phred_quality"]
        s2_phred = seqR.letter_annotations["phred_quality"]
        pos_s1 =0
        pos_s2 =0
        #Probabilty of correct call
        s1_prob =[]
        s2_prob =[]
        consensus= ['N' for x in s1]
        posterior_prob =[]
        for i in range(0,len(s1)):
            if s1[i] =='-':
                s1_prob.append(1)
            else:
                s1_prob.append(10**(-s1_phred[pos_s1]/10))
                pos_s1 +=1
    			
            if s2[i] =='-':
                s2_prob.append(1)
            else:
                s2_prob.append(10**(-s2_phred[pos_s2]/10))
                pos_s2 +=1
                
            #Bases agree
            if s1[i] == s2[i]:
                consensus[i] = s1[i]
                px = s1_prob[i]
                py = s2_prob[i]
                posterior_prob.append(((px*py)/3)/(1-px-py+((4/3)*px*py))) #Probability both are wrong given aggreemen see Edgar and Flyvberg 2015
                
            #Bases disagree
            if s1[i] != s2[i]:
                if s2_prob[i] >= s1_prob[i]: #s1 is less likely to be incorrect than s2 so use s1
                    consensus[i] = s1[i]
                    px = s1_prob[i]
                    py = s2_prob[i]
                    posterior_prob.append(px*(1-(px/3))/(px+py-(4*px*py/3))) #Probability both are wrong given aggreement see Edgar and Flyvberg 2015
                if s2_prob[i] < s1_prob[i]: #s2 is likely to be incorrect than s so use s2
                    consensus[i] = s2[i]
                    py = s1_prob[i]
                    px = s2_prob[i]
                    posterior_prob.append(px*(1-(px/3))/(px+py-(4*px*py/3))) #Probability both are wrong given aggreement see Edgar and Flyvberg 2015
            if posterior_prob[i] >0.05:# if probability of merged call being wrong is greater than 0.05 (i.e phred of 10) or lower set to N
                consensus[i] = 'N'
                #posterior_prob[i]=1.0
        consensus_phred = [-10*np.log10(x) for x in posterior_prob]
        newrecord = SeqRecord(Seq("".join(consensus)),id=name[0:3],name="",description="",letter_annotations={'phred_quality':consensus_phred})
        merged_sequence_list.append(newrecord)
        plt.plot(np.log10(posterior_prob))
        plt.xlabel("Position")
        plt.ylabel("log10(P(Error))")
        plt.savefig(folder_data + 'temp/14-needle/' + "%s%s" % (name,var3))
        plt.clf()
    SeqIO.write(merged_sequence_list, folder_data + 'temp/14-needle/Merged_Sequences.fasta', 'fasta')
merged_sequence_list = list(SeqIO.parse(folder_data + 'temp/14-needle/Merged_Sequences.fasta', 'fasta'))
ids = [x.id for x in merged_sequence_list] 
mismatch_matrix = np.zeros((len(merged_sequence_list),len(merged_sequence_list)))
calculate_mismatches(merged_sequence_list[0], merged_sequence_list[1])

print('\nnumber of sequences: ' + str(len(merged_sequence_list)))
for i in range(0,len(merged_sequence_list)):
    print('\n' + str(i))
    for j in range(0,len(merged_sequence_list)):
        print('\t' + str(j))
        mismatch_matrix[i,j] = calculate_mismatches(merged_sequence_list[i],merged_sequence_list[j])

df = pd.DataFrame(mismatch_matrix, index = ids, columns = ids)
df.to_csv(folder_data + 'temp/15-mismatch_matrix_communities.csv')	
    


