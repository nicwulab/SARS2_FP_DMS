#!/usr/bin/python
import os
import sys
import glob
from Bio import SeqIO
from collections import Counter

def ProcessMultilib(infile, sampleID, flank_5, flank_3, var_dict, length_roi):
    print ("Reading %s" % infile)
    records = SeqIO.parse(infile,"fastq")
    variants = [] 
    record_count = 0
    for record in records:
        record_count += 1
        Rseq  = record.seq
        Rroi = Rseq
        # Only include reads with correct forward/reverse primers and the correct length
        if Rroi[0:len(flank_5)] == flank_5 and Rroi[-len(flank_3):] == flank_3:
            Rroi = Rroi[len(flank_5):-len(flank_3)] # Trim forward and reverse primers
            if len(Rroi) == length_roi:
                variants.append(Rroi)
        #if record_count == 5000: break
    var_dict[sampleID] = Counter(variants)
    return var_dict 

def writing_file(outfile, var_dict, sampleIDs, muts):
    print ("Compiling results into %s" % outfile)
    outfile = open(outfile,'w')
    outfile.write('variant'+"\t"+"\t".join(sampleIDs)+"\n")
    for mut in muts:
        out = [mut]
        for sampleID in sampleIDs:
            count = 0 if mut not in var_dict[sampleID].keys() else var_dict[sampleID][mut]
            out.append(count)
        outfile.write("\t".join(map(str,out))+"\n")
    outfile.close()

def Output(var_dict, outfile):
    sampleIDs = sorted([sampleID for sampleID in list(var_dict.keys())])
    muts    = sorted(list(set([mut for sampleID in sampleIDs for mut in var_dict[sampleID].keys()])))
    writing_file(outfile, var_dict, sampleIDs, muts)

def main():
    outfile = 'result/FP_DMS_count_nuc.tsv'
    infiles = glob.glob('fastq_merged/*merged_assembled.fastq')
    length_roi = 144
    flank_5 = 'TTTGGTGGTTTTAATTTTTCACAAATATTACCA'
    flank_3 = 'AACGGCCTTACTGTTTTGCCACCTTTGCTCACAGATGAAATGATTGCTCAA'
    var_dict = {}
    for infile in infiles:
        sampleID = '_'.join(infile.rsplit('/')[1].rsplit('_')[0:2])
        var_dict = ProcessMultilib(infile, sampleID, flank_5, flank_3, var_dict, length_roi)
    Output(var_dict, outfile)

if __name__ == "__main__":
  main()
