#!/usr/bin/python
import os
import sys
import numpy as np
import pandas as pd
import operator
from Bio import SeqIO
from collections import Counter

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(map(operator.ne, str1, str2))

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "NNK":"X"}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def reading_ref(reffile):
  records = SeqIO.parse(reffile,"fasta")
  ref_dict = {}
  for record in records:
    ID  = str(record.id)
    seq = str(record.seq)
    ref_dict[ID] = seq
  return ref_dict

def analyze_merged_reads(fastq_file, ref_seq, flank_5, flank_3, outfile):
  print ('Analyzing %s' % fastq_file)
  records = SeqIO.parse(fastq_file, "fastq")
  mut_count = []
  for record in records:
    seq = str(record.seq)
    if seq[0:len(flank_5)] == flank_5 and seq[-len(flank_3):] == flank_3:
      seq = seq[len(flank_5):-len(flank_3)]
      if len(seq)==144:
        pep = (translation(seq))
        mut_num = hamming(pep, ref_seq)
        mut_num = 3 if mut_num >= 3 else mut_num
        mut_count.append(mut_num)
  outfile = open(outfile, 'w')
  outfile.write('num_mut'+"\t"+'read_count'+"\n")
  for mutaa, count in Counter(mut_count).items():
    outfile.write(str(mutaa)+"\t"+str(count)+"\n")
  outfile.close()
  
def wrapper(fastq_file, outfile):
  ref_dict  = reading_ref('Fasta/FP_ref.fa')
  ref_ID    = 'FP_ref'
  ref_seq  =  ref_dict[ref_ID]
  flank_5 = 'TTTGGTGGTTTTAATTTTTCACAAATATTACCA'
  flank_3 = 'AACGGCCTTACTGTTTTGCCACCTTTGCTCACAGATGAAATGATTGCTCAA'
  analyze_merged_reads(fastq_file, ref_seq, flank_5, flank_3, outfile)

def main():
  wrapper("fastq_merged/Rep1_ipt_merged.fastq",'result/Lib1_mut_count.tsv')
  wrapper("fastq_merged/Rep2_ipt_merged.fastq",'result/Lib2_mut_count.tsv')

if __name__ == "__main__":
  main()
