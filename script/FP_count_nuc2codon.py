#!/usr/bin/python
import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from functools import partial
from collections import defaultdict

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

def reading_barcode(barfile):
  print ('reading: %s' % barfile)
  infile = open(barfile, 'r')
  bc_dict = {}
  for line in infile.readlines():
    if 'barcode' in line: continue
    bc, ID = line.rstrip().rsplit("\t")
    if ID == 'WT': 
      bc_dict[bc] = 'WT'
    else:
      cas, cas_pos = list(map(int,ID.replace('Cassette','').rsplit('_')))
      pos = (cas-1)*8+cas_pos
      bc_dict[bc] = pos
  return (bc_dict)

def reading_ref(reffile):
  records = SeqIO.parse(reffile,"fasta")
  ref_dict = {}
  for record in records:
    ID  = str(record.id)
    seq = str(record.seq)
    ref_dict[ID] = seq
  return ref_dict

def nuc_2_codon(countfile, bc_dict, ref_dict, ref_ID, roi_length):
  print ('reading: %s' % countfile)
  infile    = open(countfile, 'r')
  count_line = 0
  data_dict  = {}
  for line in infile.readlines():
    count_line += 1
    if count_line == 1:
      header = line.rstrip().rsplit("\t")[1::]
      continue
    line = line.rstrip().rsplit("\t")
    DNA  = line[0]
    data = np.array(list(map(int,line[1::])))
    if len(DNA) != roi_length: continue
    barcode = DNA[2::3] 
    if barcode not in bc_dict.keys(): continue
    mut_pos = bc_dict[barcode]
    pep = translation(DNA)
    mut = 'WT' if mut_pos == 'WT' else str(mut_pos)+'_'+DNA[(int(mut_pos)-1)*3:int(mut_pos)*3]
    if mut not in data_dict.keys():
      data_dict[mut] = data
    else:
      data_dict[mut] += data
  infile.close()
  return header, data_dict

def writing_aa_count(data_dict, ref_dict, ref_ID, positions, header, outfile, offset):
  print ('writing: %s' % outfile)
  outfile   = open(outfile, 'w') 
  header = ['pos','codon']+header
  outfile.write("\t".join(header)+"\n")
  muts = list(data_dict.keys())
  muts.remove('WT')
  for mut in sorted(muts, key=lambda x:int(x.rsplit('_')[0]))+['WT']:
    adj_pos   = 'WT' if 'WT' in mut else str(int(mut.rsplit('_')[0])+offset)
    adj_mut   = 'WT' if 'WT' in mut else mut.rsplit('_')[1]
    data = list(map(str,data_dict[mut])) if mut in data_dict.keys() else list(map(str,[0]*(len(header)-2)))
    outfile.write(adj_pos+"\t"+adj_mut+"\t"+"\t".join(data)+"\n")
  outfile.close()

def main():
  outfile   = 'result/FP_DMS_count_codon.tsv'
  countfile = 'result/FP_DMS_count_nuc.tsv'
  positions  = list(range(808, 856))
  barfile   = 'data/barcodes.tsv'
  ref_dict  = reading_ref('Fasta/FP_ref.fa')
  ref_ID    = 'FP_ref'
  bc_dict   = reading_barcode(barfile)
  roi_length = 144
  offset = 807
  header, data_dict = nuc_2_codon(countfile, bc_dict, ref_dict, ref_ID, roi_length)
  writing_aa_count(data_dict, ref_dict, ref_ID, positions, header, outfile, offset)

if __name__ == "__main__":
  main()
