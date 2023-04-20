#!/usr/bin/python
import os
import sys
from collections import defaultdict

def read_fitfile(file_file):
  infile = open(file_file, 'r')
  fit_dict = {}
  for line in infile.readlines():
    if 'pos' in line: continue
    line = line.rstrip().rsplit()
    pos  = line[1]
    fit  = line[3]
    fit_dict[pos] = float(fit)
  infile.close()
  return (fit_dict)

def normalizing_fit(fit_dict):
  max_fit  = max(fit_dict.values())
  min_fit  = min(fit_dict.values())
  print ('fit range (pre-norm): %f to %f' % (min_fit, max_fit))
  norm_fit = []
  for pos in fit_dict.keys():
    fit = fit_dict[pos]
    fit = (fit-min_fit)/(max_fit-min_fit)*100
    fit_dict[pos] = fit
    norm_fit.append(fit)
  print ('fit range (post-norm): %f to %f' % (min(norm_fit), max(norm_fit)))
  return (fit_dict)

def add_fit_to_pdb(fit_dict, pdb_file):
  assert('.pdb' in pdb_file)
  new_pdb_file = pdb_file.replace('.pdb', '_fit.pdb')
  print ("writing: %s" % new_pdb_file)
  infile  = open(pdb_file, 'r')
  outfile = open(new_pdb_file,'w')
  for line in infile.readlines():
    if "ATOM" == line[0:4]:
      pos      = int(line[23:26])
      chain    = line[21:23].rstrip()
      b_factor = line[60:66]
      if chain == 'A' and str(pos) in fit_dict.keys():
        fit = str(round(fit_dict[str(pos)],2))
        fit = fit+'0'*(2-len(fit.rsplit('.')[-1]))
        fit = ' '*(6-len(fit))+fit
        new_line = line[0:60]+fit+line[66::]
        outfile.write(new_line)
      else:
        outfile.write(line)
    else:
      outfile.write(line)
  outfile.close()
  

def main():
  pdb_file = "PDB/7my8.pdb"
  fit_file = 'result/FP_DMS_fit_by_resi.tsv'
  fit_dict = read_fitfile(fit_file)
  fit_dict = normalizing_fit(fit_dict)
  add_fit_to_pdb(fit_dict, pdb_file)

if __name__ == "__main__":
  main() 
