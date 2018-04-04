#!/usr/bin/python2.6
#
# filters a fasta file (from standard input), pulling out sequences by name 
#
from Bio import SeqIO
import sys
 
usage = "usage :  cat some_file.fa | fasta_yank.py file_of_names "
 
if len(sys.argv) != 2:
   print usage
   sys.exit(1)
 
SeqIO.write ((r for r in SeqIO.parse(sys.stdin, "fasta") if r.name in [name.strip() for name in open(sys.argv[1],"r")]) , sys.stdout, "fasta")

#for r in SeqIO.parse(sys.stdin, "fasta"):
#   print r.name
 
