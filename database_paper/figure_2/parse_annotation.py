#!/bin/env python 

import sys
import re

# extract accession and gene from records like 
# 000110BPBA008206HT      NM_001075743.1  99.38   1e-160  N/A     N/A     N/A     N/A     Bos taurus eukaryotic translation initiation factor 2B subunit alpha (EIF2B1), mRNA


print "\t".join(("estname","accession","gene_name"))
for record in sys.stdin:
   #print record
   fields = re.split("\t",record.strip())
   est = fields[0]
   accession=fields[1]
   gene=fields[1]
   matches = re.finditer("(\(\S+\))", fields[8])
   if matches is not None:
      match = list(matches)[-1]
      gene = match.groups()[0][1:len(match.groups()[0])-1]
   print "\t".join((est, accession, gene))
   
