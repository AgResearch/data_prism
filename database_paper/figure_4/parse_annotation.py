#!/bin/env python 

import sys
import re

print "\t".join(("probename","accession","gene"))
for record in sys.stdin:
   #print record
   fields = re.split("\t",record.strip())
   probe = fields[0]
   accession=fields[1]
   gene=fields[1]
   match = re.search("\S+\s+(?:MULTISPECIES\:\s)*(.+)\[", fields[8])
   if match is not None:
      gene = match.groups()[0]
   print "\t".join((probe, accession, gene))
   
