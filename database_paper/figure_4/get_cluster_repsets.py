#!/bin/env python 

#"x"
#"EFG10DP_EFN1FCNW0_1-1000"      3514
#"EFG10DP_EFN1FCNW0_100001-101000"       7152
#"EFG10DP_EFN1FCNW0_10001-11000" 1275
#"EFG10DP_EFN1FCNW0_1001-2000"   1751
#"EFG10DP_EFN1FCNW0_100501-101500"       3215
#"EFG10DP_EFN1FCNW0_101001-102000"       6523

# simply pick the first sequence met, for each cluster
import sys
import re

repset = {}
for record in sys.stdin:
   fields=re.split("\s",record.strip())
   if len(fields) == 2:
      if fields[1] not in repset:
          repset[fields[1]] = fields[0]

for item in repset.items():
   print "%s\t%s"%(item)


