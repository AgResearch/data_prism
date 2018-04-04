#!/bin/env python

import re

#illustrious$ head cluster_repsets.txt
#3382    "EFGW7B_EFN16W8GV_106001-107000"
#3480    "EFG10DR_EFN1FCNZV_260501-261500"
#illustrious$ head probe_genes.txt
#probename       accession       gene
#EFG10DP_EFN1FCNW0_1-1000        WP_064517728.1  hypothetical protein
#EFG10DP_EFN1FCNW0_501-1500      WP_064517728.1  hypothetical protein


# read in the probe annotation
print "%s\t%s"%("cluster","gene")
with open("probe_genes.txt") as annotation:
   split_stream = ( re.split("\t", record.strip()) for record in annotation )
   #print list(split_stream)
   annotation_dict=dict( ( (item[0], item[2]) for item in split_stream ) )
   with open("cluster_repsets.txt") as cluster_repsets:
      for record in cluster_repsets:
         (cluster, accession) = re.split("\t", record.strip())
         accession = accession.replace('"','')
         if accession in annotation_dict:
            annotation = annotation_dict[accession]
            annotation = re.sub("'","prime",annotation)
            print "\"%s\"\t%s"%(cluster, annotation)
