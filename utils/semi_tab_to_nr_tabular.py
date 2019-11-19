#!/bin/env pypy
from __future__ import print_function
import itertools
import sys
import argparse
import re

def blast_result_iter():
    """ input is like this :

# BLASTN 2.7.1+
# Query: D00390:494:CDT5UANXX:1:1107:19400:10481 1:N:0:0
# Database: /dataset/GBS_Tcirc/ztmp/SQ1014_analysis/db/microbial_references
# 0 hits found
# BLASTN 2.7.1+
# Query: D00390:494:CDT5UANXX:1:1107:19667:10325 1:N:0:0
# Database: /dataset/GBS_Tcirc/ztmp/SQ1014_analysis/db/microbial_references
# Fields: query id, subject id, % identity, evalue
# 1 hits found
D00390:494:CDT5UANXX:1:1107:19667:10325 NC_018708.1     95.652  1.19e-34
# BLASTN 2.7.1+
# Query: D00390:494:CDT5UANXX:1:1107:21123:10274 1:N:0:0
# Database: /dataset/GBS_Tcirc/ztmp/SQ1014_analysis/db/microbial_references
# Fields: query id, subject id, % identity, evalue
# 1 hits found
D00390:494:CDT5UANXX:1:1107:21123:10274 NZ_CP017311.1   96.970  1.18e-23
# BLASTN 2.7.1+
# Query: D00390:494:CDT5UANXX:1:1107:1439:10676 1:N:0:0
# Database: /dataset/GBS_Tcirc/ztmp/SQ1014_analysis/db/microbial_references
# Fields: query id, subject id, % identity, evalue
# 2 hits found
D00390:494:CDT5UANXX:1:1107:1439:10676  NZ_BCUS01000017.1       96.875  1.86e-22
D00390:494:CDT5UANXX:1:1107:1439:10676  NZ_BCUS01000017.1       92.647  1.12e-19
# BLASTN 2.7.1+

"""
    record_iter = (record.strip() for record in sys.stdin if len(record.strip()) > 0)
    hit_group_iter = itertools.groupby(record_iter, lambda record: {True:"name", False:"hits"}[record[2:8] == "Query:"])
    name=None
    for (group, records) in hit_group_iter:
        if group == "name":
            name=re.split("\s+",records.next())[2]
        else:
            for record in records:
                match=re.match("#\s[123456789]{1}[0123456789]*\shits", record)
                if match is not None:
                    yield (records.next())
                    break

def run():
    parsed_hits = blast_result_iter()
    for first_hit_record in parsed_hits:  
        print(first_hit_record)

def get_options():
    # no args - this here just to provide help interface , semi_tab_to_nr_tabular.py -h 
    description = """
input is via stdin , and is a stream of semi-tabular NCBI blast output (e.g. from something like -outfmt '7 qseqid sseqid pident evalue staxids sscinames scomnames sskingdoms stitle')
, output (stdout) is "non-redundant tabular" (non-redundant = only 1 hit record reported per query. queries with no hits not reported)
    """

    long_description = """
    Example :

gunzip -c myresults.blast.gz | ./semi_tab_to_nr_tabular.py
gunzip -c /somefolder/*.results.gz | ./semi_tab_to_nr_tabular.py
blastn -query myseq -db nt -outfmt '7 qseqid sseqid pident evalue staxids sscinames scomnames sskingdoms stitle' | ./semi_tab_to_nr_tabular.py

(The results can be spot checked by comparing with something like 
 gunzip -c /somefolder/*.results.gz | egrep --no-group-separator -A 1 "hits found" | grep -v BLASTN | grep -v "hits found" | grep -v BLAST )


    """
    parser = argparse.ArgumentParser(description=description, epilog=long_description, formatter_class = argparse.RawDescriptionHelpFormatter)
    args = vars(parser.parse_args())

def main():
    get_options()
    run()

if __name__ == "__main__":
   main()

