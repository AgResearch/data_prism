#!/bin/bash

#tardis.py -hpctype condor blastn -query _condition_fasta_input_agbovine.seq.masked -task blastn -num_threads 2 -db /dataset/blastdata/active/species/bt.fna -evalue 1.0e-6 -max_target_seqs 1 -outfmt  \'7 qseqid sseqid pident evalue staxids sscinames scomnames sskingdoms stitle\' -out _condition_text_output_agbovine.seq.masked.blast 

gunzip -c agbovine.seq.masked.blast.gz | grep -A 1 "hits found" | grep -v "hits found" | grep -v "# BLASTN 2.2.28" | grep -v "# BLAST processed" | grep -v '\-\-' | ./parse_annotation.py




