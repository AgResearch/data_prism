#!/bin/bash

cat heatmap_probsets.dat | ./get_cluster_repsets.py > cluster_repsets.txt

awk '{print $2}' cluster_repsets.txt | sed 's/"//g' - > cluster_repset_names.txt

cat probes.fasta | ./fasta_yank.py cluster_repset_names.txt > cluster_repset_names.seq 

tardis.py -hpctype condor blastx -query _condition_fasta_input_cluster_repset_names.seq -num_threads 2 -db /dataset/blastdata/active/mirror/bacteria.nonredundant_protein -evalue 1.0e-6 -max_target_seqs 1 -outfmt  \'7 qseqid sseqid pident evalue staxids sscinames scomnames sskingdoms stitle\' -out _condition_text_output_cluster_repset_names.seq.blastx

gunzip -c cluster_repset_names.seq.blastx.gz | grep -A 1 "hits found" | grep -v "hits found" | grep -v "# BLASTX 2.2.28" | grep -v "# BLAST processed" | grep -v '\-\-' | ./parse_annotation.py > probe_genes.txt

./get_cluster_genes.py > cluster_genes.txt
### note ###
# Re "Misunderstood parameter of NCBI BLAST impacts the correctness of bioinformatics workflows", Nidhi Shah  Michael G Nute  Tandy Warnow  Mihai Pop
# and our use of "-max_target_seqs 1"
#
# Here we have documented our used of this parameter to extract from the database a "similar" (-e 1.0-6)  sequence to the query, so as to
# provide a more meaningful label in a plot, (name of probable gene),  than the EST name; it is not essential for our purpose here
# that the label is based on the best (i.e. lowest evalue) hit in the database; only that it is based on a sequence which is likely to be
# expressed from the same gene or gene family member.
