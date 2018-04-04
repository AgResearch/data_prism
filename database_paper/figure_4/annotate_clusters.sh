#!/bin/bash

cat heatmap_probsets.dat | ./get_cluster_repsets.py > cluster_repsets.txt

awk '{print $2}' cluster_repsets.txt | sed 's/"//g' - > cluster_repset_names.txt

cat probes.fasta | ./fasta_yank.py cluster_repset_names.txt > cluster_repset_names.seq 

tardis.py -hpctype condor blastx -query _condition_fasta_input_cluster_repset_names.seq -num_threads 2 -db /dataset/blastdata/active/mirror/bacteria.nonredundant_protein -evalue 1.0e-6 -max_target_seqs 1 -outfmt  \'7 qseqid sseqid pident evalue staxids sscinames scomnames sskingdoms stitle\' -out _condition_text_output_cluster_repset_names.seq.blastx

gunzip -c cluster_repset_names.seq.blastx.gz | grep -A 1 "hits found" | grep -v "hits found" | grep -v "# BLASTX 2.2.28" | grep -v "# BLAST processed" | grep -v '\-\-' | ./parse_annotation.py > probe_genes.txt

./get_cluster_genes.py > cluster_genes.txt




