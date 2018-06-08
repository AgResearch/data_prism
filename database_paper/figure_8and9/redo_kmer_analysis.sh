#!/bin/bash
RUN=160531_D00390_0253_AC9AB9ANXX


cd /dataset/hiseq/scratch/postprocessing/$RUN.gbs/SQ0223.processed_sample/uneak/kmer_analysis/; /dataset/hiseq/active/bin/dev/DECONVQC/kmer_entropy.py -b /dataset/hiseq/scratch/postprocessing/$RUN.gbs/SQ0223.processed_sample/uneak/kmer_analysis/ -t zipfian -k 6 -p 1 -o /dataset/hiseq/scratch/postprocessing/$RUN.gbs/SQ0223.processed_sample/uneak/kmer_analysis//kmer_summary.txt  -x /dataset/hiseq/active/bin/dev/DECONVQC/cat_tag_count.sh /dataset/hiseq/scratch/postprocessing/$RUN.gbs/SQ0223.processed_sample/uneak/KGD/../tagCounts/*.cnt 
