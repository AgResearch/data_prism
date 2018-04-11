#!/bin/sh

if [ -z "$GBS_BIN" ]; then
   echo "GBS_BIN not set - quitting"
   exit 1
fi

BUILD_ROOT=/dataset/hiseq/scratch/postprocessing 
THIS_RUN=$1

# ... to do - get the cumulativei list below by doing a find on pickle files 
# - will also need to add current run as an argument. Then add this to the 
# processing script - will call it after the make has finished

# find all the runs with summaries , excluding this run
runs=""
for file in `find /dataset/hiseq/scratch/postprocessing/*.processed/taxonomy_analysis -maxdepth 1 -name samples_taxonomy_table.txt -print`; do
   # e.g. /dataset/hiseq/scratch/postprocessing/141217_D00390_0214_BC4UEHACXX.processed/taxonomy_analysis/samples_taxonomy_table.txt
   dir=`dirname $file`
   dir=`dirname $dir`
   run=`basename $dir .processed`

   # need to exclude these runs, done after paper accepted (this list may need updating)
   #180305_D00390_0348_BCC33EANXX  180222_D00390_0347_ACC8WAANXX  180208_D00390_0346_BCC5W3ANXX  180208_D00390_0345_ACC5V9ANXX  180205_D00390_0344_ACC5YYANXX  180130_D00390_0343_BCBG7MANXX  180124_D00390_0342_ACBG7VANXX  180122_D00390_0341_BHNVMKBCXY  180119_D00390_0340_AHT555BCXY  180117_D00390_0339_BCBG8CANXX  180111_D00390_0338_ACBG7GANXX  171218_D00390_0337_BCBG3AANXX 
   grep -q $run exclusions.txt
   if [ $? == 0 ]; then
      echo "excluding $run"
      continue
   fi

   if [ $run != $THIS_RUN ]; then
      runs="$runs $run"
   fi
done

# add this run - we want this run's samples to be last in the file
runs="$runs $THIS_RUN"

set -x

for summary in "frequency" "information" ; do
   # summarise all hits
   outfile=$BUILD_ROOT/all_${summary}.txt
   $GBS_BIN/summarise_global_hiseq_taxonomy.py -t $summary  -o $outfile $runs
   outfile=$BUILD_ROOT/all_${summary}_xnohit.txt
   $GBS_BIN/summarise_global_hiseq_taxonomy.py -t $summary  -o $outfile -x nohit $runs  # excluding no hits 

   # make kingdom-specific files as well
   #for kingdom in "eukaryota" "bacteria"; do
   for kingdom in "eukaryota" ; do
      outfile=$BUILD_ROOT/${kingdom}_${summary}.txt
      $GBS_BIN/summarise_global_hiseq_taxonomy.py -s $kingdom -t $summary  -o $outfile $runs
      outfile=$BUILD_ROOT/${kingdom}_${summary}_xnohit.txt
      $GBS_BIN/summarise_global_hiseq_taxonomy.py -s $kingdom -t $summary  -o $outfile -x nohit $runs  # excluding no hits 
   done
done
set +x

# extract the required sample species table (sample_species.txt)
cd $BUILD_ROOT
psql -U agrbrdf -d agrbrdf -h invincible -f $GBS_BIN/database/extract_sample_species.psql
