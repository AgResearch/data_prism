#!/usr/bin/env python 

import os
import re
import itertools
import string
from random import random
from multiprocessing import Pool
import argparse
from prbdf import Distribution , build, bin_discrete_value, get_text_stream , get_file_type,  PROC_POOL_SIZE



def kmer_count_from_sequence(sequence, *args):
    """
    yields an interator through counts of kmers in a sequence
    - e.g.
    3 ACTAT
    1 AAAAA
    etc
    """
    from Bio import SeqIO
    import itertools

    reverse_complement = args[0]
    patterns = args[1:]

    #print "DEBUG sequence%s"%str(sequence)
    #print "DEBUG reverse_complement%s"%str(reverse_complement)
    #print "DEBUG patterns%s"%str(patterns)


    kmer_iters = tuple((re.finditer(pattern, str(sequence.seq), re.I) for pattern in patterns))
    kmer_iters = (match.group() for match in itertools.chain(*kmer_iters))
    if not reverse_complement:
        kmer_count_iter = ( ( len(list(kmer_iter)),kmer) for (kmer,kmer_iter) in itertools.groupby(kmer_iters, lambda kmer:kmer) )
    else:
        kmer_count_iter = ( ( len(list(kmer_iter)),get_reverse_complement(kmer)) for (kmer,kmer_iter) in itertools.groupby(kmer_iters, lambda kmer:kmer) )

    return kmer_count_iter


def seq_from_sequence_file(datafile, *args):
    """
    yields either all or a random sample of seqs from a sequence file
    """
    from Bio import SeqIO
    (filetype, sampling_proportion) = args[0:2]
    seq_iter = SeqIO.parse(get_text_stream(datafile), filetype)

    if sampling_proportion is not None:
        seq_iter = (record for record in seq_iter if random() <= sampling_proportion)
        
        
    return seq_iter


def build_kmer_distribution(datafile, kmer_patterns, sampling_proportion, num_processes, builddir, reverse_complement):

    if os.path.exists(get_save_filename(datafile, builddir)):
        print("build_kmer_distribution- skipping %s as already done"%datafile)
        distob = Distribution.load(get_save_filename(datafile, builddir))
        distob.summary()
    else:
        filetype = get_file_type(datafile)
        distob = Distribution([datafile], num_processes)
        distob.interval_locator_parameters = (None,)
        distob.interval_locator_funcs = (bin_discrete_value,)
        distob.assignments_files = ("kmer_binning.txt",)
        distob.file_to_stream_func = seq_from_sequence_file
        distob.file_to_stream_func_xargs = [filetype,sampling_proportion]
        distob.weight_value_provider_func = kmer_count_from_sequence
        distob.weight_value_provider_func_xargs = [reverse_complement] + kmer_patterns
        
        #distdata = build(distob, use="singlethread")
        distdata = build(distob, proc_pool_size=num_processes)
        distob.save(get_save_filename(datafile, builddir))
            
        print "Distribution %s has %d points distributed over %d intervals, stored in %d parts"%(get_save_filename(datafile, builddir), distob.point_weight, len(distdata), len(distob.part_dict))

    return get_save_filename(datafile, builddir)


def use_kmer_prbdf(picklefile):
    distob = Distribution.load(picklefile)
    distdata = distob.get_distribution()
    for (interval, freq) in distdata.items():
        print interval, freq

def get_save_filename(input_filename, builddir):
    return os.path.join(builddir,"%s.kmerdist.pickle"%(os.path.basename(input_filename)))

def get_patterns(options):
    patterns = options["kmer_regexps"]
    if options["kmer_regexps"] is None:
        patterns = [''.join(ktuple) for ktuple in list(itertools.product(*options["kmer_size"]*['ACGT']))]
    return patterns

def get_reverse_complement(kmer):
    kmer=kmer.upper()
    kmer=kmer.replace("A","t")
    kmer=kmer.replace("T","a")
    kmer=kmer.replace("C","g")
    kmer=kmer.replace("G","c")
    kmer=''.join(reversed(kmer))
    return kmer.upper()
    
    
def build_kmer_distributions(options):
        
    distribution_names = []
    for file_name in options["file_names"]:
        distribution_names.append(build_kmer_distribution(file_name, get_patterns(options), options["sampling_proportion"], \
                                                          options["num_processes"], options["builddir"], options["reverse_complement"]))

    return distribution_names


def summarise_distributions(distributions, options):

    measure = "frequency"
    if options["summary_type"] == "entropy":
        measure = "unsigned_information"

    kmer_intervals = Distribution.get_intervals(distributions, options["num_processes"])

    print "summarising %s , %s across %s"%(measure, str(kmer_intervals), str(distributions))


    sample_measures = Distribution.get_projections(distributions, kmer_intervals, measure, False, options["num_processes"])
    zsample_measures = itertools.izip(*sample_measures)
    sample_name_iter = [tuple([os.path.splitext(os.path.basename(distribution))[0] for distribution in distributions])]
    zsample_measures = itertools.chain(sample_name_iter, zsample_measures)
    interval_name_iter = itertools.chain([("kmer_pattern")],kmer_intervals)
    zsample_measures_with_rownames = itertools.izip(interval_name_iter, zsample_measures)
        
        
    outfile=open(options["output_filename"], "w")

    for interval_measure in zsample_measures_with_rownames:
        print >> outfile, "%s\t%s"%("%s"%interval_measure[0], string.join((str(item) for item in interval_measure[1]),"\t"))
    outfile.close()
    

def get_options():
    description = """
    This script summaries kmer frequencies or entropies for multiple input files. The output is a single tab-delimited text file
    with one row per kmer and one column per input file

    Input files may be fasta or fastq, compressed or uncommpressed. The format of each file is inferred from its suffix.

    A sampling proportion may be specified, in which case a random sample of that proportion of each input file will be taken.

    Multiple processes are started to analyse each file (even if only one file is being processed), with each process doing an interleaved
    read of sequences from the same file. The default number of processes started is 4 (in that case process 1 handles the
    1st, 5th, 9th, etc  sequences in the file; process 2 handles the 2nd, 6th, 10th, etc sequences in the file, etc; 
    results are merged at the end). The -p option can be used to specify more or less processes.

    The summary for each file is saved in the build folder as a seralised python object file. The name of the file is based on the
    name of the input file, with a suffix ".kmerdist.pickle" added. If a serialised summary file alrady exists, the script
    will not bother re-analysing the input file. This means tables can be incrementally built, simply by re-running
    a previous build command, with additional filenames appended.

    """
    long_description = """
examples :

invincible$ ls -l /bifo/archive/Acremonium_et_al/nzgl01804/Raw
total 45209895
-rw-rw-r-- 1 jaureguir fungi_users 3272229466 Oct  5 09:50 C7N0GANXX-1804-01-4-1_L004_R1.fastq.gz
-rw-rw-r-- 1 jaureguir fungi_users 3125568978 Oct  5 09:48 C7N0GANXX-1804-01-4-1_L004_R2.fastq.gz
-rw-rw-r-- 1 jaureguir fungi_users 3266721620 Oct  5 10:12 C7N0GANXX-1804-01-7-1_L004_R1.fastq.gz


# make a table summarising base composition of all files with .fa sufix in current folder
kmer_entropy.py -t frequency -k 1 ./*.fa

# make a table summarising 6-mer entropy for all fastq files in /data/project2
# , based on a random sample of 1/1000 seqs, split over 20 processes
kmer_entropy.py -t entropy -k 6 -p 20 -s .001 /data/project2/*.fastq.gz

# as above , but now include 2 references. If this is run in the same folder as the
# above, the script will re-use the previous results in making the table, so will only
# have to analyse the two new files listed. The two new files will not be randomly sampled
# (no -s option specified), however for the existing files the re-used results are
# based on random sampling. 
kmer_entropy.py -t entropy -k 6 -p 20  /data/project2/*.fastq.gz /references/ref1.fa /references/ref2.fa 
                                                                            
    """
    parser = argparse.ArgumentParser(description=description, epilog=long_description, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file_names', type=str, nargs='+',metavar="filename", help='list of files to summarise')
    parser.add_argument('-t', '--summary_type' , dest='summary_type', default="frequency", choices=["frequency", "entropy"],help="type of summary")
    parser.add_argument('-k', '--kmer_size' , dest='kmer_size', default=6, type=int, help="kmer size (default 6)")
    parser.add_argument('-r', '--kmer_regexp_list' , dest='kmer_regexps', default=None, type=str, help="list of regular expressions (not currently supported)")
    parser.add_argument('-b', '--build_dir' , dest='builddir', default=".", type=str, help="build folder (default '.')")
    parser.add_argument('-p', '--num_processes' , dest='num_processes', default=4, type=int, help="number of processes to start (default 4)")
    parser.add_argument('-s', '--sampling_proportion' , dest='sampling_proportion', default=None, type=float, help="proportion of sequence records to sample (default None means process all records)")
    parser.add_argument('-o', '--output_filename' , dest='output_filename', default="distributions.txt", type=str, help="name of the output file to contain table of kmer distribution summaries for each input file (default 'distributions.txt')")
    parser.add_argument('-c', '--reverse_complement' , dest='reverse_complement', action='store_true', help="for each kmer tabulate the frequency or entropy of its reverse complement (default False)")    

    
    
    args = vars(parser.parse_args())

    # checks
    if args["num_processes"] < 1 or args["num_processes"] > PROC_POOL_SIZE:
        parser.error("num_processes must be between 1 and %d"%PROC_POOL_SIZE)
        
    # either input file or distribution file should exist 
    for file_name in args["file_names"]:
        if not os.path.isfile(file_name) and not os.path.exists(get_save_filename(file_name, args["builddir"])):
            parser.error("couldnot find either %s or %s"%(file_name,get_save_filename(file_name, args["builddir"])))
        break

    # output file should not already exist
    if os.path.exists(args["output_filename"]):
        parser.error("error output file %(output_filename)s already exists"%args)
    
        
    return args


def main():

    options = get_options()
    print options 

    distributions = build_kmer_distributions(options)

    summarise_distributions(distributions, options)   
    
    return 

    
if __name__ == "__main__":
   main()



        

