#!/usr/bin/env python 

import os
import re
import itertools
import string
import exceptions 
from random import random
from multiprocessing import Pool
import subprocess
import argparse
from data_prism import prism , build, bin_discrete_value, get_text_stream , get_file_type,  PROC_POOL_SIZE


class kmer_prism_exception(exceptions.Exception):
    def __init__(self,args=None):
        super(kmer_prism_exception, self).__init__(args)

#********************************************************************
# methods for getting kmer counts from sequence files
#********************************************************************

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
    pattern_window_length = args[1]  # optional - for fixed length patterns e.g. 6-mers etc, to speed up search
    if callable(args[2]):
        weight=args[2](sequence)
    else:
        weight=args[2]
       
    patterns = args[3:]

    #print "DEBUG sequence%s"%str(sequence)
    #print "DEBUG reverse_complement%s"%str(reverse_complement)
    #print "DEBUG patterns%s"%str(patterns)

    if pattern_window_length is None:
        # search for each pattern. Note that this does not count multiple instances 
        # of a pattern that overlap - for example in TTTTTTT , the pattern TTTTTT will only count once. 
        kmer_iters = tuple((re.finditer(pattern, str(sequence.seq), re.I) for pattern in patterns))
        kmer_iters = (match.group() for match in itertools.chain(*kmer_iters))
        if not reverse_complement:
            kmer_count_iter = ( ( weight * len(list(kmer_iter)),kmer) for (kmer,kmer_iter) in itertools.groupby(kmer_iters, lambda kmer:kmer) )
        else:
            kmer_count_iter = ( ( weight * len(list(kmer_iter)),get_reverse_complement(kmer)) for (kmer,kmer_iter) in itertools.groupby(kmer_iters, lambda kmer:kmer) )
    else:
        # slide the window along the sequence and accumulate matching patterns.Note that unlike
        # the above regexp based search, this would count multiple instances of a pattern
        # that overlap - for example in TTTTTTT , the pattern TTTTTT would count twice.
        # overlap_patterns is used to emulate the regexp behaviour 
        strseq = str(sequence.seq)
        kmer_dict = {}
        overlap_patterns = pattern_window_length * [""]        
        kmer_iter = (strseq[i:i+pattern_window_length] for i in range(0,1+len(strseq)-pattern_window_length))
        for kmer in kmer_iter:
            if kmer not in overlap_patterns:
                overlap_patterns.insert(0,kmer)
            elif overlap_patterns[-1] == kmer:
                overlap_patterns.insert(0,kmer)
            else:
                overlap_patterns.insert(0,"")
            overlap_patterns.pop()

            if kmer not in overlap_patterns[1:]:
                kmer_dict[kmer] = weight + kmer_dict.setdefault(kmer,0)
                
        kmer_count_iter = ( (kmer_dict[kmer], kmer) for kmer in kmer_dict )
        
        
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

def parse_weight_from_sequence_description(sequence):
    # this is used where the fasta file is marked up with a count , and this ends up in the 
    # description of the seq_record - e.g. where the sequence is from a gbs tag count file
    # assumed the count is encoded like this in the description
    # seq_28639 count=1.004008
    weighting_match = re.search("count=(\d*\.*\d*)\s*$", sequence.description)
    if weighting_match is not None:
        weight = float(weighting_match.groups()[0])
    else:
        weight = 1.0

    return weight


#********************************************************************
# methods for getting kmer counts from tag count files 
#********************************************************************
def tag_count_from_tag_count_file(datafile, *args):
    """
    yields either all or a random sample of tag counts  from a tassel tag count file.
    This method trims each tag sequence , to remove the poly-A padding added by tassel and
    optionally any common prefix
    It returns an iterator over tuples like (tag , count).

    The tag count file contains records like
TGCAGAAGTCTTGAATTTAATTCAGGATACTCGTCTACCACGTTGTCCATGTCTCCGCAAGGGA        64      1
TGCAGAAGTCTTGAATTTAGTTCAGGATACTCGTCTACCACGTTGTCCATGTCTCCGCAAAGGA        64      1
TGCAGAAGTCTTGGCCTGAGGAGCTGAGTTGTGCATCACCCTGCAAAAAAAAAAAAAAAAAAAA        45      3
TGCAGAAGTCTTGGTGATGTTGTAAAGGTGTGTTGATGTCTCTGTGGTTGAGGACACATCATCA        64      3

- the first number indicates how much of the tag to keep , the second number
indicates how many of that tag there are

    """
    (input_driver_config, sampling_proportion) = args[0:2]

    if input_driver_config is None:
        raise kmer_prism_exception("""
must supply input driver config for .cnt files. This consists of the name of a
script which lists the contents of the tile in text. Example code using tassel3 and bash shell:
mkfifo f1
nohup run_pipeline.pl -fork1 -BinaryToTextPlugin  -i $infile -o f1 -t TagCounts -endPlugin -runfork1 >$errfile 2>&1 &
cat <$f1
""")
    else:
        remove_prefix=True  # hard coded true for now but may pass in as part of drive config at some point
        common_prefix=""
        cat_tag_count_command = [input_driver_config, "%s"%datafile]


        # if we are to remove a common prefix (e.g. TGCA in the above example), then
        # scan the tags to identify the prefix
        if remove_prefix:
            print "scanning tags for a common prefix to remove..."

            proc = subprocess.Popen(cat_tag_count_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (stdout, stderr) = proc.communicate()
            if proc.returncode == 0:
                tuple_iter = (re.split("\s+",record.strip().upper()) for record in re.split("\n", stdout))    # parse the 3 elements 
                tuple_iter = ((my_tuple[0], int(my_tuple[1]), int(my_tuple[2]))  for my_tuple in tuple_iter if len(my_tuple) == 3)  # skip the header and make ints
                tuple_iter = (my_tuple[0][0:my_tuple[1]] for my_tuple in tuple_iter)  # use the tag-length to substring the tag then throw away the numbers
                sorted_tuples = sorted(tuple_iter)
                # find the longest common start-string in the first and last elements
                common_prefix_length = 0
                match = True
                while( common_prefix_length < min(len(sorted_tuples[0]), len(sorted_tuples[-1])) ):
                    # try length + 1
                    if sorted_tuples[0][0:common_prefix_length+1] == sorted_tuples[-1][0:common_prefix_length+1]:
                        common_prefix_length += 1
                    else:
                        break
            else:
                raise kmer_prism_exception("Error encountered running %s - return code was %s, stderr:%s stdout:%s"%(" ".join(cat_tag_count_command),proc.returncode,stderr,stdout))

            if common_prefix_length > 0 :
                common_prefix = sorted_tuples[0][0:common_prefix_length]
                print "found common prefix %s - will exclude from analysis"%common_prefix
            else:
                print "(no common prefix found)"

        print "summarising tags..."
        proc = subprocess.Popen(cat_tag_count_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = proc.communicate()
        if proc.returncode == 0:
            tagcount_iter = (re.split("\s+",record.strip().upper()) for record in re.split("\n", stdout))    # parse the 3 elements 
            tagcount_iter = ((my_tuple[0], int(my_tuple[1]), int(my_tuple[2]))  for my_tuple in tagcount_iter if len(my_tuple) == 3)  # skip the header and make ints
            tagcount_iter = ( (my_tuple[0][common_prefix_length:my_tuple[1]], my_tuple[2]) for my_tuple in tagcount_iter)  # use the tag-length to substring the tag, also remove common prefix, and throw away tag length
            #for record in tagcount_iter:
            #    print "DEBUG %s"%str(record)
            #    sys.exit(1)
        else:
            raise kmer_prism_exception("Error encountered running %s"%" ".join(cat_tag_count_command))

    #print "DEBUG got tag count iter"
    return tagcount_iter

def kmer_count_from_tag_count(tag_count_tuple, *args):
    """
    yields an interator through counts of kmers in a tag - but multiplied
    up by the tag count.
    """
    reverse_complement = args[0]
    pattern_window_length = args[1]  # optional - for fixed length patterns e.g. 6-mers etc, to speed up search
    weight = args[2]  # un-used currently
    patterns = args[3:]
    (tag,tag_count) = tag_count_tuple

    #print "DEBUG args, tag count%s"%str(args, tag_count_tuple)

    if pattern_window_length is None:
        # search for each pattern. Note that this does not count multiple instances 
        # of a pattern that overlap - for example in TTTTTTT , the pattern TTTTTT will only count once. 
        kmer_iters = tuple((re.finditer(pattern, tag, re.I) for pattern in patterns))
        kmer_iters = (match.group() for match in itertools.chain(*kmer_iters))
        if not reverse_complement:
            kmer_count_iter = ( ( len(list(kmer_iter)),kmer) for (kmer,kmer_iter) in itertools.groupby(kmer_iters, lambda kmer:kmer) )
        else:
            kmer_count_iter = ( ( len(list(kmer_iter)),get_reverse_complement(kmer)) for (kmer,kmer_iter) in itertools.groupby(kmer_iters, lambda kmer:kmer) )
    else:
        # slide the window along the sequence and accumulate matching patterns.Note that unlike
        # the above regexp based search, this would count multiple instances of a pattern
        # that overlap - for example in TTTTTTT , the pattern TTTTTT would count twice.
        # overlap_patterns is used to emulate the regexp behaviour 
        strseq = tag
        kmer_dict = {}
        overlap_patterns = pattern_window_length * [""]        
        kmer_iter = (strseq[i:i+pattern_window_length] for i in range(0,1+len(strseq)-pattern_window_length))
        for kmer in kmer_iter:
            if kmer not in overlap_patterns:
                overlap_patterns.insert(0,kmer)
            elif overlap_patterns[-1] == kmer:
                overlap_patterns.insert(0,kmer)
            else:
                overlap_patterns.insert(0,"")
            overlap_patterns.pop()

            if kmer not in overlap_patterns[1:]:
                kmer_dict[kmer] = 1 + kmer_dict.setdefault(kmer,0)
                
        kmer_count_iter = ( (tag_count * kmer_dict[kmer], kmer) for kmer in kmer_dict )
        
        
    return kmer_count_iter

#********************************************************************
# general analysis / summary methods 
#********************************************************************
def build_kmer_spectrum(datafile, kmer_patterns, sampling_proportion, num_processes, builddir, reverse_complement, pattern_window_length, input_driver_config, input_filetype=None, weighting_method = None):

    if os.path.exists(get_save_filename(datafile, builddir)):
        print("build_kmer_spectrum- skipping %s as already done"%datafile)
        kmer_prism = prism.load(get_save_filename(datafile, builddir))
        kmer_prism.summary()
        
    else:
        print("build_kmer_spectrum- processing %s"%datafile)
        filetype = input_filetype
        if filetype is None:
            filetype = get_file_type(datafile)
        kmer_prism = prism([datafile], num_processes)
        kmer_prism.interval_locator_parameters = (None,)
        kmer_prism.interval_locator_funcs = (bin_discrete_value,)
        kmer_prism.assignments_files = ("kmer_binning.txt",)
        kmer_prism.file_to_stream_func = seq_from_sequence_file
        kmer_prism.file_to_stream_func_xargs = [filetype,sampling_proportion]
        kmer_prism.spectrum_value_provider_func = kmer_count_from_sequence

        if weighting_method is None:
            kmer_prism.spectrum_value_provider_func_xargs = [reverse_complement, pattern_window_length, 1] + kmer_patterns        
        elif weighting_method == "tag_count":
            kmer_prism.spectrum_value_provider_func_xargs = [reverse_complement, pattern_window_length, parse_weight_from_sequence_description] + kmer_patterns
        
        if filetype == ".cnt":
            #print "DEBUG setting methods for count file"
            kmer_prism.file_to_stream_func = tag_count_from_tag_count_file
            kmer_prism.file_to_stream_func_xargs = [input_driver_config,sampling_proportion]
            kmer_prism.spectrum_value_provider_func = kmer_count_from_tag_count 
            spectrum_data = build(kmer_prism, use="singlethread")
        else:
            spectrum_data = build(kmer_prism, proc_pool_size=num_processes)
            
        kmer_prism.save(get_save_filename(datafile, builddir))
            
        print "spectrum %s has %d points distributed over %d intervals, stored in %d parts"%(get_save_filename(datafile, builddir), kmer_prism.total_spectrum_value, len(spectrum_data), len(kmer_prism.part_dict))

    return get_save_filename(datafile, builddir)


def assemble_kmer_spectrum(kmer_list_file, sequence_file, sequence_file_type, sampling_proportion, input_driver_config = None,counts_file = None):

    # get the kmer list
    with open(kmer_list_file , "r") as kmer_stream:
        kmer_list = [ item.strip() for item in kmer_stream if len(item.strip()) > 0 ]

    # get an iter of (sequence, count)
    if sequence_file_type is None:
        filetype = get_file_type(sequence_file)
    file_to_stream_func = seq_from_sequence_file
    file_to_stream_func_xargs = [filetype,sampling_proportion]
    if filetype == ".cnt":
        file_to_stream_func = tag_count_from_tag_count_file
        file_to_stream_func_xargs = [input_driver_config,sampling_proportion]
        zsequences_counts_stream=file_to_stream_func(sequence_file, *file_to_stream_func_xargs)
    else:
        if counts_file is None:
            zsequences_counts_stream = itertools.izip(file_to_stream_func(sequence_file, *file_to_stream_func_xargs), itertools.repeat(1))
        else:
            with open(counts_file,"r") as counts_stream:
                counts_iter = (int(item.strip()) for item in counts_stream if len(item.strip()) > 0)
                zsequences_counts_stream = itertools.izip(file_to_stream_func(sequence_file, *file_to_stream_func_xargs), counts_iter)
                
    #for record in zsequences_counts_stream:
    #    print record
    #
    #return
        
    unassembled_dict = {} # key is a tuple of supporting kmers, value is count 

    pattern_window_length = max( len(kmer) for kmer in kmer_list)
    if pattern_window_length != min( len(kmer) for kmer in kmer_list):
        raise trim_exception("error -  all kmers in supporting list mustbe the same length")
    
    for (sequence, sequence_count) in zsequences_counts_stream:      # z prefix denotes a zipped stream  - returns tuple of (seq, count)
        # analyse sequence     
        # slide the window along the sequence and accumulate exact matches to members of patterns.
        #print sequence, sequence_count
        strseq = str(sequence.seq)
        supporting_run = []
        supporting_runs = []
        kmer_iter = (strseq[i:i+pattern_window_length] for i in range(0,1+len(strseq)-pattern_window_length))
        #print "assembling %s using %s"%(str(kmer_list), strseq)
        for kmer in kmer_iter:
            #print "checking %s"%kmer
            if kmer in kmer_list:
                #print "added"
                supporting_run.append(kmer) 
            else:
                if len(supporting_run) > 0:   # conclude this supporting run and append to the list of supporting runs
                    supporting_runs.append(supporting_run)
                    supporting_run = []

        # check for a supporting run that included the last kmer
        if len(supporting_run) > 0:
           supporting_runs.append(supporting_run)

        # if there are any supporting runs, store the longest
        if len(supporting_runs) > 0:
            supporting_runs = sorted(supporting_runs, lambda a,b:cmp(len(a), len(b)), None, True)
            best_supporting_run = tuple(supporting_runs[0])
            unassembled_dict[best_supporting_run] = sequence_count + unassembled_dict.setdefault(best_supporting_run,0)

    # summarise the length distribution of supporting kmer runs 
    unassembled_dist = {}
    for (kmer_tuple, count) in unassembled_dict.items():
        unassembled_dist[len(kmer_tuple)] = count + unassembled_dist.setdefault(len(kmer_tuple),0)

    

    # assemble each of the runs 
    assembled_dict = {}
    for (kmer_tuple, count)in unassembled_dict.items():
        if len(kmer_tuple) == 0:
            assembly = ""
        elif len(kmer_tuple) == 1:
            assembly = kmer_tuple[0]
        else:
            assembly = reduce(lambda x,y:x+y[-1],kmer_tuple[1:],kmer_tuple[0])
        assembled_dict[assembly] = count
    
    print "\n\n\n"
    for key in sorted(unassembled_dist.keys()):
        print "%s\t%s"%(key, unassembled_dist[key])

    print "\n\n\n"        
    for key in sorted(assembled_dict.keys(),lambda x,y:cmp(len(x),len(y)),None,True):
        print "%s\t%s"%(key, assembled_dict[key])

def use_kmer_prbdf(picklefile):
    kmer_prism = prism.load(picklefile)
    spectrum_data = kmer_prism.get_distribution()
    for (interval, freq) in spectrum_data.items():
        print interval, freq

def get_save_filename(input_filename, builddir):
    sanitised_input_filename = re.sub("[\s\$]","_", input_filename)
    return os.path.join(builddir,"%s.kmerdist.pickle"%(os.path.basename(sanitised_input_filename)))


def get_reverse_complement(kmer):
    kmer=kmer.upper()
    kmer=kmer.replace("A","t")
    kmer=kmer.replace("T","a")
    kmer=kmer.replace("C","g")
    kmer=kmer.replace("G","c")
    kmer=''.join(reversed(kmer))
    return kmer.upper()
    
    
def build_kmer_spectra(options):
        
    spectrum_names = []
    for file_name in options["file_names"]:
        spectrum_names.append(build_kmer_spectrum(file_name, options["kmer_regexps"], options["sampling_proportion"], \
                                                          options["num_processes"], options["builddir"], options["reverse_complement"], \
                                                          options["kmer_size"], options["input_driver_config"], options["input_filetype"], options["weighting_method"]))

    return spectrum_names


def summarise_spectra(distributions, options):

    measure = "raw"
    if options["summary_type"] in ["zipfian","entropy"]:
        measure = "unsigned_information"

    kmer_intervals = prism.get_intervals(distributions, options["num_processes"])

    if options["alphabet"] is not None:
        kmer_intervals1 = [ interval for interval in kmer_intervals if re.search("^[%(alphabet)s]+$"%options , interval[0], re.IGNORECASE) is not None ]
        print "(restricting kmers to those from alphabet %s , deleted %d / %d kmers)"%(options["alphabet"],len(kmer_intervals) - len(kmer_intervals1) , len(kmer_intervals))
        kmer_intervals  = kmer_intervals1
        

    #print "summarising %s , %s across %s"%(measure, str(kmer_intervals), str(distributions))
    print "summarising %s , %d kmers across %s"%(measure, len(kmer_intervals), str(distributions))


    sample_measures = prism.get_projections(distributions, kmer_intervals, measure, False, options["num_processes"])
    zsample_measures = itertools.izip(*sample_measures)
    sample_name_iter = [tuple([os.path.splitext(os.path.basename(distribution))[0] for distribution in distributions])]
    zsample_measures = itertools.chain(sample_name_iter, zsample_measures)
    interval_name_iter = itertools.chain([("kmer_pattern")],kmer_intervals)
    
    outfile=open(options["output_filename"], "w")

    if options["summary_type"] in ["entropy", "frequency"]:
        zsample_measures_with_rownames = itertools.izip(interval_name_iter, zsample_measures)
        for interval_measure in zsample_measures_with_rownames:
            print >> outfile, "%s\t%s"%("%s"%interval_measure[0], string.join((str(item) for item in interval_measure[1]),"\t"))
        outfile.close()
    elif options["summary_type"] in ["ranks", "zipfian"]:
        # duplicate interval_name_iter - needed 3 times
        interval_name_iter_dup = itertools.tee(interval_name_iter, 3)

        # triplicate zsample_measures (0 used to get ranks; 1 used to output measures; 3 used to get distances)
        zsample_measures_dup = itertools.tee(zsample_measures,3)
        ranks = prism.get_rank_iter(zsample_measures_dup[0])

        # duplicate ranks (0 used to output; 1 used to get distances)
        ranks_dup = itertools.tee(ranks, 2)
        ranks_with_rownames = itertools.izip(interval_name_iter_dup[0], ranks_dup[0])

        # output ranks
        print >> outfile , "*** ranks *** :"
        for interval_rank in ranks_with_rownames:
            print >> outfile, "%s\t%s"%("%s"%interval_rank[0], string.join((str(item) for item in interval_rank[1]),"\t"))

        # output measures
        print >> outfile , "*** entropies *** :"
        zsample_measures_with_rownames = itertools.izip(interval_name_iter_dup[1], zsample_measures_dup[1])
        for interval_measure in zsample_measures_with_rownames:
            print >> outfile, "%s\t%s"%("%s"%interval_measure[0], string.join((str(item) for item in interval_measure[1]),"\t"))

        # get distances
        print >> outfile , "*** distances *** :"
        (distance_matrix, point_names_sorted) = prism.get_zipfian_distance_matrix(zsample_measures_dup[2], ranks_dup[1])
        prism.print_distance_matrix(distance_matrix, point_names_sorted, outfile)
    else:
        print "warning, unknown summary type %(summary_type)s, no summary available"%options
        
        
        outfile.close()


def get_options():
    description = """
    This script summaries kmer frequencies or entropies for multiple input files. The output is a single tab-delimited text file
    with one row per kmer and one column per input file

    Input files may be fasta or fastq, compressed or uncompressed. The format of each file is inferred from its suffix.

    A sampling proportion may be specified, in which case a random sample of that proportion of each input file will be taken.

    Multiple processes are started to analyse each file (even if only one file is being processed), with each process doing an interleaved
    read of sequences from the same file. The default number of processes started is 4 (in that case process 1 handles the
    1st, 5th, 9th, etc  sequences in the file; process 2 handles the 2nd, 6th, 10th, etc sequences in the file, etc; 
    results are merged at the end). The -p option can be used to specify more or less processes.

    The kmer summary for each input file is cached in the build folder as a serialised python object file. The name of the file is based on the
    name of the input file, with a suffix ".kmerdist.pickle" added. If a serialised summary is already cached, the script
    will not bother re-analysing the input file. This means the all-files summary table can be incrementally built, simply by re-running
    a previous build command, with additional filenames appended.

    """
    long_description = """
examples :


# make a table summarising base composition of all files with .fa suffix in current folder
kmer_prism.py -t frequency -k 1 ./*.fa

# make a table summarising 6-mer self-information for all fastq files in /data/project2
# , based on a random sample of 1/1000 seqs, split over 20 processes
kmer_prism.py -t entropy -k 6 -p 20 -s .001 /data/project2/*.fastq.gz

# as above , but now also include 2 reference genomes. If this is run in the same folder as the
# above, the script will re-use the previously cached results, so will only
# have to analyse the kmer distribution for the two new files listed. The       if filetype is None:
            filetype = get_file_type(datafile)two new files will not
# be randomly sampled (no -s option specified), however for the existing files the cached results are
# based on a random sample.  
kmer_prism.py -t entropy -k 6 -p 20  /data/project2/*.fastq.gz /references/ref1.fa /references/ref2.fa

# obtain a text file containing self-information and ranks for 6-mers in a tag count file
./kmer_prism.py -t zipfian -k 6 -p 1 -o tag_zipfian.txt -x /dataset/hiseq/active/bin/hiseq_pipeline/cat_tag_count.sh /dataset/hiseq/scratch/postprocessing/151016_D00390_0236_AC6JURANXX.gbs/SQ0124.processed_sample/uneak/tagCounts/G88687_C6JURANXX_1_124_X4.cnt

# as above but feeed in tag count data from a text file, use "cat" to list it (but need .cnt suffix so
# program expects the stream to contain tag counts
./kmer_prism.py -t frequency -k 6 -p 1 -o test_freqs.txt -x cat tagtestdata_as_text.cnt



                                                                            
    """
    parser = argparse.ArgumentParser(description=description, epilog=long_description, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file_names', type=str, nargs='+',metavar="filename", help='list of files to process')
    parser.add_argument('-t', '--summary_type' , dest='summary_type', default="frequency", choices=["frequency", "entropy", "ranks", "zipfian", "assembly", "test"],help="type of summary")
    parser.add_argument('-k', '--kmer_size' , dest='kmer_size', default=None, type=int, help="kmer size (default None)")
    parser.add_argument('-r', '--kmer_regexp_list' , dest='kmer_regexps', default=None, type=str, help="list of regular expressions (not currently supported)")
    parser.add_argument('-b', '--build_dir' , dest='builddir', default=".", type=str, help="build folder (default '.')")
    parser.add_argument('-p', '--num_processes' , dest='num_processes', default=4, type=int, help="number of processes to start (default 4)")
    parser.add_argument('-s', '--sampling_proportion' , dest='sampling_proportion', default=None, type=float, help="proportion of sequence records to sample (default None means process all records)")
    parser.add_argument('-o', '--output_filename' , dest='output_filename', default="distributions.txt", type=str, help="name of the output file to contain table of kmer distribution summaries for each input file (default 'distributions.txt')")
    parser.add_argument('-c', '--reverse_complement' , dest='reverse_complement', action='store_true', help="for each kmer tabulate the frequency or entropy of its reverse complement (default False)")
    parser.add_argument('-x', '--input_driver_config' , dest='input_driver_config', default=None, help="this is use to configure input from custom file formats such as tassel count files")    
    parser.add_argument('-a', '--alphabet' , dest='alphabet', default=None, type=str, help="alphabet used to filter kmers when summarising distributions (not applied when building distribution)")
    parser.add_argument('-f', '--input_filetype' , dest='input_filetype', default=None, type=str, choices=["fasta", "fastq"], help="optionally specify input format ( if not specified will try to guess)")
    parser.add_argument('--kmer_listfile' , dest='kmer_listfile', default=None, type=str,  help="list of kmers for an assembly run")
    parser.add_argument('--sequence_countfile' , dest='sequence_countfile', default=None, type=str,  help="sequence count file")
    parser.add_argument('--weighting_method' , dest='weighting_method', default=None, type=str,  choices=["tag_count"], help="weighting method")
    
    
    args = vars(parser.parse_args())


    # checks
    if args["summary_type"] != "assembly": 
        if args["num_processes"] < 1 or args["num_processes"] > PROC_POOL_SIZE:
            parser.error("num_processes must be between 1 and %d"%PROC_POOL_SIZE)

        # should specify either a kmer_size, or a list of patterns (but not both)

        if args["kmer_size"] is None and args["kmer_regexps"] is None:
            parser.error("should specify either kmer_size or a list of patterns")
        elif args["kmer_size"] is not None and args["kmer_regexps"] is not None:
            parser.error("should specify either kmer_size or a list of patterns but not both")
            
        # either input file or distribution file should exist 
        for file_name in args["file_names"]:
            if not os.path.isfile(file_name) and not os.path.exists(get_save_filename(file_name, args["builddir"])):
                parser.error("could not find either %s or %s"%(file_name,get_save_filename(file_name, args["builddir"])))
            break

        # output file should not already exist
        if os.path.exists(args["output_filename"]):
            parser.error("error output file %(output_filename)s already exists"%args)


        # parse kmer_regexps
        if args["kmer_regexps"] is not None:
            args["kmer_regexps"] = re.split("\s*,\s*", args["kmer_regexps"])
        else:
            args["kmer_regexps"]= []
    
        
    return args

def test(options):
    from Bio import SeqIO
    for file_name in options["file_names"]:
        filetype = get_file_type(file_name)
        seq_iter = SeqIO.parse(get_text_stream(file_name), filetype)
        for seq in seq_iter:
            print seq.description
        


def main():

    options = get_options()
    print options

    if options["summary_type"] == "test":
        test(options)
        

    if options["summary_type"] != "assembly":
        distributions = build_kmer_spectra(options)
        summarise_spectra(distributions, options)   
    else:
        assemble_kmer_spectrum(options["kmer_listfile"], options["file_names"][0], options["input_filetype"], options["sampling_proportion"], \
                               options["input_driver_config"],options["sequence_countfile"])

    return 

    
if __name__ == "__main__":
   main()



        

