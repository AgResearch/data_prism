#!/usr/bin/env python
from __future__ import print_function

import itertools,os,re,argparse,string,sys
sys.path.append('/dataset/bioinformatics_dev/active/data_prism')
from data_prism import prism , build, bin_discrete_value, get_text_stream , get_file_type,  PROC_POOL_SIZE

def my_first_hit_provider(filename, *xargs):
    """
    provide first-hits from a standard blastn output , like
Query= NC_002944.2_1-100 Mycobacterium avium subsp. paratuberculosis str.
k10, complete genome

Length=100
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

  NC_002944.2 Mycobacterium avium subsp. paratuberculosis str. k1...  185     5e-48


> NC_002944.2 Mycobacterium avium subsp. paratuberculosis str.
k10, complete genome
Length=4829781

 Score = 185 bits (100),  Expect = 5e-48
 Identities = 100/100 (100%), Gaps = 0/100 (0%)
 Strand=Plus/Plus

Query  1   TTGGCCGATGACCCCGGTTCAAGCTTCACCACGGTGTGGAATGCGGTCGTTTCGGAGCTC  60
           ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1   TTGGCCGATGACCCCGGTTCAAGCTTCACCACGGTGTGGAATGCGGTCGTTTCGGAGCTC  60

Query  61  AACGGCGAGCCCGTCGCCGACGGCGGAGCCGCCAACCGCA  100
           ||||||||||||||||||||||||||||||||||||||||
Sbjct  61  AACGGCGAGCCCGTCGCCGACGGCGGAGCCGCCAACCGCA  100
.
.
> another_hit
etc

Th context of this is in summarising blast results against a single genome, where we are
interested in whether the query sequence has a good match in that genome or not , so that we are only interested
in the first match. Whether the sequence also has a hit to other genomes is determined by a series
of other blasts runs. 
"""
    pattern_to_match = xargs[0]
    pattern_to_not_match = xargs[1]
        
    tuple_stream = get_text_stream(filename)
    tuple_stream=(re.split("\s+", record.strip()) for record in tuple_stream)
    atuple = tuple_stream.next()

    query = None
    identity = None
    text = None
    while True:
        if len(atuple) <= 1:
            try:
                atuple = tuple_stream.next()
            except StopIteration:
                if query is not None and identity is not None:
                    yield tuple([query, identity])
                raise StopIteration                     
            continue
        if atuple[0] == "Query=":
            query=atuple[1]
            identity = None
        if atuple[0][0] == ">":
            text = (" ".join(atuple)).lower()
        elif atuple[0] == "Identities":
            if query is not None:
                identity = atuple[3]
                if pattern_to_match is not None:
                    match = re.search(pattern_to_match, text, re.IGNORECASE)
                if pattern_to_not_match is not None:
                    not_match = re.search(pattern_to_not_match, text, re.IGNORECASE)
                    

                if pattern_to_match is None and pattern_to_not_match is None:
                    yield tuple([query, identity])   # i.e. yield the first identity encountered
                                                     # (unless there is a pattern that was not matched by hit description)
                elif pattern_to_match is not None and match is not None:
                    yield tuple([query, identity])   # i.e. yield the first identity encountered
                                                     # (unless there is a pattern that was not matched by hit description)
                elif pattern_to_not_match is not None and not_match is None:
                    yield tuple([query, identity])   # i.e. yield the first identity encountered
                                                     # (unless there is a pattern that was not matched by hit description)
                
                query = None
                identity = None
                text = None
        try:
            atuple = tuple_stream.next()
        except StopIteration:
            if query is not None and identity is not None:
                yield tuple([query, identity])
            raise StopIteration



def my_good_hit_provider(filename, *xargs):
    """
    provide good-hits from a standard blastn output , like
Query= NC_002944.2_1-100 Mycobacterium avium subsp. paratuberculosis str.
k10, complete genome

Length=100
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

  NC_002944.2 Mycobacterium avium subsp. paratuberculosis str. k1...  185     5e-48


> NC_002944.2 Mycobacterium avium subsp. paratuberculosis str.
k10, complete genome
Length=4829781

 Score = 185 bits (100),  Expect = 5e-48
 Identities = 100/100 (100%), Gaps = 0/100 (0%)
 Strand=Plus/Plus

Query  1   TTGGCCGATGACCCCGGTTCAAGCTTCACCACGGTGTGGAATGCGGTCGTTTCGGAGCTC  60
           ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1   TTGGCCGATGACCCCGGTTCAAGCTTCACCACGGTGTGGAATGCGGTCGTTTCGGAGCTC  60

Query  61  AACGGCGAGCCCGTCGCCGACGGCGGAGCCGCCAACCGCA  100
           ||||||||||||||||||||||||||||||||||||||||
Sbjct  61  AACGGCGAGCCCGTCGCCGACGGCGGAGCCGCCAACCGCA  100
.
.
> another_hit
etc

The context of this is in summarising blast results against genbank , where we are
interested in summarising all "good hits" that the query has. We are only interested
in counting the taxa that are hit, so we just yield a text description 
"""
    tuple_stream = get_text_stream(filename)
    tuple_stream=(re.split("\s+", record.strip()) for record in tuple_stream)

    while True:
        atuple = tuple_stream.next()
        if len(atuple) <= 1:
            continue 
        if atuple[0][0] == ">":
            text = (" ".join(atuple[1:5]))
        if atuple[0] == "Identities":
            # only consider alignments across the length of the query , i.e 100, and with at least 95 identities 
            match=re.search("(\d+)/(\d+)", atuple[2])
            (identities,length) = (int(match.groups()[0]), int(match.groups()[1]))
            if length == 100 and identities >= 95:
                yield (text,)
                    

def my_top_hit_provider(filename, *xargs):
    """
    takes a stream which may contain multiple hits, and yields just the top hit in each group
    (in this case this is unnecessary as my_first_hit_provider only yields a single hit)
    """
    groups = itertools.groupby(my_first_hit_provider(filename, *xargs), lambda x:x[0])    
    top_hits = (group.next() for (key, group) in groups)
    return top_hits

def pctidentity_value_provider(point, *xargs):
    # top hit provider returns query, statistic  - weight value provider returns statistic, queryid
    (query, identity_string ) = point   # e.g. NC_002944.2_1-100   ,(100%)
    value=identity_string.translate(None, ",()%") 
    value_tuple = (value, query)    
    return (value_tuple,)


def test(datafile):
    #s = my_first_hit_provider(datafile, 0)
    #for t in s:
    #    print(t)
    #s=my_top_hit_provider(datafile,0)
    s=my_good_hit_provider(datafile,0)
    
    for t in s:
        print(t)


def build_spectrum(datafile, options, target="first_hit"):
    hit_prism = prism([datafile], 1)
    #distob.DEBUG=True
    if target=="first_hit":
        hit_prism.file_to_stream_func = my_top_hit_provider
        hit_prism.file_to_stream_func_xargs = [options["hit_pattern_to_match"], options["hit_pattern_to_not_match"]]
        hit_prism.spectrum_value_provider_func = pctidentity_value_provider
        hit_prism.interval_locator_funcs = [bin_discrete_value]        
    else:
        hit_prism.file_to_stream_func = my_good_hit_provider
        hit_prism.interval_locator_funcs = [bin_discrete_value]        
    
    spectrum_data = build(hit_prism,"singlethread")
    hit_prism.save("%s.pickle"%datafile)
    return spectrum_data

def list_spectrum(datafile, options):
    hit_prism = prism.load("%s.pickle"%datafile)
    spectrum_data = hit_prism.get_spectrum()
    for (interval, freq) in spectrum_data.items():
        print(interval, freq)
    

def summarise_spectra(options):

    measure=options["summary_type"]
        
    spectrum_names = ["%s.pickle"%datafile for datafile in options["filename"]]

    print("getting intervals...")
    probe_intervals = prism.get_intervals(spectrum_names, options["num_processes"])

    #print "summarising %s , %s across %s"%(measure, str(kmer_intervals), str(distributions))
    print("summarising %s , %d probes across %s"%(measure, len(probe_intervals), str(spectrum_names)))

    print("getting projections...")
    probe_measures = prism.get_projections(spectrum_names, probe_intervals, measure, False, options["num_processes"])

    print("output projections...")
    if options["summary_type"] in ["raw"]:
        prism.save_projections(spectrum_names, probe_intervals, probe_measures, options["output_filename"])
    else:
        print("warning, unknown summary type %(summary_type)s, no summary available"%options)

def query_spectra(options):
    target_spectra = prism.load_projections(options["filename"][0])
    query_spectra = prism.load_projections(options["filename"][1])
    target_names = prism.load_projections(options["filename"][0], colnames_only = True)    
    query_results = prism.query_spectra(target_spectra, query_spectra, target_names)

    # extract and sort those results with a distance of less than 0.001 from the query
    # query_results is a list of tuples , (probe_name, distance)
    # probe names are lie this : NC_002945.4_999201-999300 ,i.e. name_start_stop - we want
    # to sort by accession and position (and then later , reassemble contiguous fragments)
    close_results  = [ item[0] for item in query_results if item[1] < .001 ]

    # define a custom comparator , which knows about both accession name and start position
    def my_sorter(result1, result2):
        tokens1 = re.split("_", result1)  # e.g. splitting maybe NC_002945.4_999201-999300
        tokens2 = re.split("_", result2)

        
        (start1, stop1)  = re.split("-", tokens1[-1])
        (start2, stop2)  = re.split("-", tokens2[-1])
        
        name1 = "_".join(tokens1[0:len(tokens1)-1])
        name2 = "_".join(tokens2[0:len(tokens2)-1])
                         

        if name1 != name2:
            return cmp(name1,name2)
        else:
            return cmp(int(start1), int(start2))

    sorted_close_results = sorted(close_results, cmp=my_sorter)
    for record in sorted_close_results:
        print("%s"%record)
        
    
def get_options():
    description = """
    """
    long_description = """
examples :

summarise_probe_blasts.py -a "build_good" design/Probe7.fa.genbankhits.gz


    """

    parser = argparse.ArgumentParser(description=description, epilog=long_description, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', type=str, nargs="*",help='input file of blast hits (optionally compressed with gzip)')    
    parser.add_argument('--use', dest='use', default="use", choices=["pctidentity"],help="use")
    parser.add_argument('-a', dest='action', default="build_first", choices=["build_first", "build_good", "list", "summarise", "query"],help="action")
    parser.add_argument('-p', dest='num_processes', type=int, default=10, help="num_processes")
    parser.add_argument('-t', dest='summary_type', type=str, default="raw", choices=["raw"], help="summary_type")
    parser.add_argument('-X', dest='hit_pattern_to_not_match', type=str, default = None , help="hit regexp to not match")    
    parser.add_argument('-I', dest='hit_pattern_to_match', type=str, default = None , help="hit regexp to match")    
    parser.add_argument('-o', '--output_filename' , dest='output_filename', default="spectra.txt", type=str, help="name of the output file (default 'spectra.txt')")
    
    args = vars(parser.parse_args())

    # output file should not already exist
    if os.path.exists(args["output_filename"]):
        parser.error("error output file %(output_filename)s already exists"%args)

    
    return args

        
    
def main():
    args=get_options()

    #test(args["filename"][0])
    #return

    if args["action"] in ["build_first", "build_good", "list"]:
        for filename in  args["filename"]:

            if args["action"] == "build_first":
                dist = build_spectrum(filename, args, target="first_hit")
            if args["action"] == "build_good":
                dist = build_spectrum(filename, args, target="good_hits")                
            elif args["action"] == "list":
                dist = list_spectrum(filename, args)
    elif args["action"] == "summarise":
        summarise_spectra(args)
        #print dist
    elif args["action"] == "query":
        query_spectra(args)
    else:
        return

                                
if __name__ == "__main__":
   main()



        

