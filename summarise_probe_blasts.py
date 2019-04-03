#!/usr/bin/env python
from __future__ import print_function

import itertools,os,re,argparse,string,sys
sys.path.append('/dataset/bioinformatics_dev/active/data_prism')
from data_prism import prism , build, bin_discrete_value, get_text_stream , get_file_type,  PROC_POOL_SIZE

def my_hit_provider(filename, *xargs):
    """
    provide hits from a standard blastn output , like
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
        if atuple[0] == ">":
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

def my_top_hit_provider(filename, *xargs):
    """
    takes a stream which may contain multiple hits, and yields just the top hit in each group
    """
    groups = itertools.groupby(my_hit_provider(filename, *xargs), lambda x:x[0])    
    top_hits = (group.next() for (key, group) in groups)
    return top_hits

def my_value_provider(point, *xargs):
    # top hit provider returns query, statistic  - weight value provider returns statistic, queryid
    (query, identity_string ) = point   # e.g. NC_002944.2_1-100   ,(100%)
    value=identity_string.translate(None, ",()%") 
    value_tuple = (value, query)    
    return (value_tuple,)


def test(datafile):
    #s = my_hit_provider(datafile, 0)
    #for t in s:
    #    print(t)
    s=my_top_hit_provider(datafile,0)
    for t in s:
        print(my_value_provider(t))


def build_spectrum(datafile, options):
    hit_prism = prism([datafile], 1)
    #distob.DEBUG=True
    hit_prism.file_to_stream_func = my_top_hit_provider
    hit_prism.file_to_stream_func_xargs = [options["hit_pattern_to_match"], options["hit_pattern_to_not_match"]]
    hit_prism.spectrum_value_provider_func = my_value_provider

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
    

    #print("zipping projections...")
    #zprobe_measures = itertools.izip(*probe_measures)
    #probe_name_iter = [tuple([os.path.splitext(os.path.basename(spectrum))[0] for spectrum in spectrum_names])]
    #zprobe_measures = itertools.chain(probe_name_iter, zprobe_measures)
    #interval_name_iter = itertools.chain([("probe")],probe_intervals)
    
    ##outfile=open(options["output_filename"], "w")

    #if options["summary_type"] in ["raw"]:
    #    zprobe_measures_with_rownames = itertools.izip(interval_name_iter, zprobe_measures)
    #    for interval_measure in zprobe_measures_with_rownames:
    #        #print("%s\t%s"%("%s"%interval_measure[0], string.join((str(item) for item in interval_measure[1]),"\t")))
    #        print("%s\t%s"%("%s"%interval_measure[0], string.join((str(item) for item in interval_measure[1]),"\t")), file=open(options["output_filename"], "a"))
    #        #with open(outfile, "a") as outfile:
    #        #    outfile.writelines("%s\t%s\n"%("%s"%interval_measure[0], string.join((str(item) for item in interval_measure[1]),"\t")))
    #else:
    #    print("warning, unknown summary type %(summary_type)s, no summary available"%options)
        
        
    #   outfile.close()



def get_options():
    description = """
    """
    long_description = """
    """

    parser = argparse.ArgumentParser(description=description, epilog=long_description, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', type=str, nargs="*",help='input file of blast hits (optionally compressed with gzip)')    
    parser.add_argument('--use', dest='use', default="use", choices=["pctidentity"],help="use")
    parser.add_argument('-a', dest='action', default="build", choices=["build", "list", "summarise"],help="action")
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

    if args["action"] in ["build", "list"]:
        for filename in  args["filename"]:

            if args["action"] == "build":
                dist = build_spectrum(filename, args)
            elif args["action"] == "list":
                dist = list_spectrum(filename, args)
    elif args["action"] == "summarise":
        summarise_spectra(args)
        
        #print dist

    return

                                
if __name__ == "__main__":
   main()



        

