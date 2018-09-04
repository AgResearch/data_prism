#!/usr/bin/env python2.7

import itertools,os,re,argparse,string,sys
#sys.path.append('/usr/local/agr-scripts')
#from prbdf import Distribution , build, from_tab_delimited_file, bin_discrete_value
from data_prism import prism, build, from_tab_delimited_file, bin_discrete_value


def my_hit_provider(filename, *xargs):
    """
    transform the tab-delimited stream, to only yield the records that relates either to a hit
    or "no hit" . Note that sometimes this format reports multiple hits to the same target
    - we only want the top hit - this is provided by the next method
    """
    weighting_method = xargs[0]
    tuple_stream = from_tab_delimited_file(filename,*xargs[1:])

    atuple = tuple_stream.next()
    query = ""
    while True:
        weight = 1
        query_match = re.search("^#\s+Query:\s+(.*)$",atuple[0].strip())
        if query_match is not None:
            query = query_match.groups()[0]
        if re.search(" 0 hits",atuple[0],re.IGNORECASE) is not None:
            if weighting_method == "tag_count":
                weighting_match = re.search("count=(\d*\.*\d*)\s*$", query)
                weight = float(weighting_match.groups()[0])
            yield ((query,'No hits','No hits'),weight)
        elif atuple[1:] == (None,None):
            pass
        elif atuple[1] is None or atuple[2] is None:
            raise Exception("error - unexpected results %s from blast output - incomplete taxonomy tuple"%str(atuple))        
        else:
            if weighting_method == "tag_count":
                weighting_match = re.search("count=(\d*\.*\d*)\s*$", query)
                weight = float(weighting_match.groups()[0])
            yield (atuple, weight)
        
        atuple = tuple_stream.next()

def my_top_hit_provider(filename, *xargs):
    """
    takes a stream which may contain multiple hits, and yields just the top hit in each group
    """
    groups = itertools.groupby(my_hit_provider(filename, *xargs), lambda x:x[0][0])    
    top_hits = (group.next() for (key, group) in groups)
    return top_hits


def my_spectrum_value_provider(interval_weight, *xargs):
    """
    this takes the items from top hit provider - e.g.
    (('seq_25449', 'Eukaryota', 'Pacific oyster'), 52)
    and transforms to e.g.
    ((52, 'seq_25449', 'Eukaryota', 'Pacific oyster'),)
  
    """
    #print interval_weight
    return ((interval_weight[1],interval_weight[0][1],interval_weight[0][2]),)       

def build_tax_distribution(datafile, weighting_method = None):
    distob = prism([datafile], 1)
    distob.file_to_stream_func = my_top_hit_provider
    #distob.DEBUG = True
    distob.file_to_stream_func_xargs = [weighting_method,0,7,6] # i.e. pick out first field, then kingdom, comnames
    distob.interval_locator_funcs = [bin_discrete_value, bin_discrete_value]
    distob.spectrum_value_provider_func = my_spectrum_value_provider
    distdata = build(distob,"singlethread")
    distob.save("%s.pickle"%datafile)
    return distdata

def tax_cmp(x,y):
    ord=cmp(x[0],y[0])
    if ord == 0:
        ord = cmp(x[1], y[1])
    return ord

def get_sample_tax_distribution(sample_tax_summaries, measure,rownames):
    sample_tax_lists = [ prism.load(sample_tax_summary).get_spectrum().keys() for sample_tax_summary in sample_tax_summaries ] 
    all_taxa = set( reduce(lambda x,y:x+y, sample_tax_lists))
    all_taxa_list = list(all_taxa)
    all_taxa_list.sort(tax_cmp)

    #print all_taxa_list

    if measure == "frequency":
        if not rownames:
            sample_tax_distributions = [[re.sub("'|#","","%s\t%s"%item) for item in all_taxa_list]] + [ prism.load(sample_tax_summary).get_raw_projection(all_taxa_list) for sample_tax_summary in sample_tax_summaries]
        else:
            sample_tax_distributions = [[re.sub("'|#","","%s_%s"%item) for item in all_taxa_list]] + [ prism.load(sample_tax_summary).get_raw_projection(all_taxa_list) for sample_tax_summary in sample_tax_summaries]

    else:
        if not rownames:
            sample_tax_distributions = [[re.sub("'|#","","%s\t%s"%item) for item in all_taxa_list]] + [ prism.load(sample_tax_summary).get_unsigned_information_projection(all_taxa_list) for sample_tax_summary in sample_tax_summaries]
        else:
            sample_tax_distributions = [[re.sub("'|#","","%s_%s"%item) for item in all_taxa_list]] + [ prism.load(sample_tax_summary).get_unsigned_information_projection(all_taxa_list) for sample_tax_summary in sample_tax_summaries]



    fd_iter = itertools.izip(*sample_tax_distributions)
    if not rownames:
        heading = itertools.izip(*[["Kingdom\tFamily"]]+[[re.split("\.",os.path.basename(path.strip()))[0]] for path in sample_tax_summaries])
    else:
        heading = itertools.izip(*[["Kingdom_Family"]]+[[re.split("\.",os.path.basename(path.strip()))[0]] for path in sample_tax_summaries])

    #print heading

    fd_iter = itertools.chain(heading, fd_iter)

    for record in fd_iter:
        print string.join([str(item) for item in record],"\t")


def write_summaries(basename,tax_dist):
    """
    write out a top level summary file giving totals for each kingdom
    and then a file for each kingdom giving totals for each family in taht kingdom
    """
    tax_list = tax_dist.items()
    """ e.g.
    {('Eukaryota', 'African rice'): 2, ('Bacteria', 'Clavibacter michiganensis subsp. sepedonicus'): 4
    """
    # group by kingdom    
    tax_list.sort(lambda x,y:cmp(x[0][0],y[0][0])) # first sort by kingdom
    kingdom_iter = itertools.groupby(tax_list, lambda x:x[0][0])

    summary_path = "%s.kingdom_summary.txt"%basename
    summary_writer = open(summary_path,"w")
    print >> summary_writer, "kingdom","\t","count"
    for (kingdom, family_iter) in kingdom_iter:
        family_list = list(family_iter)
        family_list.sort(lambda x,y:cmp(x[0][1], y[0][1]))
        kingdom_total = reduce(lambda x,y : x+y, (record[1] for record in family_list))
        print >> summary_writer, kingdom,"\t", kingdom_total

        kingdom_filename="%s.family_%s_summary.txt"%(basename,re.sub("[^A-Za-z]","_",kingdom))
        kingdom_path = kingdom_filename       
        kingdom_writer = open(kingdom_path,"w")
        print >> kingdom_writer, "family","\t","count"
        for record in family_list:
            print >> kingdom_writer, record[0][1],"\t", record[1]
        kingdom_writer.close()
    summary_writer.close()

def debug(options):
    #test_iter = my_hit_provider(options["filenames"][0], *[None,0,7,6])
    test_iter = my_top_hit_provider(options["filenames"][0], *["tag_count",0,7,6])

    for item in test_iter:
        print item
        #print my_spectrum_value_provider(item, *[])

class outer_list(list):        
    def __getitem__(self, key):
        if key >= self.__len__():
            return None
        else:
            return super(outer_list,self).__getitem__(key)

def get_options():
    description = """
    """
    long_description = """

example :

./summarise_hiseq_taxonomy.py  /dataset/hiseq/scratch/postprocessing/Salmon_mixed_runs.processed_in_progress/taxonomy_in_progress/Project_Salmon_HalfVol_ApeKI_Sample_SQ0031.list.nt_blastresults.txt.gz

where input file is generated by the following blast format string :

"7 qseqid sseqid pident evalue staxids sscinames scomnames sskingdoms stitle"

- for example

# BLASTN 2.2.28+
# Query: INV-D00390:220:C6GKKANXX:8:1114:13520:88778 1:N:0:
# Database: nt
# 0 hits found
# BLASTN 2.2.28+
# Query: INV-D00390:220:C6GKKANXX:8:1114:18656:89413 1:N:0:
# Database: nt
# Fields: query id, subject id, % identity, evalue, subject tax ids, subject sci names, subject com names, subject super kingdoms, subject title
# 1 hits found
INV-D00390:220:C6GKKANXX:8:1114:18656:89413     gi|688443106|emb|LL194098.1|    100.00  9e-18   6339    Heligmosomoides polygyrus       Heligmosomoides polygyrus      Eukaryota        Heligmosomoides polygyrus genome assembly H_bakeri_Edinburgh ,scaffold HPBE_scaffold0005644
.
.
.
optionally, the query line can specify a weighting to be used as a count instead of 1 - e.g.
# Query: seq_26674 count=16

(this is used when blasting queries such as unique tags)
"""

    parser = argparse.ArgumentParser(description=description, epilog=long_description, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filenames', type=str, nargs="*",help='input file of blast hits (optionally compressed with gzip)')    
    parser.add_argument('--summary_type', dest='summary_type', default="sample_summaries", \
                   choices=["sample_summaries", "summary_table"],help="summary type (default: sample_summaries")
    parser.add_argument('--measure', dest='measure', default="frequency", \
                   choices=["frequency", "information"],help="measure (default: frequency")
    parser.add_argument('--rownames' , dest='rownames', default=False,action='store_true', help="combine kingdom and family fields to make a rowname")
    parser.add_argument('--weighting_method' , dest='weighting_method', default=None,choices=["tag_count"],help="weighting method")


    args = vars(parser.parse_args())
    return args

        
    
def main():
    args=get_options()

    #test = my_top_hit_provider(filename, 0,7,6)
    #test = my_hit_provider(filename, 0,7,6)
    #for record in test:
    #    print record

    #return
    #debug(args)

    if args["summary_type"] == "sample_summaries" :
        for filename in  args["filenames"]:
            tax_dist = build_tax_distribution(filename, weighting_method = args["weighting_method"])
            print tax_dist
            write_summaries(filename,tax_dist)
    elif args["summary_type"] == "summary_table" :
        #print "summarising %s"%str(args["filename"])
        get_sample_tax_distribution(args["filenames"], args["measure"], args["rownames"])

    

    return

                                
if __name__ == "__main__":
   main()



        

