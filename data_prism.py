from __future__ import print_function
import re
import itertools
from multiprocessing import Pool, cpu_count
import pickle
import os
import string
import sys
if sys.version_info <= (2, 8):
   from exceptions import Exception
import gzip
import csv
import math
from types import GeneratorType


PROC_POOL_SIZE=30

class data_prism_exception(Exception):
    def __init__(self,args=None):
        super(data_prism_exception, self).__init__(args)



class prism(object):
    """
    this class is used to build a "spectrum" data structure, which is often but not alwaysa discrete multivariate
    frequency/probability distribution, from large input data (for example an alignment file of nextgen sequence data),
    containing mixed continuous and discrete multivariate data (the continuous variables are binned).
    """

    DEBUG = False
    input_streams = None
    """this allows a user to pass in one or more "pre-cooked" python value providers - i.e.
    arrays of tuples. The tuples provided must be consistent with the tuples that the
    file_to_stream_func function yields when applied to the input files. Note however that
    if you pass in an input stream, you can't use the multithreaded build method,
    as streams can't be pickled"""

    interval_locator_parameters = []
    """array with either length zero or same length as the number of dimensions. (Specify this if you use the
    built-in interval locators)"""
    
    

    def __init__(self,input_filenames, part_count = 1, input_streams=None):
        """
        prism constructor
       
        :param input_filenames: a list of filenames
        :param part_count: usually 1 but may be greater than 1 for large spectra
       
        """
       
        super(prism, self).__init__()

        self.name = "noname"
        self.input_filenames = input_filenames
        if self.input_filenames is None:
            self.input_filenames = []

        self.input_streams = input_streams
                                        
        self.part_count = part_count
        self.total_spectrum_value = 0
        self.approximate_zero = None
        self.interval_locator_parameters = []   
        
        self.interval_locator_funcs = []  # an array of functions suitable for locating which multi-dimensional interval a given value-tuple
                                          # should be assigned to. There should be same length as interval_locator_parameters
        self.assignments_files = [] # should be same length as interval_locator_parameters
        self.file_to_stream_func = None   # this should yield the tuples that are the domain of the spectrum
        self.file_to_stream_func_xargs = []        
        self.spectrum_value_provider_func = default_spectrum_value_provider
        self.spectrum_value_provider_func_xargs = []
        self.spectrum = None


    def summary(self, details = False):
        """

        print a summary of the spectrum
        
        :param details: true for more details, false for less
        
        """
        if details:
            print("""
               spectrum summary:
               
               total interval spectrum_value =     %s
               dimension count =        %s
               part_count =             %s

               details :
     
               input filenames:         %s
               dimension bins:          %s
               interval locators:       %s
               assignments files:       %s
               value provider:          %s
               value provider xargs:    %s  
            """%(self.total_spectrum_value, len(self.interval_locator_funcs), self.part_count,\
                 str(self.input_filenames), str(self.interval_locator_parameters), str(self.interval_locator_funcs),\
                 str(self.assignments_files),str(self.file_to_stream_func),\
                 str(self.file_to_stream_func_xargs)))
        else:
            print("""
               spectrum summary:
               
               total interval spectrum_value =     %s
               dimension count =        %s
               part_count =             %s
               approx_zero =            %s
               input filenames:         %s
               
            """%(self.total_spectrum_value, len(self.interval_locator_funcs), self.part_count,self.approximate_zero,\
                 str(self.input_filenames)))
            
    def list(self):
        spectrumdata = self.get_spectrum()
        for (interval, value) in spectrumdata.items():
            print(interval, value)
            
             
    def check_settings(self):
        if len(self.input_filenames) > 0 and self.file_to_stream_func is None:
            raise data_prism_exception("""
            Error - you have specified one or more input files, so you also need to set a file_to_stream_func - e.g.
            
            from data_prism import from_tab_delimited_file
            my_prism.file_to_stream_func = from_tab_delimited_file    # e.g. read a tab-delimited file 
            my_prism.file_to_stream_func_xargs = [1,3,6]              # e.g. use fields 1 3 and 6 
            """)
            
        if len(self.interval_locator_funcs) != len(self.interval_locator_parameters) and len(self.interval_locator_parameters) > 0:
            raise data_prism_exception("Error - there are %d interval locator functions specified, but %d dimensions"%(len(self.interval_locator_funcs),len(self.interval_locator_parameters)))
        if len(self.assignments_files) > 0:
            if len(self.assignments_files) != len(self.interval_locator_funcs):
                raise data_prism_exception("Error - there are %d interval assignments files specified, but %d dimensions"%(len(self.assignments_files),len(self.interval_locator_funcs)))
        
        
    def get_spectrum_values_stream(self):
        file_streams = [self.file_to_stream_func(filename, *self.file_to_stream_func_xargs) for filename in self.input_filenames]
        if self.input_streams is not None:
            #print "(using additional %d input readers with record counts : %s)"%(len(self.inpucount_providert_streams), str(map(len, self.input_streams)))
            print("(using additional %d input readers)"%len(self.input_streams))

            all_streams = file_streams + self.input_streams
            all_streams = itertools.chain(*all_streams)
        else:
            all_streams = itertools.chain(*file_streams)
        return all_streams

    def get_containing_interval(self, interval_value):
        return tuple(map(lambda f,value,interval:f(value,interval), self.interval_locator_funcs, interval_value, self.interval_locator_parameters))
        
    def get_partial_spectrum(self, slice_number):
        """
        this obtains a partial spectrum for one slice of the input file - this allows multiprocessing
        of an input file , as well as reducing intermediate dictionary sizes
        """
        #print "DEBUG in partial spectrum"

        self.check_settings()
        
        spectrum_values = self.get_spectrum_values_stream()

        myslice=itertools.islice(spectrum_values, slice_number, None, self.part_count) 

        # in this block we do a raw summary of the input data - without at this stage
        # assigning it to the final spectrum intervals. This minimises the number of calls
        # to the interval locator function. 
        sparse_data_summary = {}
        for interval in myslice:
            if self.DEBUG:
                print("DEBUG raw value : %s"%str(interval))

            for spectrum_value_tuple in self.spectrum_value_provider_func(interval, *self.spectrum_value_provider_func_xargs):
                if self.DEBUG:
                    print("DEBUG spectrum_value_tuple = %s length %d"%(str(spectrum_value_tuple), len(spectrum_value_tuple)))

                if len(spectrum_value_tuple) != 1+len(self.interval_locator_funcs):
                    raise data_prism_exception("Error - I received a spectrum_value_tuple %s but there are %d interval locators for locating the value"%(str(spectrum_value_tuple), len(self.interval_locator_funcs)))
                spectrum_value = float(spectrum_value_tuple[0])
                spectrum_interval = spectrum_value_tuple[1:]
                if self.DEBUG:
                    print("interval tuple = %s spectrum_value = %s"%(spectrum_interval, spectrum_value))
                if type(spectrum_interval) != tuple:
                    raise data_prism_exception("Error - I got %s from the spectrum_value provider: interval should be a tuple , instead it is %s (%s)"%(str(spectrum_value), str(spectrum_interval), type(spectrum_interval)))
                
                sparse_data_summary[spectrum_interval] = spectrum_value  +  sparse_data_summary.setdefault(spectrum_interval,0)
                

        # now we process the raw summary from the above, to accumulate the spectrum in the
        # final spectrum intervals
        partial = {}
        total_spectrum_value = 0
        assignment_writers = []
        if len(self.assignments_files) > 0:
            assignment_writers = [open("%s.%d"%(assignments_file,slice_number) ,"w") for assignments_file in self.assignments_files]

        for (sparse_key, sparse_total) in sparse_data_summary.items():

            if len(sparse_key) != len(self.interval_locator_funcs):
                print("warning  - interval to map (%s) is %d dimensional but %d locators are specified"%(str(sparse_key), len(sparse_key), len(self.interval_locator_funcs)))
                continue
                
            interval = self.get_containing_interval(sparse_key)
            
            if len(self.assignments_files) > 0:
                map(lambda raw, assignment, writer: writer.write("%s\t%s\n"%(raw, assignment)), sparse_key , interval, assignment_writers)

            partial[interval] = partial.setdefault(interval,0) + sparse_total
            
            total_spectrum_value += sparse_total
        if len(assignment_writers) > 0:
            map(lambda writer:writer.close(), assignment_writers)
        
        return (slice_number, partial)

    def get_spectrum(self):
        """
        this method obtains the complete spectrum, by adding together the partial spectra
        """
        self.spectrum = {}
        self.total_spectrum_value = 0
        for part in self.part_dict.values():
            for interval in part:
                self.spectrum[interval] = self.spectrum.setdefault(interval,0) + part[interval]
                self.total_spectrum_value += part[interval]

        # calulate an "approximate_zero"  - it is half the minimum raw value in any interval of the spectrum
        values = self.spectrum.values()
        if len(values) > 0:
            #self.approximate_zero = max(0.5, min(self.spectrum.values())/2.0)
            self.approximate_zero = 0.01
        else:
            #self.approximate_zero = 0.5
            self.approximate_zero = 0.01
            self.total_spectrum_value = self.approximate_zero 

        return self.spectrum

    
    def save(self, filename):
        print("saving prism object to %s"%filename)
        #print "object contains : %s"%dir(self)
        
        input_streams = self.input_streams
        
        self.input_streams = None 
        
        pwriter = open(filename, "wb")
        pickle.dump(self, pwriter)
        pwriter.close()
        
        self.input_streams = input_streams
        
    @staticmethod
    def load(filename):
        return p_load(filename)
            
    def get_spectrum_value(self,interval, default_spectrum_value = 0, return_interval = False ):
        interval = self.get_containing_interval(interval)
        if return_interval:
            return (self.spectrum.get(interval, default_spectrum_value), interval)
        else:
            return self.spectrum.get(interval, default_spectrum_value)
   
    def get_information(self,interval, method = None):
        spectrum_value = self.get_spectrum_value(interval,0)
        if spectrum_value > 0:
            return -1.0 * math.log(self.get_spectrum_value(interval) / float(self.total_spectrum_value))
        else:
            return -1.0 * math.log(1.0 / (1.0 + float(self.total_spectrum_value)))
            
    def get_information_projection(self, intervals, return_intervals = False):
        return p_get_information_projection((self, intervals, return_intervals))

    def get_signed_information_projection(self, intervals, return_intervals = False):
        return p_get_signed_information_projection((self, intervals, return_intervals))

    def get_unsigned_information_projection(self, intervals, return_intervals = False):
        return p_get_unsigned_information_projection((self, intervals, return_intervals))

    def get_raw_projection(self, intervals, return_intervals = False):
        return p_get_raw_projection((self, intervals, return_intervals))


    @staticmethod
    def get_intervals(spectrum_names,proc_pool_size=PROC_POOL_SIZE):
        """
        this method gets a union of all the intervals from a number of spectra.
        (Since it is returned as a list, this should be in a consistent order from call to call)
        """
        pool = Pool(proc_pool_size)
        spectra = pool.map(p_load,spectrum_names)
        intervals = set()
        for spectrum in spectra:
            intervals |= set(spectrum.get_spectrum().keys())

        return sorted(list(intervals))
        

    @staticmethod
    def get_projections(spectrum_names, intervals, projection_type, return_intervals = False, proc_pool_size = PROC_POOL_SIZE):
        """
        this method gets projections of a set of intervals across multiple spectra, returning
        the projections as a list of lists of (value1, value2, value3, etc) tuples 
        """
        print("distributing projections across %d processes"%proc_pool_size)        
        pool = Pool(proc_pool_size)
        print("get_projections : loading %s"%str(spectrum_names))
        spectra = pool.map(p_load,spectrum_names)
        print("done loading %d spectra"%len(spectra))

        args = zip(spectra, [intervals for spectrum in spectra], [return_intervals for spectrum in spectra])

        if projection_type == "raw": 
            projections = pool.map(p_get_raw_projection, args)
        elif projection_type == "unsigned_information":
            projections = pool.map(p_get_unsigned_information_projection, args)
        elif projection_type == "signed_information":
            projections = pool.map(p_get_signed_information_projection, args)
        elif projection_type == "information":
            projections = pool.map(p_get_information_projection, args)
        else:
            raise data_prism_exception("projection type %s not supported"%projection_type)

        return projections

    @staticmethod
    def save_projections(spectrum_names, intervals, projections, filename):
       """
       this method saves projections (a list of lists of tuples) as a tab-delimited text file including
       row and column names. There is one column per element of spectrum_names, and one row per member of
       intervals. The interval list is assumed to be in the same order as the project tuple list.
       The projections would typically be obtained by a call to get_projections( . . .return_intervals = False)
       """

       columnname_iter = [tuple([spectrum_name for spectrum_name in spectrum_names])]
       
       # this yields a row iterator, with the first row being column headings 
       row_iter = itertools.chain(columnname_iter, itertools.izip(*projections))

       # make a rowname iterator , including column header
       rowname_iter = itertools.chain([("interval",)],intervals)

       row_and_rowname_iter = itertools.izip(rowname_iter, row_iter)

       with open(filename, "w") as outfile:
          for row in row_and_rowname_iter:
             print("%s\t%s"%(":".join(row[0]), string.join((str(item) for item in row[1]),"\t")), file=outfile)

    @staticmethod
    def load_projections(filename, colnames_only = False):
       """
       this method loads from a saved file of projections, stipping column and row names and
       returns a row iterator
       """
       with open(filename,"r") as instream:
          if colnames_only:
             record_tuples = (re.split("\t", record.strip())[0] for record in instream)
             record_tuples.next()   
             return list(record_tuples)             
          else:  
             record_tuples = (re.split("\t", record.strip())[1:] for record in instream)
             record_tuples.next()   
             return list(record_tuples)


    @staticmethod
    def query_spectra(target_spectra, query_spectra, target_names):
       """ this method emits the euclidean distance between each spectrum in
       target_spectra, with each spectrum on query_spectra. This can be used to query
       a set of target spectra, with a given query spetrum - the output would be
       sorted and the "matches" are those with the shortest distance 
       """
       query_results = []
       name_iter = (name for name in target_names)
       for target in target_spectra:
          target_name = name_iter.next() 
          for query  in query_spectra:
             query_results.append((target_name, math.sqrt(reduce(lambda x,y:x+y, map(lambda q,t:(float(q)-float(t))**2, query, target)))))

       
       return sorted(query_results, lambda x,y:cmp(x[1],y[1]))
                   
             
      
#################################################
# distance-related  methods
# these methods don't really belong in this class here as a convenience.
# They calculate distances between column vectors of a matrix acccording
# to various metrics related to the spectrum_value distributions that are the
# subject of this class
#################################################

    @staticmethod
    def get_rank_iter(space_iter):
        """
        This method takes a matrix and replaces each value in the column with its rank in the column,
        and returns the matrix
        """
        interval_names = space_iter.next()
        if prism.DEBUG:
            print("**** DEBUG ranking")            
            print(str(interval_names))
            
        space_tuples = list(space_iter)
        if prism.DEBUG:
            print("**** DEBUG ranking")
            print(str(space_tuples))

        ranked_columns=[]
        for iinterval in range(0,len(interval_names)):
            index = range(0,len(space_tuples))

            if prism.DEBUG:
                print("**** DEBUG ranking")
                print(index)
                for i in range(0,len(space_tuples)):
                    print(iinterval, i, space_tuples[i][iinterval])
                
                
            ordered_index = sorted(index, lambda index1, index2: cmp(space_tuples[index1][iinterval], space_tuples[index2][iinterval]))
            ranks = sorted(index, lambda index1, index2: cmp(ordered_index[index1], ordered_index[index2]))

            #make 1-based ranks
            ranks=map(lambda x:1+x,ranks)
            
            ranked_columns.append(ranks)

        if prism.DEBUG:
            print("**** DEBUG ranking")
            print(ranked_columns)

        return itertools.chain([interval_names], itertools.izip(*ranked_columns)) 

    @staticmethod
    def get_distance_matrix(space_iter, method="euclidean"):
        """
        this method doesn't really belong in this class but is here as a convenience.
        This calculates the distances between column vectors of a matrix , where each
        column is probably an entropyprojecttion - but doens't have to be.
        Each column is headed up by a name of the column.
        """

        if method == "euclidean":
            return get_euclidean_distance_matrix(space_iter)
        elif method == "zipfian":
            return get_zipfian_distance_matrix(space_iter)
        else:
            raise data_prism_exception("unknown distance method %s"%method)
        
    @staticmethod
    def get_euclidean_distance_matrix(space_iter):
        """
        this method doesn't really belong in this class but is here as a convenience.
        This calculates the distances between column vectors of a matrix , where each
        column is probably an entropy projecttion - but doens't have to be.
        Each column is headed up by a name of the column.
        """
        distance_matrix = {}
        
        interval_names = space_iter.next()
        pair_names = [pair for pair in itertools.combinations(interval_names,2)]
        dimension_count = 0
        for dimension in space_iter:  # each row is a dimension
            dimension_count += 1
            #print "DEBUG %s"%str(dimension)
            dimension_pairs = [ dimension_pair for dimension_pair in itertools.combinations(map(float, dimension),2)]
            dimension_distances = [ dimension_distance for dimension_distance in map(lambda x:(x[1]-x[0])**2, dimension_pairs)] 
            for i in range(0,len(pair_names)):
                distance_matrix[pair_names[i]] = distance_matrix.get(pair_names[i],0) + dimension_distances[i]

        for pair_name in distance_matrix:
            distance_matrix[pair_name]  = math.sqrt( distance_matrix[pair_name] )

        # make the matrix symmetric
        for pair_name in pair_names:
            distance_matrix[(pair_name[1], pair_name[0])] = distance_matrix[pair_name]
        # fill in diagonal
        for interval_name in interval_names:
            distance_matrix[(interval_name, interval_name)] = 0

        interval_names_sorted = sorted(interval_names)

        return (distance_matrix, interval_names_sorted)


    @staticmethod
    def get_zipfian_distance_matrix(space_iter, rank_iter):
        """
        this method doesn't really belong in this class but is here as a convenience.
        This calculates the distances between "zipfian functions" , which relate 
        frequency-rank to frequency. The metric is based on a standard inner product
        for pairs of functions, and is based on log(frequency) rather than frequency.
        It looks similar to Euclidean distance, except that each term in the sum
        is inverse weighted by rank.

        Here we regard each pair of columns J from space_iter and rank_iter as
        defining a function, consisting of the ordered pairs (rank_iter[i,J], space_iter[i,J])

        space_iter and rank_iter are isomorphic data structures, with rank_iter containing
        the ranks of each element of space_iter
        """
        distance_matrix = {}
        
        interval_names = space_iter.next()
        rank_iter.next() # throw away headings in the rank clone
        
        pair_names = [pair for pair in itertools.combinations(interval_names,2)]


        # obtain the space of zipfian functions. Each function is represented by a dictionary, with key=rank,
        # value = entropy. The space is a list of these dictionaries
        space_aslist = list(space_iter)
        rank_aslist = list(rank_iter)

        function_space = [dict(map(lambda a,b:(a[j], b[j]), rank_aslist, space_aslist)) for j in range(0,len(interval_names))]

        #print "DEBUG"
        #for key in function_space[0].keys():
        #    print function_space[0][key], function_space[1][key], function_space[2][key]


        def zipf_distance_func(f1,f2):
            #print "DEBUG1"
            #total = 0.0
            #for key in f1.keys():
            #    print key, (f1[key]-f2[key])**2 / float(key)
            #    total += (f1[key]-f2[key])**2 / float(key)
            #print "manual total = %s"%total
                
            d = math.sqrt(reduce(lambda dsum,rank: dsum + ((f1[rank]-f2[rank])**2 / float(rank))   ,f1.keys(), 0.0))
            return d

        # get distances between pairs of intervals            
        for pair_name in pair_names:
            pair_index = (interval_names.index(pair_name[0]), interval_names.index(pair_name[1]))
            #print "DEBUG"
            #print "processing %s ( %s )"%(str(pair_name), str(pair_index))
            distance_matrix[pair_name] = zipf_distance_func(function_space[pair_index[0]],function_space[pair_index[1]])
            #print "got %s"%distance_matrix[pair_name]

        # make the matrix symmetric
        for pair_name in pair_names:
            distance_matrix[(pair_name[1], pair_name[0])] = distance_matrix[pair_name]
        # fill in diagonal
        for interval_name in interval_names:
            distance_matrix[(interval_name, interval_name)] = 0

        interval_names_sorted = sorted(interval_names)

        return (distance_matrix, interval_names_sorted)
                                                            
        

    @staticmethod
    def print_distance_matrix(distance_matrix, interval_names_sorted, outfile=sys.stdout):   
        print(string.join([""] + interval_names_sorted, "\t"), file=outfile)
        for row_interval in interval_names_sorted:
            print(string.join([row_interval]+[str(distance_matrix[(row_interval, col_interval)]) for col_interval in interval_names_sorted], "\t"),file=outfile)



#################################################
# top level versions of spectrum methods
# for use in multiprocessing context
#################################################
def p_get_raw_projection(arg_tuple):

    (spectrum, intervals, return_intervals) = arg_tuple 

    projection = len(intervals) * [None]
    if return_intervals:
        intervals = len(intervals) * [None]
    
    index = 0
    for interval in intervals:
        if not return_intervals:
            projection[index] = spectrum.get_spectrum_value(interval,0)
        else:
            (projection[index], intervals[index]) = spectrum.get_spectrum_value(interval,0, True)
        
        #projection[total_spectrum_value] = spectrum.get_spectrum_value(interval,0)
        #total_spectrum_value += 1
        index += 1

    if not return_intervals:
        return projection
    else:
        return (projection, intervals)


def p_get_information_projection(arg_tuple):
    (spectrum, intervals, return_intervals) = arg_tuple
    projection = len(intervals) * [None]
    if return_intervals:
        intervals = len(intervals) * [None]
        
    missing_spectrum_value = 0
    index = 0
    total_spectrum_value = 0
    for interval in intervals:
        (spectrum_value, interval)= spectrum.get_spectrum_value(interval,0, True)
        if not return_intervals:
            projection[index] = spectrum_value
        else:
            (projection[index], intervals[index]) = (spectrum_value, interval)
            
        if projection[index] == 0:
            missing_spectrum_value += 1
            
        index += 1
        total_spectrum_value += spectrum_value

    missing_spectrum_value = float(missing_spectrum_value)

    # all intervals projected onto missing get spectrum_value "missing_spectrum_value"
    # calculate information content based on the basis frequencies grossed up for 
    # missing intervals 
    for index in range(len(projection)):
        if projection[index] == 0:
            projection[index] = missing_spectrum_value

        projection[index] = -1.0 * math.log(projection[index] / float(total_spectrum_value + missing_spectrum_value), 2.0)

    if not return_intervals:
        return projection
    else:
        return (projection, intervals)


def p_get_signed_information_projection(arg_tuple):
    (spectrum, intervals, return_intervals) = arg_tuple
    projection = len(intervals) * [None]
    if return_intervals:
        intervals = len(intervals) * [None]
        
    missing_spectrum_value = 0
    index = 0
    total_spectrum_value = 0
    for interval in intervals:
        (spectrum_value, interval)= spectrum.get_spectrum_value(interval,0, True)
        if not return_intervals:
            projection[index] = spectrum_value
        else:
            (projection[index], intervals[index]) = (spectrum_value, interval)
            
        if projection[index] == 0:
            missing_spectrum_value += 1
            
        index += 1
        total_spectrum_value += spectrum_value

    missing_spectrum_value = float(missing_spectrum_value)

    #if missing_spectrum_value == len(intervals):
    #    raise data_prism_exception("no interval maps - unable to calculate projection")

    # as above (unsigned) all intervals projected onto missing get spectrum_value "missing_spectrum_value"
    # but now the missing intervals do not affect the projection of other intervals. 
    # Instead missing intervals get a negative information measure calculated as
    # P/1-P log P
    # where P is the probability that a projection interval maps - i.e. is not missing
    P = (len(intervals) - missing_spectrum_value) / float(len(intervals))
    for index in range(len(projection)):
        if projection[index] == 0:
            projection[index] =  (P / float(1-P) ) * math.log(P,2.0)
        else:
            projection[index] = -1.0 * math.log(projection[index] / total_spectrum_value, 2.0)

    if not return_intervals:
        return projection
    else:
        return (projection, intervals)


def p_get_unsigned_information_projection(arg_tuple):
    (spectrum, intervals, return_intervals) = arg_tuple
    projection = len(intervals) * [None]
    if return_intervals:
        intervals = len(intervals) * [None]
        
    missing_spectrum_value = 0
    index = 0
    total_spectrum_value = 0
    for interval in intervals:
        (spectrum_value, interval)= spectrum.get_spectrum_value(interval,0, True)
        if not return_intervals:
            projection[index] = spectrum_value
        else:
            (projection[index], intervals[index]) = (spectrum_value, interval)
            
        if projection[index] == 0:
            missing_spectrum_value += 1
            
        index += 1
        total_spectrum_value += spectrum_value

    missing_spectrum_value = float(missing_spectrum_value)

    #if missing_spectrum_value == len(intervals):
    #    raise data_prism_exception("no interval maps - unable to calculate projection")

    # all intervals projected onto missing get spectrum_value of approximate_zero (e.g. often .5) 
    # (but do not affect the projection of other intervals)

    if total_spectrum_value == 0:
        total_spectrum_value = spectrum.approximate_zero

    for index in range(len(projection)):
        if projection[index] == 0:
            projection[index] = -1.0 * math.log(spectrum.approximate_zero / float(total_spectrum_value), 2.0)
        else:
            projection[index] = -1.0 * math.log(projection[index] / float(total_spectrum_value), 2.0)

    if not return_intervals:
        return projection
    else:
        return (projection, intervals)



def p_load(filename):
    preader = open(filename, "rb")
    pinstance = pickle.load(preader)
    preader.close()

    if not isinstance(pinstance, prism):
        raise data_prism_exception("%s is not a spectrum object"%filename )
    return pinstance

        
#################################################
# built-in locator functions, for locating which interval a value
# is in. These are usually too slow and it is usually best to supply your
# own locator function. However useful if intervals are of varying
# sizes, and there is not much data
#################################################
def bin_continuous_value(value, intervals):
    """
    method for locating a scalar numeric value within an array of
    numeric scalar intervals. The intervals array consists simply
    of an array of the lower bound of each interval - i.e.
    it contains contiguous intervals
    """
    if value is None:
        return None
    
    interval = None

    nvalue = float(value)

    for i in range(0,len(intervals)):
        if  intervals[i] > nvalue :
            if i > 0:
                interval = intervals[i-1]
            break
        
    return interval
            
def bin_discrete_value(value, intervals):
    """
    method for locating a discrete value within an array of
    discrete intervals. This simplest method simply returns the
    interval passed in, if it is contained in intervals (or nothing in
    intervals), else None
    """
    if value is None:
        return None

    if intervals is None:
        return value
    elif len(intervals) == 0:
        return value
    elif value in intervals:
        return value
    else:
        return None

#################################################
# convenience value provider functions, for providing a tuple of
# values to be accumulated in the density function. A factory method
# returns a value-provider closure
#################################################

class outer_list(list):        
    def __getitem__(self, key):
        if key >= self.__len__():
            return None
        else:
            return super(outer_list,self).__getitem__(key)

def default_spectrum_value_provider(interval, *xargs):
    """
    the default spectrum_value_provider function is suitable when used with a stream
    provider (e.g file_to_stream or a user-provided steam)  that yields the domain intervals of the
    spectrum.
    """
    return (((1,)+interval),)     # this needs to provide a tuple of 1 tuples (not just a tuple)
                               # (because some provider functions provide tuples of many tuples)

def get_file_type(file_path):
    """
    Make best efforts to infer  file type from file name for some common file formats, and return a canonical
    type name. For other file formats, just return the suffix (e.g.
    foo.bar and foo.BAR would both return .bar in that case). May not infer correctly for 
    some names (e.g. blah.fastq.fasta etc)
    """
    real_path=file_path

    if os.path.exists(file_path):
       real_path=os.path.realpath(file_path)
       
    if re.search("(\.fasta|\.fa|\.fna|\.faa|\.seq)(\.|$)", real_path, re.I) != None:
        return "fasta"
    elif re.search("(\.fastq|\.fq)(\.|$)", real_path, re.I) != None:
        return "fastq"
    else:
        return os.path.splitext(file_path)[1].lower()
    

def get_text_stream(file_path):
    """
    this is used by most of the value providers below that use text files
    (tab seperated, csv etc)
    """

    real_path=file_path

    if os.path.exists(file_path):
       real_path=os.path.realpath(file_path)
    
    text_stream = None
    if re.search("\.gz$", real_path, re.I) != None:
        text_stream = gzip.open(real_path, 'rb')
    else:
        text_stream = open(real_path,"r")
    
    return text_stream
    
    
def from_flat_file(file_name, *xargs):
    """
    basic method - just split each record on white space and return tuples
    - will need to provide dimension interval specs and locators for each
    element of the tuple
    """
    if len(xargs) == 0:
        return (tuple(re.split("\s+", record.strip())) for record in get_text_stream(file_name))
    else:
        return (tuple([outer_list(re.split("\s+", record.strip()))[index] for index in xargs])  for record in get_text_stream(file_name))

def from_tab_delimited_file(file_name, *xargs):
    """
    basic method - just split each record on white space and return tuples
    - will need to provide dimension interval specs and locators for each
    element of the tuple.

    Note that this may return "ragged" records as it strips leading and trailing whitespace (i.e. including tabs)
    """
    if len(xargs) == 0:
        return (tuple(re.split("\t", record.strip())) for record in get_text_stream(file_name))
    else:
        return (tuple([outer_list(re.split("\t", record.strip()))[index] for index in xargs])  for record in get_text_stream(file_name))


def from_nonragged_tab_delimited_file(file_name, *xargs):
    """
    basic method - just split each record on white space and return tuples
    - will need to provide dimension interval specs and locators for each
    element of the tuple
    """
    if len(xargs) == 0:
        return (tuple(re.split("\t", record.translate(None,"\n\r"))) for record in get_text_stream(file_name))
    else:
        return (tuple([outer_list(re.split("\t", record.translate(None,"\n\r")))[index] for index in xargs])  for record in get_text_stream(file_name))

    
def from_csv_file(file_name, *xargs):
    """
    basic method - just split each record on white space and return tuples
    - will need to provide dimension interval specs and locators for each
    element of the tuple
    """
    if len(xargs) == 0:
        return csv.reader(get_text_stream(file_name))
    else:
        return (tuple([outer_list(record)[index] for index in xargs])  for record in csv.reader(get_text_stream(file_name)))


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
    spectrum_value = args[2] # un-used currently 
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
            kmer_count_iter = ( ( len(list(kmer_iter)),kmer) for (kmer,kmer_iter) in itertools.groupby(kmer_iters, lambda kmer:kmer) )
        else:
            kmer_count_iter = ( ( len(list(kmer_iter)),get_reverse_complement(kmer)) for (kmer,kmer_iter) in itertools.groupby(kmer_iters, lambda kmer:kmer) )
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
                kmer_dict[kmer] = 1 + kmer_dict.setdefault(kmer,0)
                
        kmer_count_iter = ( (kmer_dict[kmer], kmer) for kmer in kmer_dict )
        
        
    return kmer_count_iter
        



#################################################
# build methods 
#################################################


def build_part(arg_tuple):
    (spectrum_instance, slice_number) = arg_tuple
    print("build_part is building part %d"%slice_number)
    return spectrum_instance.get_partial_spectrum(slice_number)

def build(spectrum_instance, use="multithreads", proc_pool_size = PROC_POOL_SIZE):

    if use == "cluster":
        args = [(spectrum_instance, slice_number)  for slice_number in range(0,spectrum_instance.part_count)]
        spectrum_instance.part_dict = dict(pool.map(build_part,args))

        return spectrum_instance.get_spectrum()
        
    elif use == "multithreads":
        pool = Pool(proc_pool_size)
        
        args = [(spectrum_instance, slice_number)  for slice_number in range(0,spectrum_instance.part_count)]

        print("mapping %s build parts to a pool of size %d"%(len(args), proc_pool_size))
        spectrum_instance.part_dict = dict(pool.map(build_part,args))

        return spectrum_instance.get_spectrum()

    elif use == "singlethread":
        
        args = [(spectrum_instance, slice_number)  for slice_number in range(0,spectrum_instance.part_count)]

        results = []
        for arg in args:
            results.append(build_part(arg))
            
        spectrum_instance.part_dict = dict(results)

        return spectrum_instance.get_spectrum()
    
    else:
        raise data_prism_exception("error - unknown resource specified for build : %s"%use)
    

        

