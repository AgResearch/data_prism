import re
import itertools
from multiprocessing import Pool, cpu_count
import pickle
import os
import string
import exceptions
import gzip
import csv
import math

PROC_POOL_SIZE=30

class prbdfException(exceptions.Exception):
    def __init__(self,args=None):
        super(prbdfException, self).__init__(args)



class Distribution(object):
    """
    this class is used to build a discrete multivariate frequency/probability Distribution, from very
    large input data (for example an alignment file of nextgen sequence data), containing
    mixed continuous and discrete multivariate data (the continuous variables are binned)    
    """

    
    def __init__(self,input_filenames, part_count = 1, input_streams=None):
        super(Distribution, self).__init__()

        self.name = "noname"
        self.input_filenames = input_filenames
        if self.input_filenames is None:
            self.input_filenames = []
            
        self.input_streams = input_streams # this allows a user to pass in one or more "pre-cooked" python value providers - i.e.
                                           # arrays of tuples. The tuples provided must be consistent with the tuples that the
                                           # file_to_stream_func function yields when applied to the input files. Note however that
                                           # if you pass in an input stream, you can't use the multithreaded build method,
                                           # as streams can't be pickled 
                                        
        self.part_count = part_count
        self.point_count = 0
        self.approximate_zero_frequency = None
        self.interval_locator_parameters = []   # Array with either length zero or same length as the number of dimensions. (Specify this if you use the
                                          # built-in interval locators)
        
        self.interval_locator_funcs = []  # an array of functions suitable for locating which multi-dimenisonal interval a given value-tuple
                                          # should be assigned to. There should be same length as interval_locator_parameters
        self.assignments_files = [] # should be same length as interval_locator_parameters
        self.file_to_stream_func = None
        self.file_to_stream_func_xargs = []        
        self.value_provider_func = default_value_provider
        self.count_provider_func = default_count_provider
        self.frequency_distribution = None
        self.DEBUG = False


    def summary(self, details = False):
        if details:
            print """
               distribution summary:
               
               point count =            %s
               dimension count =        %s
               part_count =             %s

               details :
     
               input filenames:         %s
               dimension bins:          %s
               interval locators:       %s
               assignments files:       %s
               value provider:          %s
               value provider xargs:    %s  
            """%(self.point_count, len(self.interval_locator_funcs), self.part_count,\
                 str(self.input_filenames), str(self.interval_locator_parameters), str(self.interval_locator_funcs),\
                 str(self.assignments_files),str(self.file_to_stream_func),\
                 str(self.file_to_stream_func_xargs))
        else:
            print """
               distribution summary:
               
               point count =            %s
               dimension count =        %s
               part_count =             %s
               approx_zero =            %s
               input filenames:         %s
               
            """%(self.point_count, len(self.interval_locator_funcs), self.part_count,self.approximate_zero_frequency,\
                 str(self.input_filenames))
            
    def list(self):
        distdata = self.get_distribution()
        for (interval, freq) in distdata.items():
            print interval, freq
            
             
    def check_settings(self):
        if len(self.input_filenames) > 0 and self.file_to_stream_func is None:
            raise prbdfException("""
            Error - you have specified one or more input files, so you also need to set a file_to_stream_func - e.g.
            
            from prbdf import from_tab_delimited_file
            my_dist_ob.file_to_stream_func = from_tab_delimited_file    # e.g. read a tab-delimited file 
            my_dist_ob.file_to_stream_func_xargs = [1,3,6]              # e.g. use fields 1 3 and 6 
            """)
            
        if len(self.interval_locator_funcs) != len(self.interval_locator_parameters) and len(self.interval_locator_parameters) > 0:
            raise prbdfException("Error - there are %d interval locator functions specified, but %d dimensions"%(len(self.interval_locator_funcs),len(self.interval_locator_parameters)))
        if len(self.assignments_files) > 0:
            if len(self.assignments_files) != len(self.interval_locator_funcs):
                raise prbdfException("Error - there are %d interval assignments files specified, but %d dimensions"%(len(self.assignments_files),len(self.interval_locator_funcs)))
        
        
    def get_raw_values_stream(self):
        file_streams = [self.file_to_stream_func(filename, *self.file_to_stream_func_xargs) for filename in self.input_filenames]
        if self.input_streams is not None:
            #print "(using additional %d input readers with record counts : %s)"%(len(self.inpucount_providert_streams), str(map(len, self.input_streams)))
            print "(using additional %d input readers)"%len(self.input_streams)

            #reader_iters =  [(self.value_provider_func(record) for record in reader) for reader in self.input_streams]
            all_streams = file_streams + self.input_streams
            all_streams = itertools.chain(*all_streams)
        else:
            all_streams = itertools.chain(*file_streams)
        return all_streams

    def get_containing_interval(self, point_value):
        return tuple(map(lambda f,value,interval:f(value,interval), self.interval_locator_funcs, point_value, self.interval_locator_parameters))
        
    def get_partial_distribution(self, slice_number):
        #print "DEBUG in partial distribution"

        self.check_settings()
        
        raw_values = self.get_raw_values_stream()

        myslice=itertools.islice(raw_values, slice_number, None, self.part_count) 
        
        sparse_data_summary = {}

        count = 0
        
        for point in myslice:
            if self.DEBUG:
                print "DEBUG raw value : %s"%str(point)
            sparse_data_summary[self.value_provider_func(point)] = self.count_provider_func(point) + sparse_data_summary.setdefault(self.value_provider_func(point),0)
            count += 1

        partial = {}
        count = 0
        point_count = 0
        assignment_writers = []
        if len(self.assignments_files) > 0:
            #print "DEBUG opening writers : %s"%str(self.assignments_files)
            assignment_writers = [open("%s.%d"%(assignments_file,slice_number) ,"w") for assignments_file in self.assignments_files]

        #print "DEBUG running assignments using %s"%str(self.interval_locator_funcs)
        for (sparse_key, sparse_total) in sparse_data_summary.items():
            #print "DEBUG3 %s %s"%(sparse_key, str(self.interval_locator_funcs))

            if len(sparse_key) != len(self.interval_locator_funcs):
                #raise prbdfException("error - interval to map is %d dimensional but %d locators are specified"%(len(sparse_key), len(self.interval_locator_funcs)))
                print "warning  - interval to map (%s) is %d dimensional but %d locators are specified"%(str(sparse_key), len(sparse_key), len(self.interval_locator_funcs))
                continue
                
            interval = self.get_containing_interval(sparse_key)
            
            if len(self.assignments_files) > 0:
                map(lambda raw, assignment, writer: writer.write("%s\t%s\n"%(raw, assignment)), sparse_key , interval, assignment_writers)
            
            count += 1

            partial[interval] = partial.setdefault(interval,0) + sparse_total
            
            point_count += sparse_total
        if len(assignment_writers) > 0:
            map(lambda writer:writer.close(), assignment_writers)
        
        return (slice_number, partial)

    def get_distribution(self):
        """
        returns a merge of the part pd's
        """
        self.frequency_distribution = {}
        self.point_count = 0
        for part in self.part_dict.values():
            for interval in part:
                self.frequency_distribution[interval] = self.frequency_distribution.setdefault(interval,0) + part[interval]
                self.point_count += part[interval]

        # calulate an "approximate_zero" frequency - it is half the minimum frequency in any interval
        self.approximate_zero_frequency = min(self.frequency_distribution.values())/2.0

        return self.frequency_distribution

    
    def save(self, filename):
        print "saving prbdf object to %s"%filename
        #print "object contains : %s"%dir(self)
        
        input_streams = self.input_streams
        #value_provider_filters = self.value_provider_filters
        #value_provider = self.value_provider
        #interval_locator_funcs = self.interval_locator_funcs
        
        self.input_streams = None 
        #self.value_provider_filters = {}
        #self.value_provider = None
        #self.interval_locator_funcs = []
        
        pwriter = open(filename, "wb")
        pickle.dump(self, pwriter)
        pwriter.close()
        
        self.input_streams = input_streams
        #self.value_provider_filters = value_provider_filters
        #self.value_provider = value_provider
        #self.interval_locator_funcs = interval_locator_funcs
        
    @staticmethod
    def load(filename):
        return p_load(filename)
        
    def get_frequency(self,point, default_frequency = 0):
        interval = self.get_containing_interval(point)
        return self.frequency_distribution.get(interval, default_frequency)

    def get_density(self,point, method = None):
        return self.get_frequency(point) / float(self.point_count)
    
    def get_information(self,point, method = None):
        frequency = self.get_frequency(point,0)
        if frequency > 0:
            return -1.0 * math.log(self.get_frequency(point) / float(self.point_count))
        else:
            return -1.0 * math.log(1.0 / (1.0 + float(self.point_count)))
            
    def get_information_projection(self, points):
        return p_get_information_projection((self, points))

    def get_signed_information_projection(self, points):
        return p_get_signed_information_projection((self, points))

    def get_unsigned_information_projection(self, points):
        return p_get_unsigned_information_projection((self, points))

    def get_frequency_projection(self, points):
        return p_get_frequency_projection((self, points))

    @staticmethod
    def get_projections(distribution_names, points, projection_type):
        """
        this method gets frequency projections of a set of points across multiple distributions
        """
        pool = Pool(PROC_POOL_SIZE)
        print "DEBUG about to load %s"%str(distribution_names)
        distributions = pool.map(p_load,distribution_names)
        print "DEBUG done loading "

        args = zip(distributions, [points for distribution in distributions])
        if projection_type == "frequency": 
            projections = pool.map(p_get_frequency_projection, args)
        elif projection_type == "unsigned_information":
            projections = pool.map(p_get_unsigned_information_projection, args)
        elif projection_type == "signed_information":
            projections = pool.map(p_get_signed_information_projection, args)
        elif projection_type == "information":
            projections = pool.map(p_get_information_projection, args)
        else:
            raise prbdfException("projection type %s not supported"%projection_type)

        return projections


    @staticmethod
    def get_distance_matrix(space_iter):
        """
        this method doesn't really belong in this class but is here as a convenience.
        This calculates the distances between column vectors of a matrix , where each
        column is probably a frequency projecttion - but doens't have to be.
        Each column is headed up by a name of the column.
        """
        distance_matrix = {}

        #print "DEBUG DIM"
        #print space_iter.next()
        #print space_iter.next()
        
        point_names = space_iter.next()
        pair_names = [pair for pair in itertools.combinations(point_names,2)]
        dimension_count = 0
        for dimension in space_iter:
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
        for point_name in point_names:
            distance_matrix[(point_name, point_name)] = 0

        point_names_sorted = sorted(point_names)

        return (distance_matrix, point_names_sorted)

    @staticmethod
    def print_distance_matrix(distance_matrix, point_names_sorted):   
        print string.join([""] + point_names_sorted, "\t")
        for row_point in point_names_sorted:
            print string.join([row_point]+[str(distance_matrix[(row_point, col_point)]) for col_point in point_names_sorted], "\t")



#################################################
# top level versions of Distribution methods
# for use in multiprocessing context
#################################################
def p_get_frequency_projection((distribution, points)):
    projection = len(points) * [None]
    point_count = 0
    for point in points:
        projection[point_count] = distribution.get_frequency(point,0)
        point_count += 1

    return projection


def p_get_information_projection((distribution, points)):
    projection = len(points) * [None]
    missing_count = 0
    index = 0
    for point in points:
        projection[index] = distribution.get_frequency(point,0)
        if projection[index] == 0:
            missing_count += 1
        index += 1

    missing_count = float(missing_count)

    # all points projected onto missing get frequency "missing_count"
    # calculate information content based on the basis frequencies grossed up for 
    # missing points 
    for index in range(len(projection)):
        if projection[index] == 0:
            projection[index] = missing_count

        projection[index] = -1.0 * math.log(projection[index] / float(distribution.point_count + missing_count), 2.0)

    return projection


def p_get_signed_information_projection((distribution, points)):
    projection = len(points) * [None]
    missing_count = 0
    index = 0
    for point in points:
        projection[index] = distribution.get_frequency(point,0)
        if projection[index] == 0:
            missing_count += 1
        index += 1

    missing_count = float(missing_count)

    if missing_count == len(points):
        raise prbdfException("no point maps - unable to calculate projection")

    # as above (unsigned) all points projected onto missing get frequency "missing_count"
    # but now the missing points do not affect the projection of other points. 
    # Instead missing points get a negative information measure calculated as
    # P/1-P log P
    # where P is the probability that a projection point maps - i.e. is not missing
    P = (len(points) - missing_count) / float(len(points))
    for index in range(len(projection)):
        if projection[index] == 0:
            projection[index] =  (P / float(1-P) ) * math.log(P,2.0)
        else:
            projection[index] = -1.0 * math.log(projection[index] / float(distribution.point_count), 2.0)

    return projection


def p_get_unsigned_information_projection((distribution, points)):
    projection = len(points) * [None]
    missing_count = 0
    index = 0
    for point in points:
        projection[index] = distribution.get_frequency(point,0)
        if projection[index] == 0:
            missing_count += 1
        index += 1

    missing_count = float(missing_count)

    if missing_count == len(points):
        raise prbdfException("no point maps - unable to calculate projection")

    # all points projected onto missing get frequency of approximate_zero_frequency (e.g. often .5) 
    # (but do not affect the projection of other points)
    for index in range(len(projection)):
        if projection[index] == 0:
            projection[index] = -1.0 * math.log(distribution.approximate_zero_frequency / float(distribution.point_count), 2.0)
        else:
            projection[index] = -1.0 * math.log(projection[index] / float(distribution.point_count), 2.0)

    return projection

def p_load(filename):
    preader = open(filename, "rb")
    pinstance = pickle.load(preader)
    preader.close()

    if not isinstance(pinstance, Distribution):
        raise prbdfException("%s is not a distribution object"%filename )
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


def default_count_provider(point):
    return 1

def default_value_provider(point):
    return point


def get_file_type(file_name):
    """
    infer  file type from file name for some commons file formats.
    """
    if re.search("(\.fasta|\.fa|\.seq)(\.|$)", file_name, re.I) != None:
        return "fasta"
    elif re.search("(\.fastq|\.fq)(\.|$)", file_name, re.I) != None:
        return "fastq"
    else:
        return "?"
    

def get_text_stream(file_name):
    """
    this is used by most of the value providers below that use text files
    (tab seperated, csv etc)
    """
    text_stream = None
    if re.search("\.gz$", file_name, re.I) != None:
        text_stream = gzip.open(file_name, 'rb')
    else:
        text_stream = open(file_name,"r")
    
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
    element of the tuple
    """
    if len(xargs) == 0:
        return (tuple(re.split("\t", record.strip())) for record in get_text_stream(file_name))
    else:
        return (tuple([outer_list(re.split("\t", record.strip()))[index] for index in xargs])  for record in get_text_stream(file_name))
    
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
        






#################################################
# build methods 
#################################################


def build_part((distob, slice_number)):
    print "build_part is building part %d"%slice_number
    return distob.get_partial_distribution(slice_number)

def build(distob, use="multithreads"):

    if use == "cluster":
        args = [(distob, slice_number)  for slice_number in range(0,distob.part_count)]
        #print "DEBUG : " + str(args)
        distob.part_dict = dict(pool.map(build_part,args))

        #print distob.part_dict
        return distob.get_distribution()
        
    elif use == "multithreads":
        pool = Pool(PROC_POOL_SIZE)
        
        args = [(distob, slice_number)  for slice_number in range(0,distob.part_count)]
        #print "DEBUG : " + str(args)

        print "mapping %s build parts to a pool of size %d"%(len(args), PROC_POOL_SIZE)
        #print "details of first part prbdf object : %s"%dir(args[0])
        distob.part_dict = dict(pool.map(build_part,args))

        #print distob.part_dict
        return distob.get_distribution()

    elif use == "singlethread":
        
        args = [(distob, slice_number)  for slice_number in range(0,distob.part_count)]
        #print "DEBUG : " + str(args)

        results = []
        for arg in args:
            results.append(build_part(arg))
            
        distob.part_dict = dict(results)

        #print distob.part_dict
        return distob.get_distribution()
    
    else:
        raise prbdfException("error - unknown resource specified for build : %s"%use)
    

        

