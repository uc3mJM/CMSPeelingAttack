from ast import Break
from pydoc import doc
import random
from re import M
import sys
import getopt

from numpy import ubyte
from CountingMinSkecth import CountingMinSkecth
from LogScreen import LogScreen
from LogFile import LogFile
from GenericHashFunctionsSHA512 import GenericHashFunctionsSHA512
import ipaddress
from Peeling_cpp import peeling
from cms.cms_module import *

# Main method
def main():

    # Default values for
    # Number of counters in the filter
    m = 4096
    # Number of arrays in the filter
    d = 3
    # Number of correct insertions to be done into the filter
    n = 2000
    # Size of the universe of elements -> (Cambiar por mascara)
    u = 10000
    # Step to increase the universe size
    # Try u, then u+step, u+2*step...
    step = 10000
    # bit width per word
    # k = 5
    # Hash function to be used (md5 by default)
    hash_f = 'md5'
    # Number of iterations
    it = 100
    # Maximum size of the universe
    end_u = -1

    # Retrieve 2the option values from command line
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hm:d:a:i:f:")
    except getopt.GetoptError:
        print('argv[0] -m <size> -d <array_number> -n <insertions> -u <universe_size> -s <step_size> -a <hash> -i <iterations>')
        sys.exit(2)

    for opt, arg in opts:
        # Help option. Print help and leave.
        if opt == '-h':
            print('argv[0] -m <size> -d <array_number> -n <insertions> -u <universe_size> -s <step_size> -a <hash> -i <iterations>')
            sys.exit()
        # -m option for setting the number of counters in the filter
        elif opt == "-m":
            m = int(arg)
        # -k option to set the number of hash elements to select the bits to be set
        elif opt == "-d":
            d = int(arg)
        # -n option for setting the number of correct elements that will be inserted into the filter
        elif opt == "-n":
            n = int(arg)
        # -u option for setting the size of the universe
        elif opt == "-u":
            u = int(arg)
        # -s option for setting the step to increase the size
        elif opt == "-s":
            step = int(arg)
        # -e option for setting the step to increase the size
        elif opt == "-e":
            end_u = int(arg)
        # -a option to change the default md5 hash to other ("sha512" supported)
        elif opt == "-a":
            hash_f = arg
        # -i option to change the number of iterations
        elif opt == "-i":
            it = int(arg)

    # Pass the parameters to the run function
    #run(m, mask, step, it, mask_end, d, hash_f,f)
    # Pass the parameters to the run function
    run(m, d, n, u, step, end_u, it, hash_f)

    return

# Function to find all elements from the universe that returns a estimate value bigger than 1 from CMS
# cbf is the Count Min Sketch
# max_val is the maximum integer value. Universe will include elements from 1 to max_val
def find_p(cms, max_val):
    # Create the list P of (true and false) positive elements
    p = list()
    # Check all elements of the universe, from 1 to max_val4
    j= 0
    for i in range(max_val+1):
        # If one of the positions is 0, then it is a negative
        # Otherwise, add it to P
        if cms.estimate(i):
            p.append(i)
            j = j +1

    return p

def generate_random_elements(num, cms=None, ds=None, max_val=1000000000, exclude=None):
    # Elements are added to the set to check for repetitions
    if exclude is None:
        exclude = set()
    s = set()
    s.clear()

    # Keeps data of stored elements
    stored = 0
    # Generate elements until the value "stored" is reached
    while stored < num:
        # Generate integers between 1 and max_val
        entry = random.randint(1, max_val)
        # if entry was already selected or is in the exclude set,
        # go to next iteration
        if entry in s or entry in exclude:
            continue
        # When an CountingMinSketch is received
        if cms is not None:
            # Add the entry to the filter
            cms.update(entry)
        # When a list is received
        if ds is not None:
            # Add the element to the list
            ds.append(entry)
        # Another element has been stored
        stored = stored + 1
        # Add it to the set so they are not repeated
        s.add(entry)

    return stored

def increase_random_elements(num, cms, ds, stored):

    # Crear una lista con  
    random_values = [] 
    mu = num/2
    sigma = num/8

    inc = num * 10

    while len(random_values) < inc:
        sample = random.gauss(mu, sigma)
        if sample >= 0 and sample < num:
            random_values.append(int(sample))

    for values in random_values:
        cms.add(ds[values])
        stored += 1
    return stored


# Run the actual experiment using the parameters received
def run(m=65536, d=3, n=10000, u=240000, step=30000, u_end=-1, iters=10, hash_f='md5'):
    
    max_val = u
    sc = LogScreen()
    i = 0

    # Directory to write the logs
    directory = './logs/'

    # Definition of the name of the output files.
    log_output = 'result_m%s_d%s_n%s_h%s_i%s' % (m, d, n, hash_f, iters)
    log = LogFile(directory + log_output, "w")

    # CountMinSketch file
    cms = None

    # Number of times to execute each experiment to get an average
    total_iterations = iters

    # Message printing the parameters used for the experiment
    info = "Initializing parameters: Arrays=%d, Counters=%d, hash_f=%s, Insertions=%d" % (
        d, m, hash_f, n)
    sc.write(info)

    # Build the CMS using CPP module
    cms = CMSF(m, d) # type: ignore

    # list to store the true positive elements to check if they were retrieved correctly
    ds = list()

    info = "universe;false_pos;total_inserted;total_stored;success_iter_perc;average_retrieved;worst_retrieved;iterations"
    log.write(info + "\n")

    # Continue the execution until the number of FP does not allow any iteration to retrieve the content
    while True:

        # Message printing the parameters used for the experiment
        info = "Peeling test with Universe size = %d" % (max_val)
        sc.write(info)

        # Number of total elements from P (true and false positives)
        total_elements = 0
        # Number of total elements introduced (true positives)
        total_inserted = 0
        # Number of total elements stored
        total_stored = 0
        # Number of iterations completed successfully
        total_completed = 0
        # Number of total elements retrieved from the CMS
        total_extracted = 0
        # Number of elements retrieved in the worst case iteration
        worst_extracted = n

        # Run for the number of iterations expected
        for i in range(total_iterations):

            # Clear the CMS and ds list for next iteration
            cms.clear()
            ds.clear()

            # Generate a set of n elements between 1 and max_val.
            # Include the values in CMS and add them to ds.
            
            stored = generate_random_elements(n, cms, ds, max_val)
            #stored = increase_random_elements(n, cms, ds, stored)
            
            # Find all elements of the universe 1 to max_val that gets a positive result from CBF
            p = cms.findp(max_val)
     
            # Accumulate the elements in total_elements for later processing.
            total_elements += len(p)
            total_inserted += len(ds)
            total_stored += stored
            
            # Retrieve all the elements from CMS that you are capable of
            #info = "Peeling algoritm..."
            # sc.write(info)
            positives = peeling(m, d, cms, p, sc)

            # If ds and positives have the same size and there is no difference between them
            if len(positives) == len(ds) and len(positives.difference(ds)) == 0:
                # The iteration was completed successfully
                total_completed = total_completed + 1
            elif len(positives.difference(ds)) > 0:
                sc.write("Error: false positives extracted")
                sc.write(positives.difference(ds))
                sc.write(ds)
                exit(1)
            # Accumulate the number of elements retrieved in total_extracted for later processing.
            total_extracted += len(positives)
            # Check for the worst case iteration and store the number of elements extracted
            if len(positives) < worst_extracted:
                worst_extracted = len(positives)
            # peeling(m, k, cbf, ds, sc)

        # Update the average of elements extracted by dividing the accumulated sum by the number of iterations
        total_extracted = (total_extracted / total_iterations)
        # Number of false positives by getting the average of elements from P and subtracting the real positives
        false_positives = (total_elements - total_inserted) / total_iterations

        sc.write("False positives: %d and success in %d %% iterations" % (false_positives,
                                                                          (total_completed / total_iterations * 100)))
        sc.write("Average extracted: %d; Worst case: %d" % (total_extracted,
                                                            worst_extracted))
        # Store the information into a file
        info = "%s;%d;%d;%d;%s;%s;%s;%s" % (max_val, false_positives, total_inserted, total_stored, (total_completed / total_iterations * 100),
                                      total_extracted, worst_extracted, total_iterations)
        log.write(info + "\n")
        log.flush()

        # If total_extracted is less than 1%, finish the process (at least one run was executed)
        if total_extracted < 0.01 * n:
            break

        max_val += step

        if u < u_end < max_val:
            break

    return

if __name__ == "__main__":
    main()
