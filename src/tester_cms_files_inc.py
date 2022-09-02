from ast import Break
import random
from re import M
import sys
import getopt
from numba import jit

from numpy import ubyte
from CountingMinSkecth import CountingMinSkecth
from LogScreen import LogScreen
from LogFile import LogFile
from GenericHashFunctionsSHA512 import GenericHashFunctionsSHA512
import ipaddress
from Peeling import peeling
import timeit

# Main method
def main():

    # Default values for
    # Number of counters in the filter
    m = 2048
    # Number of arrays in the filter
    d = 4
    # Number of correct insertions to be done into the filter
    n = 150000
    # Size of the universe of elements -> (Cambiar por mascara)
    #u = 240000
    # Mask for the IPv4 netwwork to use as an universe of elements.For example -> /16, equivalent to u = 65536
    mask = 12
    # Step to increase the universe size, decreasing the mask value
    # Try mask, then mask - step, mask - 2*step...
    step = 5000
    # bit width per word
    # k = 5
    # Hash function to be used (md5 by default)
    hash_f = 'md5'
    # Number of iterations
    it = 1
    # Maximum size of the insertions
    end_n = -1
    # file path to read
    #f=[r".\traces\equinix-chicago.dirB.20140619-132600.txt",r".\traces\equinix-chicago.dirA.20140619-130900.txt",r".\traces\equinix-sanjose.dirA.20140320-130400.txt"]
    f = r".\traces\equinix-sanjose.dirA.20140320-130400.txt"
    # Retrieve the option values from command line
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hm:d:a:u")
    except getopt.GetoptError:
        print('argv[0] -m <size> -d <array_number> -n <insertions> -a <hash> -mask <network_size> -step <step_size> -end_n <max_insertions> -i <iterations>')
        sys.exit(2)

    for opt, arg in opts:
        # Help option. Print help and leave.
        if opt == '-h':
            print('argv[0] -m <size> -d <array_number> -n <insertions> -a <hash> -mask <network_size> -step <step_size> -end_n <max_insertions> -i <iterations>')
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
        # -a option to change the default md5 hash to other ("sha512" supported)
        elif opt == "-a":
            hash_f = arg
        # -s option for setting the step to increase the size
        elif opt == "-s":
            step = int(arg)
        # -e option for setting the step to increase the size
        elif opt == "-end_n":
            end_n = int(arg)
        # -mask option for setting the first size of the network (universe size)
        elif opt == "-mask":
            mask = int(arg)
        # -i option to change the number of iterations
        elif opt == "-i":
            it = int(arg)

    # Pass the parameters to the run function
    run(m, d, n, mask, step, end_n, it, hash_f, f)

    return

def check_biggest_subnets(path_to_file,mask):
    
    nets = {}

    with open(path_to_file) as f:
        lines = f.readlines()
        for line in lines:
            # Get the source IP
            source_ip = line.split(' ')[0]
            if source_ip != '':
                net = ipaddress.ip_network(source_ip + '/' + str(mask), strict=False)
                if net in nets:
                    nets[net] += 1
                else:
                    nets[net] = 1

    # Printing all element with its count of Occurrence
    #for key,value in nets.items():
    #    print("The Occurrence of {0} in the sample list is: {1}".format(key, value))
    
    max_key = max(nets, key=nets.get)
    print("The Occurrence of {0} in the sample list is: {1}".format(max_key, nets[max_key]))

    return max_key


def find_p_ipaddress(cms, net):

    # Create the list P of (true and false) positive elements
    p = list()
    # Check all elements of the universe, from 1 to max_val
    # for ip in ipaddress.IPv4Network(net + '/' + mask):
    #    print(ip)
    
    #addr = [int(ip) for ip in ipaddress.IPv4Network(net)]
    p = [p.append(int(ip)) for ip in ipaddress.IPv4Network(net) if cms.estimate(int(ip)) > 0 ]
    """
    for ip in addr:
        if cms.estimate(ip) > 0:
            p.append(ip)
    """
    #return p



def read_CAIDA(path_to_file, net, mask, max_n, exclude=None):
    # Elements are added to the set to check for repetitions

    if exclude is None:
        exclude = set()
    s = set()
    s.clear()
    
    # IP to be introduced in the CMS
    ips = [] 

    # Keeps data of stored elements
    stored = 0

    with open(path_to_file) as f:
        lines = f.readlines()

        for line in lines:
            source_ip = line.split(' ')[0]
            if source_ip != '':

                if source_ip in exclude:
                    continue

                ip = ipaddress.ip_network(source_ip + '/' + str(mask) , strict=False)
                # IF the IP is in the network stored it
                if net == ip:

                    ips.append(int(source_ip))
                    stored += 1
            # Generate elements until the value "stored" is reached
            #if stored >= max_n:
            #    break

    return ips

def initialize_values(ips, max_n, cms=None, ds=None, exclude=None):
    # Elements are added to the set to check for repetitions

    if exclude is None:
        exclude = set()
    s = set()
    s.clear()
    
    # Keeps data of stored elements
    stored = 0

    random.shuffle(ips)
    for source_ip in ips:
        # go to next iteration
        # When an CountingMinSketch is received
        if cms is not None:
            # Add the entry to the filter
            cms.add(source_ip)
        # Another element has been stored
        stored = stored + 1
        # is value already in the cms
        if source_ip in s:
            continue
        # When a list is received
        if ds is not None:
            # Add the element to the list
            ds.append(source_ip)
        # Add it to the set so they are not repeated
        s.add(source_ip)

        # Generate elements until the value "stored" is reached
        if stored >= max_n:
            break

    return stored


# Run the actual experiment using the parameters received
def run(m=65536, d=6, n=100, mask = 16, step = 100, n_end=-1, iters=10, hash_f='md5', f=None):
    
    max_n = n
    sc = LogScreen()
    i = 0

    #universe size
    u = pow(2, (32-mask))
    # Number of times to execute each experiment to get an average
    total_iterations = iters

    # Directory to write the logs
    directory = './logs/'

    # Definition of the name of the output files.
    log_output = 'result_m%s_d%s_u%s_h%s_i%s_f' % (m, d, u, hash_f, total_iterations)
    log = LogFile(directory + log_output, "w")

    # CountMinSketch file
    cms = None

    # Message printing the parameters used for the experiment
    info = "Initializing parameters: Arrays=%d, Counters=%d, hash_f=%s, Universe size= %d" % (
        d, m, hash_f, u)
    sc.write(info)

    # Build the filter passing a SHA512 hash function
    if hash_f == 'sha512':
        sha = GenericHashFunctionsSHA512(k=m, nhash=d)
        cms = CountingMinSkecth(m=m, d=d, hash_f=sha)
        # Otherwise build it using the default MD5 hash
    else:
        cms = CountingMinSkecth(m=m, d=d)

    # list to store the true positive elements to check if they were retrieved correctly
    ds = list()

    info = "universe;false_pos;total_inserted;total_stored;success_iter_perc;average_retrieved;worst_retrieved;iterations"
    log.write(info + "\n")

    # Message printing the parameters used for the experiment
    info = "Initializing CMS with file \"%s\"" % (f)
    sc.write(info)

    # Message printingt
    info = "Checking the subnets with more adresses"
    sc.write(info)

    #net = check_biggest_subnets(f,mask)  
    #net = ipaddress.ip_network("68.142.0.0/16", strict=False)
    #net = ipaddress.ip_network("68.142.0.0/15", strict=False)
    net = ipaddress.ip_network("68.128.0.0/12", strict=False)


    # Continue the execution until the number of FP does not allow any iteration to retrieve the content
    while True:

        #universe size
        u = pow(2, (32-mask))
        i = 0

        # Message printing the parameters used for the experiment
        #info = "Peeling test with Mask = %d, Universe size = %d" % (mask, u)
        #sc.write(info)
        
        # Message printing the parameters used for the experiment
        info = "Peeling test with %d insertions" % (max_n)
        sc.write(info)

        # Number of total elements from P (true and false positives)
        total_elements = 0
        # Number of total unique elements introduced (true positives)
        total_inserted = 0
        # Number of total elements stored
        total_stored = 0
        # Number of iterations completed successfully
        total_completed = 0
        # Number of total elements retrieved from the CMS
        total_extracted = 0
        # Number of elements retrieved in the worst case iteration
        worst_extracted = u

        ips = read_CAIDA(f, net, mask, max_n)

        # Run for the number of iterations expected
        for i in range(total_iterations):

            # Clear the CMS and ds list for next iteration
            cms.clear()
            ds.clear()

            stored = initialize_values(ips, max_n, cms, ds)
            #Check the subnet with more adresses
            # net = check_biggest_subnets(f,mask)0
            # Find all elements of the universe 1 to max_val that gets a positive result from CBF
            p = find_p_ipaddress(cms, net)

            # Accumulate the elements in total_elements for later processing.
            total_elements += len(p)
            total_inserted += len(ds)
            total_stored += stored

            # Retrieve all the elements from CMS that you are capable of
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
            i += 1

        # Update the average of elements extracted by dividing the accumulated sum by the number of iterations
        total_extracted = (total_extracted / total_iterations)
        # Number of false positives by getting the average of elements from P and subtracting the real positives

        false_positives = (total_elements - total_inserted) / total_iterations

        sc.write("False positives: %d and success in %d %% iterations" % (false_positives,
                                                                          (total_completed / total_iterations * 100)))
        sc.write("Average extracted: %d; Worst case: %d" % (total_extracted,
                                                            worst_extracted))
        # Store the information into a file
        info = "%s;%d;%d;%d;%s;%s;%s;%s" % (u, false_positives, total_inserted, total_stored, (total_completed / total_iterations * 100),
                                      total_extracted, worst_extracted, total_iterations)
        log.write(info + "\n")
        log.flush()

        max_n += step

        if n < n_end < max_n:
            break

    return

if __name__ == "__main__":
    main()
