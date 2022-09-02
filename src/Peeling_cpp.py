import random
import sys
import getopt

from numpy import rint
from CountingMinSkecth import CountingMinSkecth
from LogScreen import LogScreen
from LogFile import LogFile
from GenericHashFunctionsSHA512 import GenericHashFunctionsSHA512


# Get the set of elements from the Universe (1 to maxVal) that may be included into the filter.
# m is the number of positions (counters) of each of the CMS arrays
# d is the number of arrays
# cms is the Count Min Sketch
# p is the P array with all the elements from the universe that returned positive from CBF
# sc is the printing element (LogScreen o LogFile) to print information
def peeling(m, d, cms, p, sc):
    # Retrieve the hash function used
    #hashf = cms.get_hashval()
    # T array with m positions to store the elements that are mapped to them
    # elements = [None] * m
    elements = [([None] * m) for i in range(d)]
    # Count of elements mapped to each position
    count = [([0] * m) for i in range(d)]

    # For all the positions in p
    for i in range(len(p)):
        # for the d hash functions
        for j in range(d):
            # Get the position mapped for the element p[i] and the jth hash function
            pos = cms.gethashval(p[i], j)
            # Retrieve the position pos of the T array            
            list_pos = elements[j][pos]
            # If no elements are assigned to that position, create a list and assign it
            if list_pos is None:
                list_pos = list()
                elements[j][pos] = list_pos
            # Include the element into the list of elements mapped to the position
            list_pos.append(p[i])
            # Increase the count of elements mapped to the position
            count[j][pos] += 1

    # Set that will store the positives that were extracted from the filter
    positives = set()

    # Values for the CMS counters
    #counters = cms.get_counters_array()
    
    while True:
        # If we found an element that could be extracted in this iteration
        found = False
        # Go through the m positions
        # TODO: exclude positions that were 0 in the previous iterations to speed up the process
        for i in range(d): 
            for j in range(m):
                # when the ith position does not have values, go to next position (para ccp darle la vuelta ai y j)
                if count[i][j] != 1:
                    continue
                # the position has values => it is not empty
                # we only want those positions where we have the same number of elements in T and CMS
                # if count[i][j] != counters[i][j]:
                #    continue
                # we found at least one position
                found = True
                # add the elements of the ith position to the positives
                positives.update(elements[i][j])

                # all these elements must be removed as well, but not from elements
                removers = elements[i][j].copy()

                # call the function that clears the removers and related false positives
                # pass True as last parameter as they are real positives
                clear_positions(elements, removers, cms, cms.getcountervalue(j, i), count, cms, d, sc, True)
            # if no new positives were found in the iteration, we should finish the algorithm

        if not found:
            break

    # return elements that were retrieved from the CBF
    return positives

# Function that clears the element from its positions in T and also clears all the related false positives
# elements is the T array
# positives is the list of elements to be removed
# count_cms is the list of counters from the CMS
# count_cms_decrement is the number of values to decrement
# count is the list of counters from T
# hashf is the hash function used in the CMS
# d is the number of arrays
# sc is the printing element (LogScreen o LogFile) to print information
# is_positive indicates if it is a real positive (true) or a false positive (false)
def clear_positions(elements, positives, count_cms, count_cms_decrement, count, cms, d, sc, is_positive):

    # Additional elements to be removed
    additional = list()
    # and iterate over them
    num = len(positives)
    # Traverse the positives list
    for i in range(num):
        # get next element to be removed
        next_positive = positives[i]
        # for the k hash functions
        for j in range(d):
            # Get the position mapped for the element and the jth hash function
            jpos = cms.gethashval(next_positive, j)
            # Element might have been removed in a different level of recursion
            if elements[j][jpos].count(next_positive) == 0:
                break
            # Remove the element from the position
            elements[j][jpos].remove(next_positive)
            # Reduce the T counter for that position
            count[j][jpos] -= 1
            # Reduce the CMS counter only when it is a real positive
            if is_positive:
                count_cms.decrementcountervalue(jpos, j, count_cms_decrement)
            # If no more elements are mapped to this position in the CMS
            # we can remove all the pending elements from T as they are false positives
            #if count_cms.getcountervalue(jpos, j) == 0 and count_cms.getcountervalue(jpos, j) != 0:
            if count_cms.getcountervalue(jpos, j) == 0 and count[j][jpos] != 0:
                # Add those elements to additional list
                additional.extend(elements[j][jpos])

    # Recursive call to remove the false positive elements
    if len(additional) > 0:
        # Pass False as last parameter as they are false positives
        clear_positions(elements, additional, count_cms, 0, count, cms, d, sc, False)

    return 0