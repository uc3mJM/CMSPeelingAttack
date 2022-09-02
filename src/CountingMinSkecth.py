import hashlib
import math
import random
import string
from typing import List
import timeit


from LogScreen import LogScreen
from GenericHashFunctionsMD5 import GenericHashFunctionsMD5


# Adaptive CountMinSketch filter
class CountingMinSkecth:
    #cms_structure: list[int]

    def __init__(self, m=65536, d=5, hash_f=None):
        # number of counters
        self.m = m
        # number of arrays
        self.d = d

        # the structure is stored as a series of arrays
        self.cms_structure_arrays = [([0] * m) for i in range(d)]

        # the hash class used to generate the functions
        if hash_f is None:
            self.hash = GenericHashFunctionsMD5(m, d)
        else:
            self.hash = hash_f
        # Depth o the Count Min Sketch, this is the number of counters arrays, each one wiht a unique hash function
        

    # clear the list of counters
    def clear(self):
        self.cms_structure_arrays = [([0] * len(self.cms_structure_arrays[0])) for i in range(len(self.cms_structure_arrays))]

    # Change the hash object that generates the function
    def set_hash(self, hash_object):
        if hash_object is not None:
            self.hash = hash_object
        return

    # Get the hash object that generates the function
    def get_hash(self):
        return self.hash

    # method to add an element into the filter
    def add(self, data):
        # extract a position from each hash to set the bit in the selected word
        for i in range(self.d):
            # position for the ith hash of the group_i function 
            # the final position in the array is a combination
            # of word index and bit index
            idx = self.hash.getbit_idx(data, i)
            # idx = int(bitidx)

            # set the appropriate bit
            self.cms_structure_arrays[i][idx] += 1

        return

    # method to delete an element from the filter
    def remove(self, data):
        # extract a position from each hash to set the bit in the selected word
        for i in range(self.d):
            # position for the ith hash of the group_i function 
            # the final position in the array is a combination
            # of word index and bit index
            idx = self.hash.getbit_idx(data, i)
            # idx = int(bitidx)

            # set the appropriate bit
            self.cms_structure_arrays[i][idx] -= 1

        return

    # check the cms filter for the specified data
    def check(self, data, threshold=1):
        # extract a position from each hash to get the bit from the selected word
        for i in range(self.d):
            # position for the ith hash of the group_i function
            # the final position in the array is a combination
            # of word index and bit index
            idx = self.hash.getbit_idx(data, i)

            # check if the bit is at least the threshold for the specified word
            # if not, the data is a negative
            if self.cms_structure_arrays[i][idx] < threshold:
                return False
        return True

    def estimate(self, data):
        estimate_value = 2147483647
        # extract a position from each hash to get the bit from the selected word
        #estimate_values = [self.cms_structure_arrays[i][self.hash.getbit_idx(data, i)] for i in range(self.d)]

        
        for i in range(self.d):
            # position for the ith hash of the group_i function
            # the final position in the array is a combination
            # of word index and bit index
            idx = self.hash.getbit_idx(data, i)

            # check if the bit is at least the threshold for the specified word
            # if not, the data is a negative
            if self.cms_structure_arrays[i][idx] < estimate_value:
                estimate_value = self.cms_structure_arrays[i][idx]

        return estimate_value

    # Retrieve the value of a counter
    def get_counter_array(self, position, d):
        if len(self.cms_structure_arrays[d]) < position:
            return 0
        return self.cms_structure_arrays[d][position]

    # Retrieve the structure of counters
    def get_counters_array(self):
        return self.cms_structure_arrays
        
    # Debug function for printing content
    def printme(self):
        sc = LogScreen()
        for i in range(self.d):
            for j in range(self.m):
                info = "Row=%d, Array=%d, Count=%d" % (i, self.cms_structure_arrays[i][j])
                sc.write(info)
        return
