# -*- coding: utf-8 -*-
"""
Created on Tue Jun 06 15:17:57 2021

@author: markus
"""
import time
import gengraph
from contextlib import contextmanager
import sys, os

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

a=time.time()
with suppress_stdout():

    testfile1 = '.\\alignments\\aln_south_africa_subset.fa'
    testfile2 = '.\\alignments\\aln_subset.fa'
    testfile3 = '.\\alignments\\temp_darren.fa'

    new_subgraph,Testlist = gengraph.fasta_alignment_to_subnet(testfile1)

print (sys.getsizeof(Testlist))
print (sys.getsizeof(new_subgraph))
print((time.time()-a))



# testfile1 suppress_stdout() on: 76s
# testfile1 suppress_stdout() off: 398s
