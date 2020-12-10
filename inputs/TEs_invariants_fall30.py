#!/usr/bin/python

import os
import sys
import subprocess

te_file = open(sys.argv[1], 'r').readlines() # /Users/pogo/Documents/Maria/Baypass/Bergland/definitive_files_2019/Bergland_2014_for_geno.txt
new_file = open(sys.argv[2], 'w') # /Users/pogo/Documents/Maria/Baypass/Bergland/definitive_files_2019/Bergland_2014_for_geno_nofixed.txt

for t in te_file:
    line = (t.split("\n")[0]).split("\t")
    chrom = line[0]
    te = line[1]
    pos = line[2]
    alternative =[]
    reference = []
    # Change here according to the number of populations VERY IMPORTANT!!!!!
    for n in range(3,23):
        if n % 2 != 0:
            reference.append(float(line[n]))
        else:
            alternative.append(float(line[n]))
    if sum(alternative) != 0 and sum(reference) != 0:
        print>>new_file, "\t".join(t.split("\n"))
    else:
        print te
        print reference
        print alternative