#!/usr/bin/python

# 16th March 2018
# Creates geno file needed as input for running Baypass
# Using VCF file from Bergland 2014

import os
import sys
import subprocess

vcf = open("/Users/pogo/Documents/Maria/Baypass/DrosEU/data/DrosEU-PoolSnp_full-filtered_ann_exmulti_nohash.vcf", 'r').readlines()

geno_file = []
#columns = [9,11,13,15,16,17,18,21,23,25,29,31,32,33,36,37,41,42,44,47,49,50,56] # Selction only of populations we are going to use for spring
#columns = [10,12,14,19,20,22,28,34,35,38,39,40,46,48,51,52,53,54,55] # Selction only of populations we are going to use for autum
columns = [19,20,22,28,34,35,38,39,40,46] # Selction only of populations we are going to use for fall for only the 30 samples
#columns = [9,11,13,15,16,17,18,19,20,21,23,25,29,31,32,33,34,35,36,37,39,41,42,44,45,47,49,50,52,53,55,56] # For running all DrosEU
#columns = [18,19,20,21,23,25,29,31,32,33,34,35,36,37,39,41,42,44,45,47] # For running all DrosEU only for the 30 strains

flag = 0
geno_output = open("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/Droseu14_fall30_exmulti.geno", 'w')
for line in vcf:
    geno_file_line = []
    chromosome = line.split("\t")[0]
    position = line.split("\t")[1]
    geno_file_line.append(chromosome)
    geno_file_line.append(position)
    alternative_list = []
    reference_list = []
    for x in columns:
        population = line.split("\t")[x]
        alternative = population.split(":")[2]
        original = population.split(":")[1]
        #original = int(total) - int(alternative)
        geno_file_line.append(str(original))
        geno_file_line.append(str(alternative))
        if original == "." or alternative == ".":
            reference_list.append(0)
            alternative_list.append(0)
        else:    
            reference_list.append(int(original))
            alternative_list.append(int(alternative))
    ref_total = sum(reference_list)
    alt_total = sum(alternative_list)
    if ref_total != 0 and alt_total != 0:
        geno_file.append(geno_file_line)
        print>> geno_output, "\t".join(geno_file_line)
    flag = flag + 1
    print(flag)

print("longitud: " + str(len(geno_file)))

# geno_output = open("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/Droseu14_all30_exmulti.geno", 'w')
# for geno in geno_file:
#     print >> geno_output, "\t".join(geno)