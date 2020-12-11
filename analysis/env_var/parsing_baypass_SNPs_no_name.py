#!/usr/bin/python

# 31st May 2019
# Similar to /Users/pogo/Documents/Maria/Scripts/Baypass/parsing_baypass_results_new.py but for SNPs instead only FBti

# Do kind of the simmilar thing of the median but comparing how many of the 3 falls in the same category and put a 1,2 or 3 deppending on how many of them coincides and then order them to check how many 3, 2 or 1 are.

import sys
from numpy import median

seeds = int(sys.argv[2])*8

final_file = open(sys.argv[3], 'w')

with open(sys.argv[1]) as condition_file:
    for f in condition_file:
        c = '\t'.join(f.split())
        condition = c.split("\t")[0]
        snp = c.split("\t")[8]
        bfactor_list = []
        for x in range(4,seeds,8):
            bfactor_list.append(float(c.split("\t")[x]))
        bf_median = median(bfactor_list)
        print >> final_file, condition + "\t" + str(bf_median) + "\t" + snp



        # for x in range(0,seeds):
        #     flag = flag + 1
        #     c = '\t'.join(c1.split())
        #     cond = c.split("\t")[0]
        #     bfactor = float(c.split("\t")[5])
        #     TE_name_n = c.split("\t")[6]
        #     TE_name = TE_name_n.split("\n")[0]