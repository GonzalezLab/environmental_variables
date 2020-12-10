#!/usr/bin/python

# Simmilar to the way frequencies are calculated from Tlex
# Output will be a dataset to be used in BAYPASS
# Usage: python pool_frequencies.py Tresults name_output

import sys
import csv
import os


# In this case results will include the data from different datasets
results = list(csv.reader(open(sys.argv[1], 'r'), delimiter='\t'))
# We would need to separate them in columns by specific order (order stablished by the data order)
#file_result = sys.argv[2]
# File with annotation from TEs
annotation = list(csv.reader(open(sys.argv[2], 'r'), delimiter='\t'))

names = []
tes = []
freqs = []
strain = []
def_list = []


for r in results:
    if r[0] not in strain:
        strain.append(r[0])
    if r[1] not in tes:
        tes.append(r[1])

for t in tes:
    def_list = []
    coordenate = '';
    for r in results:
        if t == r[1]:
            NRL = int(r[12]) # 12 for Tlex2.3 and 11 for Tlex2
            NRR = int(r[18]) # 18 for Tlex2.3 and 17 for Tlex2
            NRA = int(r[5])
            if NRL+NRR > 0:
                if NRL != 0 and NRR == 0:
                    presence_reads = NRL
                    absence_reads = NRA
                    def_list.append(str(presence_reads))
                    def_list.append(str(absence_reads))
                elif NRL == 0 and NRR != 0:
                    presence_reads = NRR
                    absence_reads = NRA
                    def_list.append(str(presence_reads))
                    def_list.append(str(absence_reads))
                elif NRL != 0 and NRR != 0:
                    presence_reads = (int(NRL)+int(NRR))/2
                    absence_reads = NRA
                    def_list.append(str(presence_reads))
                    def_list.append(str(absence_reads))
            else:
                absence_reads = NRA
                presence_reads = 0
                def_list.append(str(presence_reads))
                def_list.append(str(absence_reads))
    for a in annotation:
        if a[0] == t:
            coordenate = a[2]
            chromosome = a[1]
    print str(chromosome) + "\t" + t + "\t" + str(coordenate) + "\t" +  "\t".join(def_list)
