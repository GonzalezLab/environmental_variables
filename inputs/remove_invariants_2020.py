
import sys


#77 3 120 1 218 3 52 3 63 0 107 6 28 0 98 13 105 3 53 12 100 17
#0  1   2 3 4   5 6 7   8 9 10  11 12 13 14 15 16 17 18 19 20 21

#'2L 5403 C G\t57 19 125 12 200 36 62 2 29 14 98 31 49 5 115 20 98 5 77 3 82 2\n'

input_file = open(sys.argv[1], 'r').readlines()
map_file = open(sys.argv[2], 'w')
geno_file = open(sys.argv[3], 'w')

for i in input_file:
    info = i.split("\t")[0]
    i_list = ((i.split("\n")[0]).split("\t")[1]).split(" ")
    if i_list[0] == "NA":
        continue
    else:
        present = int(i_list[0]) + int(i_list[2]) + int(i_list[4]) + int(i_list[6]) + int(i_list[8]) + int(i_list[10]) + int(i_list[12]) + int(i_list[14]) + int(i_list[16]) + int(i_list[18]) + int(i_list[20])
        absent = int(i_list[1]) + int(i_list[3]) + int(i_list[5]) + int(i_list[7]) + int(i_list[9]) + int(i_list[11]) + int(i_list[13]) + int(i_list[15]) + int(i_list[17]) + int(i_list[19]) + int(i_list[21])
    if present != 0 and absent != 0:
        print >> map_file, i.split("\n")[0]
        print >> geno_file, (i.split("\n")[0]).split("\t")[1]


map_file.close()
geno_file.close()