import sys


xtx_file = open(sys.argv[1], 'r').readlines()
var_file = open(sys.argv[2], 'r').readlines()
final_file = open(sys.argv[3], 'w')

for x in xtx_file:
    v_list = (x.split("\n")[0]).split("\t")
    fbgn = v_list[0]
    chrom = v_list[1]
    pos = v_list[2]
    score = v_list[3]
    qvalue = v_list[4]
    snp_name = v_list[5]
    flag = "NO"
    for v in var_file:
        var_value = v.split("\n")[0]
        if var_value == snp_name:
            flag = "YES"

    print >> final_file, str(fbgn) + "\t" + str(chrom) + "\t" + str(pos) + "\t" + str(score) + "\t" + str(qvalue) + "\t" + str(flag)


final_file.close()
