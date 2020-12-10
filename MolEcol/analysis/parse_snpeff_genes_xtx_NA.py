import sys


snpeff_file = open(sys.argv[1], 'r').readlines()
xtx_file = open(sys.argv[2], 'r').readlines()
var_file = open(sys.argv[3], 'r').readlines()
final_file = open(sys.argv[4], 'w')

for s in snpeff_file:
    snp = (s.split("\n")[0]).split("\t")[3]
    fbgn = s.split("\t")[0]
    flag = "NO"
    def_chrom = ""
    def_pos = ""
    def_score = ""
    def_qvalue = ""
    for v in xtx_file:
        v_list = (v.split("\n")[0]).split("\t")
        chrom = v_list[0]
        pos = v_list[1]
        score = v_list[2]
        qvalue = v_list[3]
        name = v_list[4]
        if snp == name:
            def_chrom = chrom
            def_pos = pos
            def_score = score
            def_qvalue = qvalue
            for a in var_file:
                snp_name = a.split("\n")[0]
                if snp == snp_name:
                    flag = "YES"
    print >> final_file, str(fbgn) + "\t" + str(def_chrom) + "\t" + str(def_pos) + "\t" + str(def_score) + "\t" + str(def_qvalue) + "\t" + str(flag)


final_file.close()
