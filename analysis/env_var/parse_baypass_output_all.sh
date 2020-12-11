#### Created 24th September 2019

# Local in: 
for i in $(seq 50 $END); do awk '{print $1"\t"$2}' all30_droseu14_noinv_auto.geno.map.${i} > all30_droseu14_auto.geno.names.${i}; done
rsync

# Interactive mode
for f in $(seq 50 $END); do for i in $(seq 68 $END); do awk -v var=$i '$1 == var {print $0}' all30IS_${f}_summary_betai_reg.out  > ${i}_all30IS_${f}_betai_reg.out; done; done
for f in $(seq 50 $END); do for i in $(seq 68 $END); do awk -v var=$i '$1 == var {print $0}' all30ISs1500_${f}_summary_betai_reg.out  > ${i}_all30ISs1500_${f}_betai_reg.out; done; done
for f in $(seq 50 $END); do for i in $(seq 68 $END); do awk -v var=$i '$1 == var {print $0}' all30ISs567_${f}_summary_betai_reg.out  > ${i}_all30ISs567_${f}_betai_reg.out; done; done

for f in $(seq 50 $END); do awk -v var=69 '$1 == var {print $0}' all30IS_${f}_summary_betai_reg.out  > 69_all30IS_${f}_betai_reg.out; done
for f in $(seq 50 $END); do awk -v var=69 '$1 == var {print $0}' all30ISs1500_${f}_summary_betai_reg.out  > 69_all30ISs1500_${f}_betai_reg.out; done
for f in $(seq 50 $END); do awk -v var=69 '$1 == var {print $0}' all30ISs567_${f}_summary_betai_reg.out  > 69_all30ISs567_${f}_betai_reg.out; done

# Interactive ran out of time:

for f in {45..50}; do for i in $(seq 68 $END); do awk -v var=$i '$1 == var {print $0}' all30ISs1500_${f}_summary_betai_reg.out  > ${i}_all30ISs1500_${f}_betai_reg.out; done; done
for f in {45..50}; do for i in $(seq 68 $END); do awk -v var=$i '$1 == var {print $0}' all30ISs567_${f}_summary_betai_reg.out  > ${i}_all30ISs567_${f}_betai_reg.out; done; done

# Add coordinates to each SNP
for f in $(seq 50 $END); do for i in $(seq 69 $END); do paste /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/all30/names/all30_droseu14_auto.geno.names.${f} ${i}_all30IS_${f}_betai_reg.out > ${i}_all30IS_${f}_betai_names.out; done; done
for f in $(seq 50 $END); do for i in $(seq 69 $END); do paste /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/all30/names/all30_droseu14_auto.geno.names.${f} ${i}_all30ISs1500_${f}_betai_reg.out > ${i}_all30ISs1500_${f}_betai_names.out; done; done
for f in $(seq 50 $END); do for i in $(seq 69 $END); do paste /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/all30/names/all30_droseu14_auto.geno.names.${f} ${i}_all30ISs567_${f}_betai_reg.out > ${i}_all30ISs567_${f}_betai_names.out; done; done

# Paste three seeds output
for f in $(seq 50 $END); do for i in $(seq 69 $END); do paste ${i}_all30IS_${f}_betai_names.out ${i}_all30ISs1500_${f}_betai_names.out ${i}_all30ISs567_${f}_betai_names.out > ${i}_all30IS_${f}_observed_names.out; done; done

sbatch median_script.sh

### On 25th September 2019
mv *_median_observed_names.out median_observed_names/
# In: /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/all30/median_observed_names
for i in $(seq 68 $END); do cat ${i}_all30IS_*_median_observed_names.out > ${i}_all30IS_observed_names.out; done

for i in $(seq 68 $END); do awk '$3=="2L" {print $0}' ${i}_all30IS_observed_names.out > 2L_${i}_all30IS_observed_names.out; done
for i in $(seq 68 $END); do awk '$3=="2R" {print $0}' ${i}_all30IS_observed_names.out > 2R_${i}_all30IS_observed_names.out; done
for i in $(seq 68 $END); do awk '$3=="3R" {print $0}' ${i}_all30IS_observed_names.out > 3R_${i}_all30IS_observed_names.out; done
for i in $(seq 68 $END); do awk '$3=="3L" {print $0}' ${i}_all30IS_observed_names.out > 3L_${i}_all30IS_observed_names.out; done

for i in $(seq 68 $END); do grep "FBti" ${i}_all30IS_observed_names.out > FBti_${i}_all30IS_observed_names.out; done

for i in $(seq 68 $END); do awk -F '\t' 'NR==FNR {id[$1]; next} $3 in id' /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/all30/TEs/2L_droseu14_all30_noinv_TEs_c FBti_${i}_all30IS_observed_names.out > 2L_FBti_${i}_all30IS_observed_names.out; done
for i in $(seq 68 $END); do awk -F '\t' 'NR==FNR {id[$1]; next} $3 in id' /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/all30/TEs/2R_droseu14_all30_noinv_TEs_c FBti_${i}_all30IS_observed_names.out > 2R_FBti_${i}_all30IS_observed_names.out; done
for i in $(seq 68 $END); do awk -F '\t' 'NR==FNR {id[$1]; next} $3 in id' /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/all30/TEs/3R_droseu14_all30_noinv_TEs_c FBti_${i}_all30IS_observed_names.out > 3R_FBti_${i}_all30IS_observed_names.out; done
for i in $(seq 68 $END); do awk -F '\t' 'NR==FNR {id[$1]; next} $3 in id' /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/all30/TEs/3L_droseu14_all30_noinv_TEs_c FBti_${i}_all30IS_observed_names.out > 3L_FBti_${i}_all30IS_observed_names.out; done

for i in $(seq 68 $END); do cat 2L_${i}_all30IS_observed_names.out 2L_FBti_${i}_all30IS_observed_names.out > 2L_FBti_${i}_all30IS_observed_names_unsorted.out; done
for i in $(seq 68 $END); do cat 2R_${i}_all30IS_observed_names.out 2R_FBti_${i}_all30IS_observed_names.out > 2R_FBti_${i}_all30IS_observed_names_unsorted.out; done
for i in $(seq 68 $END); do cat 3R_${i}_all30IS_observed_names.out 3R_FBti_${i}_all30IS_observed_names.out > 3R_FBti_${i}_all30IS_observed_names_unsorted.out; done
for i in $(seq 68 $END); do cat 3L_${i}_all30IS_observed_names.out 3L_FBti_${i}_all30IS_observed_names.out > 3L_FBti_${i}_all30IS_observed_names_unsorted.out; done

for i in $(seq 68 $END); do sort -n -k4,4 2L_FBti_${i}_all30IS_observed_names_unsorted.out > 2L_FBti_${i}_all30IS_observed_names_sorted.out; done
for i in $(seq 68 $END); do sort -n -k4,4 2R_FBti_${i}_all30IS_observed_names_unsorted.out > 2R_FBti_${i}_all30IS_observed_names_sorted.out; done
for i in $(seq 68 $END); do sort -n -k4,4 3R_FBti_${i}_all30IS_observed_names_unsorted.out > 3R_FBti_${i}_all30IS_observed_names_sorted.out; done
for i in $(seq 68 $END); do sort -n -k4,4 3L_FBti_${i}_all30IS_observed_names_unsorted.out > 3L_FBti_${i}_all30IS_observed_names_sorted.out; done

for i in $(seq 68 $END); do cat 2L_FBti_${i}_all30IS_observed_names_sorted.out 2R_FBti_${i}_all30IS_observed_names_sorted.out 3L_FBti_${i}_all30IS_observed_names_sorted.out 3R_FBti_${i}_all30IS_observed_names_sorted.out > ${i}_all30IS_observed_complete.out; done

# On 26th September
# Missing variable 69
for f in $(seq 50 $END); do python /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/parsing_baypass_SNPs_with_name.py 69_all30IS_${f}_observed_names.out 3 69_all30IS_${f}_median_observed_names.out; done
mv 69_all30IS_*_median_observed_names.out median_observed_names/
# In: /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/all30/median_observed_names
cat 69_all30IS_*_median_observed_names.out > 69_all30IS_observed_names.out

awk '$3=="2L" {print $0}' 69_all30IS_observed_names.out > 2L_69_all30IS_observed_names.out
awk '$3=="2R" {print $0}' 69_all30IS_observed_names.out > 2R_69_all30IS_observed_names.out
awk '$3=="3R" {print $0}' 69_all30IS_observed_names.out > 3R_69_all30IS_observed_names.out
awk '$3=="3L" {print $0}' 69_all30IS_observed_names.out > 3L_69_all30IS_observed_names.out

grep "FBti" 69_all30IS_observed_names.out > FBti_69_all30IS_observed_names.out

awk -F '\t' 'NR==FNR {id[$1]; next} $3 in id' /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/all30/TEs/2L_droseu14_all30_noinv_TEs_c FBti_69_all30IS_observed_names.out > 2L_FBti_69_all30IS_observed_names.out
awk -F '\t' 'NR==FNR {id[$1]; next} $3 in id' /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/all30/TEs/2R_droseu14_all30_noinv_TEs_c FBti_69_all30IS_observed_names.out > 2R_FBti_69_all30IS_observed_names.out
awk -F '\t' 'NR==FNR {id[$1]; next} $3 in id' /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/all30/TEs/3R_droseu14_all30_noinv_TEs_c FBti_69_all30IS_observed_names.out > 3R_FBti_69_all30IS_observed_names.out
awk -F '\t' 'NR==FNR {id[$1]; next} $3 in id' /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/all30/TEs/3L_droseu14_all30_noinv_TEs_c FBti_69_all30IS_observed_names.out > 3L_FBti_69_all30IS_observed_names.out


cat 2L_69_all30IS_observed_names.out 2L_FBti_69_all30IS_observed_names.out > 2L_FBti_69_all30IS_observed_names_unsorted.out
cat 2R_69_all30IS_observed_names.out 2R_FBti_69_all30IS_observed_names.out > 2R_FBti_69_all30IS_observed_names_unsorted.out
cat 3R_69_all30IS_observed_names.out 3R_FBti_69_all30IS_observed_names.out > 3R_FBti_69_all30IS_observed_names_unsorted.out
cat 3L_69_all30IS_observed_names.out 3L_FBti_69_all30IS_observed_names.out > 3L_FBti_69_all30IS_observed_names_unsorted.out


sort -n -k4,4 2L_FBti_69_all30IS_observed_names_unsorted.out > 2L_FBti_69_all30IS_observed_names_sorted.out
sort -n -k4,4 2R_FBti_69_all30IS_observed_names_unsorted.out > 2R_FBti_69_all30IS_observed_names_sorted.out
sort -n -k4,4 3R_FBti_69_all30IS_observed_names_unsorted.out > 3R_FBti_69_all30IS_observed_names_sorted.out
sort -n -k4,4 3L_FBti_69_all30IS_observed_names_unsorted.out > 3L_FBti_69_all30IS_observed_names_sorted.out


cat 2L_FBti_69_all30IS_observed_names_sorted.out 2R_FBti_69_all30IS_observed_names_sorted.out 3L_FBti_69_all30IS_observed_names_sorted.out 3R_FBti_69_all30IS_observed_names_sorted.out > 69_all30IS_observed_complete.out
