#!/bin/sh

module load Python/2.7.15-foss-2018b

for f in $(seq 50 $END); do for i in $(seq 68 $END); do python /homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/analysis/parsing_baypass_SNPs_with_name.py ${i}_all30IS_${f}_observed_names.out 3 ${i}_all30IS_${f}_median_observed_names.out; done; done
