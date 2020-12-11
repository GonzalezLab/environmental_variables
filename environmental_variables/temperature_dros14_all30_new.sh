#!/bin/sh


module load eccodes/2.8.2-foss-2016b

module load Python/2.7.12-foss-2016b

python /homes/users/mbogaerts/scratch/Baypass/environmental_variables/scripts/grib_to_table_time_complete_fast.py /homes/users/mbogaerts/scratch/Baypass/environmental_variables/dros14_all_30pop/Dros14_all_30pop.txt /homes/users/mbogaerts/scratch/Baypass/environmental_variables/2013_2014_2015_2016_temperature.grib /homes/users/mbogaerts/scratch/Baypass/environmental_variables/dros14_all_30pop/temperature/ temperature

