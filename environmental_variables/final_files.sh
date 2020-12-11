#!/bin/sh

module load Python/2.7.12-foss-2016b


python /homes/users/mbogaerts/scratch/Baypass/environmental_variables/scripts/all_variables_def.py \
temperature/all_populations_temperature.txt \
precipitation/all_populations_precipitation.txt \
evaporation/all_populations_evaporation.txt \
solar_radiation/all_populations_solar_radiation.txt \
soil/all_populations_soil.txt \
wind_u/all_populations_wind_u.txt \
wind_v/all_populations_wind_v.txt \
lighthours/Dros14_all_30pop_summary_daylight_def.txt \
/homes/users/mbogaerts/scratch/Baypass/environmental_variables/dros14_all_30pop/population_order.txt \
/homes/users/mbogaerts/scratch/Baypass/environmental_variables/dros14_all_30pop \
20 \
/homes/users/mbogaerts/scratch/Baypass/environmental_variables/dros14_all_30pop/droseu14_all30pop_environmental_variables.txt

