#!/usr/bin/python

# Created 28th August 2018 (finished 6th September)
# It parses the data variables obtained from the ERA5 dataset - Copernicus
# Usage: python parse_grib_neutral.py coordinates_file grib_file /Users/pogo/Documents/Maria/Baypass/environmental_variables/precipitation/ temperature

import os
import sys
import numpy
import time
from datetime import date
from datetime import timedelta
from dateutil import rrule
from calendar import monthrange
from datetime import datetime
import numpy as np

# https://docs.python.org/2/library/datetime.html#date-objects
# First we will need to have the information from each place
# Extract with the coordinates from each year

coordinates = open(sys.argv[1],'r').readlines()
grib_file = sys.argv[2]
path = sys.argv[3]
variable = sys.argv[4]

#for c in coordinates:
#	number = c.split("\t")[0]
#	code = c.split("\t")[1]
#	country = c.split("\t")[2]
#	place = c.split("\t")[3]
#	year = str(c.split("\t")[4])[:4] # Collection year
#	latitude = c.split("\t")[5]
#	longitude_not_parsed = c.split("\t")[6]
#	longitude = longitude_not_parsed.split("\n")[0]
#	print("grib_ls -l " + str(latitude) + "," + str(longitude) + " " + grib_file + " > " + path + number + "_" + code + "_" +  place + "_" + variable + "_" + year + ".txt")
#	os.system("grib_ls -l " + str(latitude) + "," + str(longitude) + " " + grib_file + " > " + path + number + "_" + code + "_" +  place + "_" + variable + "_" + year + ".txt")

general_file = open(path + "all_populations_" + variable + ".txt" ,'w')
print >> general_file, "Population\tCollection date\tAverage " + variable + " \tMaximum\tMinimum\tStandard deviation"
for c in coordinates:
    number = c.split("\t")[0]
    code = c.split("\t")[1]
    country = c.split("\t")[2]
    place = c.split("\t")[3]
    collection = c.split("\t")[4] # Collection date
    # if str(collection[6:8]) == "31":
    #     collection = str(collection[:6]) + "30"
    print("Definitive collection date: " + collection)
    year = str(collection)[:4] # Collection year
    grib_file = open(path + number + "_" + code + "_" +  place + "_" + variable + "_" + year + ".txt", 'r').readlines()
    grib_line = '\t'.join(grib_file[-6].split())
    grib_point = grib_line.split("\t")[3]
    grib_value = grib_point.split("#")[1]
    coll_date = date(int(str(collection)[:4]), int(str(collection)[4:6]), int(str(collection)[6:8]))
    deltatime = coll_date - date(int(str(collection)[:4])-1, int(str(collection)[4:6]), int(str(collection)[6:8]))
    start_date = date(int(str(collection)[:4]), int(str(collection)[4:6]), int(str(collection)[6:8]))
    end_date = date(int(str(collection)[:4])-1, int(str(collection)[4:6]), int(str(collection)[6:8]))
    print("Definitive end_date: " + str(end_date))
    max_temp = {}
    min_temp = {}
    months = ""
    days = ""
    years = {}
    max_year = {}
    min_year = {}
    total_year = {}
    final_file = open(path + number + "_" + code + "_" + place + "_" + year + "_" + variable + "_summary.txt", 'w')
    print >> final_file, "Month\tAverage\tMax\tMin\tStandard deviation\tTotal per month"
    while end_date < start_date:
        os.system("grep '" + str(end_date.strftime("%Y%m%d")) + "' " + path + number + "_" + code + "_" +  place + "_" + variable + "_" + year + ".txt > " + path + place + "_" + str(number) + "_" + variable + "_" + str(end_date.strftime("%Y%m%d") + ".txt"))
        day_file = open(path + place + "_" + str(number) + "_" + variable + "_" + str(end_date.strftime("%Y%m%d")) + ".txt",'r').readlines()
        print(str(end_date.strftime("%Y%m%d")))
        day = []
        # Average variable
        value_year = years.get(str(end_date.strftime("%Y")), "empty")
        if value_year == "empty":
            years[str(end_date.strftime("%Y"))] = {}
            value_month = years[str(end_date.strftime("%Y"))].get(str(end_date.strftime("%m")), "empty")
            if value_month == "empty":
                years[str(end_date.strftime("%Y"))][str(end_date.strftime("%m"))] = {}
        else:
            value_month = years[str(end_date.strftime("%Y"))].get(str(end_date.strftime("%m")), "empty")
            if value_month == "empty":
                years[str(end_date.strftime("%Y"))][str(end_date.strftime("%m"))] = {}
        # Maximum variable
        value_max_year = max_year.get(str(end_date.strftime("%Y")), "empty")
        if value_max_year == "empty":
            max_year[str(end_date.strftime("%Y"))] = {}
            value_max_month = max_year[str(end_date.strftime("%Y"))].get(str(end_date.strftime("%m")), "empty")
            if value_max_month == "empty":
                max_year[str(end_date.strftime("%Y"))][str(end_date.strftime("%m"))] = {}
        else:
            value_max_month = max_year[str(end_date.strftime("%Y"))].get(str(end_date.strftime("%m")), "empty")
            if value_max_month == "empty":
                max_year[str(end_date.strftime("%Y"))][str(end_date.strftime("%m"))] = {}
        # Minimum variable
        value_min_year = min_year.get(str(end_date.strftime("%Y")), "empty")
        if value_min_year == "empty":
            min_year[str(end_date.strftime("%Y"))] = {}
            value_min_month = min_year[str(end_date.strftime("%Y"))].get(str(end_date.strftime("%m")), "empty")
            if value_min_month == "empty":
                min_year[str(end_date.strftime("%Y"))][str(end_date.strftime("%m"))] = {}
        else:
            value_min_month = min_year[str(end_date.strftime("%Y"))].get(str(end_date.strftime("%m")), "empty")
            if value_min_month == "empty":
                min_year[str(end_date.strftime("%Y"))][str(end_date.strftime("%m"))] = {}
        # Total variable
        total_value_year = total_year.get(str(end_date.strftime("%Y")), "empty")
        if total_value_year == "empty":
            total_year[str(end_date.strftime("%Y"))] = {}
            total_value_month = total_year[str(end_date.strftime("%Y"))].get(str(end_date.strftime("%m")), "empty")
            if total_value_month == "empty":
                total_year[str(end_date.strftime("%Y"))][str(end_date.strftime("%m"))] = {}
        else:
            total_value_month = total_year[str(end_date.strftime("%Y"))].get(str(end_date.strftime("%m")), "empty")
            if total_value_month == "empty":
                total_year[str(end_date.strftime("%Y"))][str(end_date.strftime("%m"))] = {}
        for d in day_file:
            line = '\t'.join(d.split())
            if grib_value == "1":
                day.append(float(line.split("\t")[10]))
            elif grib_value == "2":
                day.append(float(line.split("\t")[11]))
            elif grib_value == "3":
                day.append(float(line.split("\t")[12]))
            elif grib_value == "4":
                day.append(float(line.split("\t")[13]))
        years[str(end_date.strftime("%Y"))][str(end_date.strftime("%m"))][str(end_date.strftime("%d"))] = (sum(float(g) for g in day))/(len(day))
        max_year[str(end_date.strftime("%Y"))][str(end_date.strftime("%m"))][str(end_date.strftime("%d"))] = float(max(day))
        # https://stackoverflow.com/questions/27966757/find-min-value-in-array-0
        minval = min(i for i in day if i > 0)
        min_year[str(end_date.strftime("%Y"))][str(end_date.strftime("%m"))][str(end_date.strftime("%d"))] = minval
        total_year[str(end_date.strftime("%Y"))][str(end_date.strftime("%m"))][str(end_date.strftime("%d"))] = sum(float(g) for g in day)
        print(day)
        print(max(day))
        print(min(day))
        max_temp[str(end_date.strftime("%Y%m%d"))] = max(day)
        min_temp[str(end_date.strftime("%Y%m%d"))] = min(day)
        os.system("rm " + path + place + "_" + str(number) + "_" + variable + "_" + str(end_date.strftime("%Y%m%d") + ".txt"))
        end_date = end_date + timedelta(days=1)
    month_list = []
    collection_year = int(str(collection)[:4])
    collection_month = int(str(collection)[4:6])
    collection_day = int(str(collection)[6:8])
    for month_number in range(0,13):
        new_year = collection_year
        new_month = collection_month - month_number
        new_day = collection_day
        if new_month < 1:
            new_year = collection_year - 1
            new_month = 12 - (month_number - collection_month)
        max_day_number = monthrange(new_year, new_month)[1]
        if collection_day > max_day_number:
            new_day = max_day_number
        month_list.append(datetime(new_year, new_month, new_day, 0, 0))
    month_list = list(reversed(month_list))
    all_months = []
    max_month = []
    min_month = []
    for m in range(0,12):
        day_list = []
        max_list = []
        min_list = []
        total_list = []
        for dt in rrule.rrule(rrule.DAILY, dtstart = month_list[m], until=month_list[m+1]):
            if month_list[m] <= dt < month_list[m+1]:
                print(dt.strftime("%Y%m%d"))
                day_list.append(years[str(dt.strftime("%Y"))][str(dt.strftime("%m"))][str(dt.strftime("%d"))])
                max_list.append(max_year[str(dt.strftime("%Y"))][str(dt.strftime("%m"))][str(dt.strftime("%d"))])
                print(max_list)
                min_list.append(min_year[str(dt.strftime("%Y"))][str(dt.strftime("%m"))][str(dt.strftime("%d"))])
                total_list.append(total_year[str(dt.strftime("%Y"))][str(dt.strftime("%m"))][str(dt.strftime("%d"))])
            else:
                if len(day_list) == 0:
                    continue
                else:
                    max_month.append(numpy.mean(max_list))
                    print(dt.strftime("%Y%m"))
                    print("Max list")
                    print(max_month)
                    min_month.append(numpy.mean(min_list))
                    all_months.append(numpy.mean(day_list))
                    print >> final_file, dt.strftime("%Y%m") + "\t" + str(numpy.mean(day_list)) + "\t" + str(numpy.mean(max_list)) + "\t" + str(numpy.mean(min_list)) + "\t" + str(numpy.std(day_list)) + "\t" + str(sum(total_list))
                    print(total_list)
    print >> general_file, number + "_" + code + "_" + place + "_" + year + "\t" + collection + "\t" + str(numpy.mean(all_months)) + "\t" + str(max(max_month)) + "\t" + str(min(min_month)) + "\t" + str(numpy.std(all_months))
