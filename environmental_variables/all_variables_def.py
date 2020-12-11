#!/usr/bin/python
# Created on 19th September 2018
# Parses the data downloaded and already modified from ERA5 database (Copernicus)
# Creates bio-variables for temperature similar to the ones in Worldclim
# USAGE: python /Users/pogo/Documents/Maria/Scripts/Baypass/environmental_variables/temperature_quarters.py temperature/all_populations_temperature.txt precipitation/all_populations_precipitation.txt temperature/population_order.txt /Users/pogo/Documents/Maria/Baypass/pruebas_cluster/dros14_fall 19 /Users/pogo/Documents/Maria/Baypass/pruebas_cluster/dros14_fall/final_file.txt

import sys
import os
import numpy
import operator

# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (* 100)
# BIO4 = Temperature Seasonality (standard deviation *100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter

# Crear archivin con los nombres de las poblaciones que aparecen en el nombre del archivo (ej "12_UK_Market+Harborough_2014" + ..._temperature_summary.txt)


all_pop_tmp = open(sys.argv[1],'r').readlines() # all_populations_temperature.txt
all_pop_prec = open(sys.argv[2],'r').readlines() # all_populations_precipitation.txt
all_pop_eva = open(sys.argv[3],'r').readlines() # all_populations_evaporation.txt
all_pop_rad = open(sys.argv[4],'r').readlines() # all_populations_solar_radiation.txt
all_pop_soil = open(sys.argv[5],'r').readlines() # all_populations_soil.txt 
#sys.argv[5]
all_pop_windu = open(sys.argv[6],'r').readlines() # all_populations_wind_u.txt
all_pop_windv = open(sys.argv[7],'r').readlines() # all_populations_wind_v.txt
all_pop_hours = open(sys.argv[8],'r').readlines() # Dros14_fall_summary_daylight_def.txt
#sys.argv[8]
names_list = open(sys.argv[9],'r').readlines() # Archivo con el orden de las variables
path = sys.argv[10] # Path including the final / where all the summary.txt are
num_pop = sys.argv[11] # Number of populations for setting up the matrix
final_file = open(sys.argv[12],'w')

w, h = int(num_pop)+1, 60;
var_matrix = [[0 for x in range(w)] for y in range(h)]

flag = 0
for names in names_list:
    flag = flag + 1
    n = names.split("\n")[0]
    #print(n)
    ## VARIABLE TEMPERATURE
    for popt in all_pop_tmp:
        name = popt.split("\t")[0]
        if n == name:
            sum_file = open(path + "/temperature/" + popt.split("\t")[0] + "_temperature_summary.txt",'r').readlines()
            avg_list_t = [] 
            min_list_t =[]
            max_list_t = []
            total_list = []
            bio2_t = []  
            for s in sum_file:
                if s.split("\t")[0] == "Month":
                    continue
                else:
                    avg_list_t.append(float(s.split("\t")[1]))
                    min_list_t.append(float(s.split("\t")[3]))
                    max_list_t.append(float(s.split("\t")[2]))
                    bio2_t.append(float(s.split("\t")[2]) - float(s.split("\t")[3])) #For BIO2
                    total_list.append(float(s.split("\t")[5]))
            max_tmp_t = max(max_list_t) # BIO5
            min_tmp_t = min(min_list_t) # BIO6
            ann_range_t = float(max_tmp_t) - float(min_tmp_t) #BIO7
            diurnal_range_t = sum(bio2_t) / float(len(bio2_t)) #BIO2
            isothermality_t = (diurnal_range_t/ann_range_t)*100 #BIO3
            # Putting all the bio-variables inside the matrix
            var_matrix[0][flag] = popt.split("\t")[0] #name
            var_matrix[1][flag] = popt.split("\t")[2] #BIO1
            var_matrix[1][0] = "BIO1"
            var_matrix[2][flag] = diurnal_range_t #BIO2
            var_matrix[2][0] = "BIO2"
            var_matrix[3][flag] = isothermality_t #BIO3
            var_matrix[3][0] = "BIO3"
            var_matrix[4][flag] = (float(popt.split("\t")[5])*100)/(float(popt.split("\t")[2])) #BIO4
            var_matrix[4][0] = "BIO4"
            var_matrix[5][flag] = max_tmp_t #BIO5
            var_matrix[5][0] = "BIO5"
            var_matrix[6][flag] = min_tmp_t #BIO6
            var_matrix[6][0] = "BIO6"
            var_matrix[7][flag] = ann_range_t #BIO7
            var_matrix[7][0] = "BIO7"
            # Calculating the warmest and coldest quarter:
            quarter_dict_t = {}
            for r in range(0,11):
                if r == 10:
                    quarter_t = float(avg_list_t[r]) + float(avg_list_t[r+1]) + float(avg_list_t[0])
                    quarter_dict_t["11-12-1"] = quarter_t
                elif r == 11:
                    quarter_t = float(avg_list_t[r]) + float(avg_list_t[0]) + float(avg_list_t[1])
                    quarter_dict_t["12-1-2"] = quarter_t
                else:
                    quarter_t = float(avg_list_t[r]) + float(avg_list_t[r+1]) + float(avg_list_t[r+2])
                    quarter_name_t = str(r) + "-" + str(r+1) + "-" + str(r+2)
                    quarter_dict_t[quarter_name_t] = quarter_t
            # Quarters
            max_quarter_t = max(quarter_dict_t.iteritems(), key=operator.itemgetter(1))[0]
            min_quarter_t = min(quarter_dict_t.iteritems(), key=operator.itemgetter(1))[0]
            var_matrix[10][flag] = (quarter_dict_t[max_quarter_t])/3 #BIO10
            var_matrix[10][0] = "BIO10"
            var_matrix[11][flag] = (quarter_dict_t[min_quarter_t])/3 #BIO11
            var_matrix[11][0] = "BIO11"
    ## VARIABLE PRECIPITATION
    for popp in all_pop_prec:
        name = popp.split("\t")[0]
        if n == name:
            sum_file_p = open(path + "/precipitation/" + popp.split("\t")[0] + "_precipitation_summary.txt",'r').readlines()
            total_prec = []  
            avg_list = []
            min_list = []
            max_list = []
            for sp in sum_file_p:
                if sp.split("\t")[0] == "Month":
                    continue
                else:
                    avg_list.append(float(sp.split("\t")[1]))
                    min_list.append(float(sp.split("\t")[3]))
                    max_list.append(float(sp.split("\t")[2]))
                    total_prec.append(float((sp.split("\t")[5]).split("\n")[0]))
            var_matrix[12][flag] = sum(total_prec) #BIO12
            var_matrix[12][0] = "BIO12"
            var_matrix[13][flag] = max(total_prec) # BIO13 (wettest month)
            var_matrix[13][0] = "BIO13"
            var_matrix[14][flag] = min(total_prec) # BIO14 (driest month)
            var_matrix[14][0] = "BIO14"
            var_matrix[15][flag] = (numpy.std(total_prec))/(1+((sum(total_prec))/12))*100 # BIO15
            var_matrix[15][0] = "BIO15"
            # Calculating the wettest and driest quarter:
            quarter_dict_p = {}
            for g in range(0,11):
                if g == 10:
                    quarter_p = float(total_prec[g]) + float(total_prec[g+1]) + float(total_prec[0])
                    quarter_dict_p["11-12-1"] = quarter_p
                elif g == 11:
                    quarter_p = float(total_prec[g]) + float(total_prec[0]) + float(total_prec[1])
                    quarter_dict_p["12-1-2"] = quarter_p
                else:
                    quarter_p = float(total_prec[g]) + float(total_prec[g+1]) + float(total_prec[g+2])
                    quarter_name_p = str(g) + "-" + str(g+1) + "-" + str(g+2)
                    quarter_dict_p[quarter_name_p] = quarter_p
            # Quarter
            max_quarter_p = max(quarter_dict_p.iteritems(), key=operator.itemgetter(1))[0]
            min_quarter_p = min(quarter_dict_p.iteritems(), key=operator.itemgetter(1))[0]
            var_matrix[16][flag] = (quarter_dict_p[max_quarter_p]) #BIO16
            var_matrix[16][0] = "BIO16"
            var_matrix[17][flag] = (quarter_dict_p[min_quarter_p]) #BIO17
            var_matrix[17][0] = "BIO17"
            var_matrix[18][flag] = (quarter_dict_p[max_quarter_t]) #BIO18
            var_matrix[18][0] = "BIO18"
            var_matrix[19][flag] = (quarter_dict_p[min_quarter_t]) #BIO19
            var_matrix[19][0] = "BIO19"
            var_matrix[8][flag]= (quarter_dict_t[max_quarter_p])/3 #BIO8
            var_matrix[8][0] = "BIO8"
            var_matrix[9][flag]= (quarter_dict_t[min_quarter_p])/3 #BIO9
            var_matrix[9][0] = "BIO9"
    ## VARIABLE EVAPORATION
    for pope in all_pop_eva:
        name = pope.split("\t")[0]
        if n == name:
            sum_file_e = open(path + "/evaporation/" + pope.split("\t")[0] + "_evaporation_summary.txt",'r').readlines()
            avg_list_e = [] 
            min_list_e =[]
            max_list_e = []
            total_eva = []
            bio2_e = []  
            for e in sum_file_e:
                if e.split("\t")[0] == "Month":
                    continue
                else:
                    avg_list_e.append(float(e.split("\t")[1]))
                    min_list_e.append(float(e.split("\t")[3]))
                    max_list_e.append(float(e.split("\t")[2]))
                    bio2_e.append(float(e.split("\t")[2]) - float(e.split("\t")[3])) #For BIO2
                    total_eva.append(e.split("\t")[5])
            max_tmp_e = min(max_list_e) #BIO5
            min_tmp_e = max(min_list_e) #BIO6
            ann_range_e = float(max_tmp_e) - float(min_tmp_e) #BIO7
            diurnal_range_e = sum(bio2_e) / float(len(bio2_e)) #BIO2
            isothermality_e = (diurnal_range_e/ann_range_e)*100 #BIO3
            # Putting all the bio-variables inside the matrix
            var_matrix[20][flag] = pope.split("\t")[2] #BIO1
            var_matrix[20][0] = "BIO1E"
            var_matrix[21][flag] = diurnal_range_e #BIO2
            var_matrix[21][0] = "BIO2E"
            var_matrix[22][flag] = isothermality_e #BIO3
            var_matrix[22][0] = "BIO3E"
            var_matrix[23][flag] = float(pope.split("\t")[5]) #BIO4
            var_matrix[23][0] = "BIO4E"
            var_matrix[24][flag] = min(max_list_e) #BIO5
            var_matrix[24][0] = "BIO5E"
            var_matrix[25][flag] =max(min_list_e) #BIO6
            var_matrix[25][0] = "BIO6E"
            var_matrix[26][flag] = ann_range_e #BIO7
            var_matrix[26][0] = "BIO7E"
             # Calculating the warmest and coldest quarter:
            quarter_dict_e = {}
            for e in range(0,11):
                if e == 10:
                    quarter_e = float(avg_list_e[e]) + float(avg_list_e[e+1]) + float(avg_list_e[0])
                    quarter_dict_e["11-12-1"] = quarter_e
                elif e == 11:
                    quarter_e = float(avg_list_e[e]) + float(avg_list_e[0]) + float(avg_list_e[1])
                    quarter_dict_t["12-1-2"] = quarter_t
                else:
                    quarter_e = float(avg_list_e[e]) + float(avg_list_e[e+1]) + float(avg_list_e[e+2])
                    quarter_name_e = str(e) + "-" + str(e+1) + "-" + str(e+2)
                    quarter_dict_e[quarter_name_e] = quarter_e
            var_matrix[27][flag]= (quarter_dict_e[max_quarter_p])/3 #BIO8
            print("quarter problematico")
            print(max_quarter_p)
            print(quarter_dict_e[max_quarter_p])
            var_matrix[27][0] = "BIO8E"
            var_matrix[28][flag]= (quarter_dict_e[min_quarter_p])/3 #BIO9
            var_matrix[28][0] = "BIO9E"        
            var_matrix[29][flag] = (quarter_dict_e[max_quarter_t])/3 #BIO10
            var_matrix[29][0] = "BIO10E"
            var_matrix[30][flag] = (quarter_dict_e[min_quarter_t])/3 #BIO11
            var_matrix[30][0] = "BIO11E"
    ## VARIABLE SOLAR RADIATION
    for popr in all_pop_rad:
        name = popr.split("\t")[0]
        if n == name:
            sum_file_r = open(path + "/solar_radiation/" + popr.split("\t")[0] + "_solar_radiation_summary.txt",'r').readlines()
            avg_list_r = [] 
            min_list_r =[]
            max_list_r = []
            total_rad = []
            bio2_r = []  
            for r in sum_file_r:
                if r.split("\t")[0] == "Month":
                    continue
                else:
                    avg_list_r.append(float(r.split("\t")[1]))
                    min_list_r.append(float(r.split("\t")[3]))
                    max_list_r.append(float(r.split("\t")[2]))
                    bio2_r.append(float(r.split("\t")[2]) - float(r.split("\t")[3])) #For BIO2
                    total_rad.append(r.split("\t")[5])
            max_tmp_r = max(max_list_r)
            min_tmp_r = min(min_list_r)
            ann_range_r = float(max_tmp_r) - float(min_tmp_r) #BIO7
            diurnal_range_r = sum(bio2_r) / float(len(bio2_r)) #BIO2
            isothermality_r = (diurnal_range_r/ann_range_r)*100 #BIO3
            # Putting all the bio-variables inside the matrix
            var_matrix[31][flag] = popr.split("\t")[2] #BIO1
            var_matrix[31][0] = "BIO1S"
            var_matrix[32][flag] = diurnal_range_r #BIO2
            var_matrix[32][0] = "BIO2S"
            var_matrix[33][flag] = isothermality_r #BIO3
            var_matrix[33][0] = "BIO3S"
            var_matrix[34][flag] = (float(popr.split("\t")[5])*100)/(float(popr.split("\t")[2])) #BIO4
            var_matrix[34][0] = "BIO4S"
            var_matrix[35][flag] = max_tmp_r #BIO5
            var_matrix[35][0] = "BIO5S"
            var_matrix[36][flag] = min_tmp_r #BIO6
            var_matrix[36][0] = "BIO6S"
            var_matrix[37][flag] = ann_range_r #BIO7
            var_matrix[37][0] = "BIO7S"
    ## VARIABLE WIND U
    for popw in all_pop_windu:
        name = popw.split("\t")[0]
        if n == name:
            sum_file_w = open(path + "/wind_u/" + popw.split("\t")[0] + "_wind_u_summary.txt",'r').readlines()
            avg_list_w = [] 
            min_list_w =[]
            max_list_w = []
            total_windu = []
            bio2_w = []  
            for w in sum_file_w:
                if w.split("\t")[0] == "Month":
                    continue
                else:
                    avg_list_w.append(abs(float(w.split("\t")[1])))
                    min_list_w.append(abs(float(w.split("\t")[3])))
                    max_list_w.append(abs(float(w.split("\t")[2])))
                    bio2_w.append(abs(float(w.split("\t")[2])) - abs(float(w.split("\t")[3]))) #For BIO2
                    total_windu.append(abs(float(s.split("\t")[5])))
            max_tmp_w = max(max_list_w)
            min_tmp_w = min(min_list_w)
            ann_range_w = float(max_tmp_w) - float(min_tmp_w) #BIO7
            diurnal_range_w = sum(bio2_w) / float(len(bio2_w)) #BIO2
            isothermality_w = (diurnal_range_w/ann_range_w)*100 #BIO3
            # Putting all the bio-variables inside the matrix
            var_matrix[38][flag] = popw.split("\t")[2] #BIO1
            var_matrix[38][0] = "BIO1WU"
            var_matrix[39][flag] = diurnal_range_w #BIO2
            var_matrix[39][0] = "BIO2WU"
            var_matrix[40][flag] = isothermality_w #BIO3
            var_matrix[40][0] = "BIO3WU"
            var_matrix[41][flag] = (abs(float(popw.split("\t")[5])*100))/abs(float(popw.split("\t")[2])) #BIO4
            var_matrix[41][0] = "BIO4WU"
            var_matrix[42][flag] = max_tmp_w #BIO5
            var_matrix[42][0] = "BIO5WU"
            var_matrix[43][flag] = min_tmp_w #BIO6
            var_matrix[43][0] = "BIO6WU"
            var_matrix[44][flag] = ann_range_w #BIO7
            var_matrix[44][0] = "BIO7WU"   
    ## VARIABLE WIND V
    for popv in all_pop_windv:
        name = popv.split("\t")[0]
        if n == name:
            sum_file_v = open(path + "/wind_v/" + popv.split("\t")[0] + "_wind_v_summary.txt",'r').readlines()
            avg_list_v = [] 
            min_list_v =[]
            max_list_v = []
            total_windv = []
            bio2_v = []  
            for v in sum_file_v:
                if v.split("\t")[0] == "Month":
                    continue
                else:
                    avg_list_v.append(abs(float(v.split("\t")[1])))
                    min_list_v.append(abs(float(v.split("\t")[3])))
                    max_list_v.append(abs(float(v.split("\t")[2])))
                    bio2_v.append(abs(float(v.split("\t")[2])) - abs(float(v.split("\t")[3]))) #For BIO2
                    total_windv.append(abs(float(v.split("\t")[5])))
            max_tmp_v = max(max_list_v)
            min_tmp_v = min(min_list_v)
            ann_range_v = float(max_tmp_v) - float(min_tmp_v) #BIO7
            diurnal_range_v = sum(bio2_v) / float(len(bio2_v)) #BIO2
            isothermality_v = (diurnal_range_v/ann_range_v)*100 #BIO3
            # Putting all the bio-variables inside the matrix
            var_matrix[45][flag] = popv.split("\t")[2] #BIO1
            var_matrix[45][0] = "BIO1WV"
            var_matrix[46][flag] = diurnal_range_v #BIO2
            var_matrix[46][0] = "BIO2WV"
            var_matrix[47][flag] = isothermality_v #BIO3
            var_matrix[47][0] = "BIO3WV"
            var_matrix[48][flag] = (abs(float(popv.split("\t")[5])*100))/abs(float(popv.split("\t")[2])) #BIO4
            var_matrix[48][0] = "BIO4WV"
            var_matrix[49][flag] = max_tmp_v #BIO5
            var_matrix[49][0] = "BIO5WV"
            var_matrix[50][flag] = min_tmp_v #BIO6
            var_matrix[50][0] = "BIO6WV"
            var_matrix[51][flag] = ann_range_v #BIO7
            var_matrix[51][0] = "BIO7WV"
    ## VARIABLE SOIL
    for popl in all_pop_soil:
        name = popl.split("\t")[0]
        if n == name:
            if popl.split("\t")[0] == "Population":
                continue
            else:
                if float(popl.split("\t")[5]) == 0:
                    var_matrix[52][flag] = popl.split("\t")[2]
                else:
                    var_matrix[52][flag] = "X"
    var_matrix[52][0] = "SOILTYPE"
    ## DAYLIGHT HOURS
    for poph in all_pop_hours:
        name = poph.split("\t")[0]
        if n == name:
            sum_file_h = open(path + "/lighthours/" + poph.split("\t")[0] + "_lighthours_summary.txt",'r').readlines()
            avg_list_h = [] 
            min_list_h =[]
            max_list_h = []
            total_hours = []
            bio2_h = []  
            for h in sum_file_h:
                avg_list_h.append(float(h.split("\t")[1]))
                min_list_h.append(float(h.split("\t")[3]))
                max_list_h.append(float(h.split("\t")[2]))
                bio2_h.append(float(h.split("\t")[2]) - float(h.split("\t")[3])) #For BIO2
            max_tmp_h = max(max_list_h)
            min_tmp_h = min(min_list_h)
            ann_range_h = float(max_tmp_h) - float(min_tmp_h) #BIO7
            diurnal_range_h = sum(bio2_h) / float(len(bio2_h)) #BIO2
            isothermality_h = (diurnal_range_h/ann_range_h)*100 #BIO3
             # Putting all the bio-variables inside the matrix
            var_matrix[53][flag] = poph.split("\t")[2] #BIO1
            var_matrix[53][0] = "BIO1H"
            var_matrix[54][flag] = diurnal_range_h #BIO2
            var_matrix[54][0] = "BIO2H"
            var_matrix[55][flag] = isothermality_h #BIO3
            var_matrix[55][0] = "BIO3H"
            var_matrix[56][flag] = (float(poph.split("\t")[5])*100)/(float(poph.split("\t")[2])) #BIO4
            var_matrix[56][0] = "BIO4H"
            var_matrix[57][flag] = max_tmp_h #BIO5
            var_matrix[57][0] = "BIO5H"
            var_matrix[58][flag] = min_tmp_h #BIO6
            var_matrix[58][0] = "BIO6H"
            var_matrix[59][flag] = ann_range_h #BIO7
            var_matrix[59][0] = "BIO7H"
print(var_matrix)

for v in var_matrix:
    print >> final_file, v
