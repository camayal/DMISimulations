#!/usr/bin/env python3
# -*- coding: utf-8 -*- 

import sys
import json
from pandas import DataFrame
import numpy as np
import matplotlib 
matplotlib.use('agg') #avoid tkinter problem, that generate images without having a window appear 
import matplotlib.pyplot as plt
from pathlib import Path

"""
Count number of DI and MI in the json files 
"""
 

# [{"s":"1564200540644","generation":1001,"dmi":"Ba-CA-CB","abcaBc":0,"abcAbc":0,"aBcabc":0,"abCabc":0,"abcabC":0,"aBcaBc":0,"abCabC":0,"Abcabc":0,"AbcAbc":0,"AbcaBc":0,"ABcabc":0,"abcABc":0,"abCAbc":0,"AbcabC":0},

file = sys.argv[1] #get first arugment as file path
detailsFile = sys.argv[2] #get second argument need to be the resultTable.json file to extract the total number of hybrid individuals
plotname = Path(file).name.split("infoHybs")[0]

print("Plotting " + file + "...")

result = []
hybridpopsize = 1

#open resultTable.json to extract the total hybrid individuals 
with open(detailsFile) as details_json_file: #open file
    details = json.load(details_json_file)

with open(file) as json_file: #open file
    data = json.load(json_file) #parse json file
    for simulation, objects in enumerate(data,1): #for each object in json do:
        seed = objects["s"] #extract the seed to report
        generation = objects["generation"] #extract generation to report [important!]
        genotypes = list(objects.keys())[3:] #put in genotypes all values after 3 (excluded)
        individualDMIs = objects["dmi"].split("-") #get each dmi indidivually
        dmiMatrix = {}
        diMatrix = {}
        miMatrix = {}

        for index, genotype in enumerate(genotypes,3): #enter in a loop for each genotype in that object (aka generation)
            if int(list(objects.values())[index]) > 0: #check if the genotype has more than 0 individuals
                countdmi=[]
                countmi=[]
                countdi=[]

                for dmi in individualDMIs: #loop for each dmi
                    if dmi[0] in genotype and dmi[1] in genotype: #detect if each gene is in current genotype
                        countdmi.append(int(list(objects.values())[index])) #if both genes are  count as dmi

                        if dmi[0].islower() or dmi[1].islower(): #if one gene is lower, means muller type
                            countmi.append(int(list(objects.values())[index]))
                        else: #else means dob type
                            countdi.append(int(list(objects.values())[index]))

                #feed the dictionary
                if len(countdmi) > 0:
                    _previous =  dmiMatrix.get(len(countdmi), 0) #try to get the value x (number of dmis that have this genotype) in {x: y}, if is not return a zero
                    dmiMatrix[len(countdmi)] = _previous + sum(countdmi) / len(countdmi) #add the x to the dictionary and add the effective number of hybrid

                #see logic explaination in dmi part
                if len(countdi) > 0:
                    _previous =  diMatrix.get(len(countdi), 0)
                    diMatrix[len(countdi)] = _previous + sum(countdi) / len(countdi)

                if len(countmi) > 0:
                    _previous =  miMatrix.get(len(countmi), 0)
                    miMatrix[len(countmi)] = _previous + sum(countmi) / len(countmi)

        # print (seed)
        # print ("dmi: " + str(dmiMatrix))
        # print ("di: " + str(diMatrix))
        # print ("mi: " + str(miMatrix))
        avgdmis = 0
        avgdis = 0
        avgmis = 0 
        
        for sim in details:
            if sim['Seed'] == str(seed) and sim['Generation'] == generation:
                if sim['Pop4_Size']:
                    hybridpopsize = sim['Pop4_Size']

        # print(hybridpopsize)

        for dmis in dmiMatrix:
            avgdmis += ((dmis * dmiMatrix[dmis])/hybridpopsize)
        for dis in diMatrix:
            avgdis += ((dis * diMatrix[dis])/hybridpopsize)
        for mis in miMatrix:
            avgmis += ((mis * miMatrix[mis])/hybridpopsize)
        # print (avgdmis)
        # print (avgdis)
        # print (avgmis)
 
        # result.append([seed, generation, avgdmis, avgdis, avgmis]) #add as list to a final result list
        result.append([seed, generation, avgdis, avgmis]) #add as list to a final result list
        # print (result)



##Plot part

# #convert in result list dataframe
# df = DataFrame.from_records(result, columns=['Seed', 'Generation', 'avgDMIs', 'avgDIs', 'avgMIs']) 
df = DataFrame.from_records(result, columns=['Seed', 'Generation', 'D-type', 'M-type']) 

#Group data by Generation getting means
df_grouped = df.groupby(['Generation']).mean()

df_grouped_errors = df.groupby(['Generation']).agg(['mean','sem']) #group df using generation column but using only avgDMIsperHyb column. The aggregators are passed in tuple to add new columns.

# print(df_grouped)

#stack plot
plot = df_grouped.plot(kind='bar',stacked=True)

#Save in png file and adjust some plot parameters
plt.ylabel('Avg. incompat. per hybrid') #ylabel title
plt.xlabel('Generations')
plt.tight_layout() #fit everything in the canvas
# plt.xlim(xmin=0) #
fig = plot.get_figure()
fig.savefig(sys.argv[1].split(".")[0] + "_AVG_DnMis_intime.png")


plot = df_grouped.plot(kind='bar',stacked=False)
plt.ylabel('Avg. incompat. per hybrid') #ylabel title
plt.xlabel('Generations')
plt.tight_layout()
nonstacked = "_nonStacked_"
fig = plot.get_figure()
fig.savefig(sys.argv[1].split(".")[0] + nonstacked + "_AVG_DnMis_intime.png")


##Export to csv dfs
df_grouped.to_csv(sys.argv[1].split(".")[0] + "_AVG_DnMis_intime.csv", index=True, index_label="generation_DD-DA_" + plotname, header=["D-type_" + plotname, "M-type_" + plotname])
df.to_csv(sys.argv[1].split(".")[0] + "_AVG_RawDataFrame_DnMis_intime.csv", index=True)
df_grouped_errors.to_csv(sys.argv[1].split(".")[0] + "_AVG_DnMis_intime_withErrors.csv", index=True,index_label="generation_DD-DA_" + plotname, header=["D-type_" + plotname, "+-","M-type_" + plotname, "+-"])


