#!/usr/bin/env python3
# -*- coding: utf-8 -*- 

import sys
import pandas as pd
import numpy as np
import matplotlib 
matplotlib.use('agg') #avoid tkinter problem, that generate images without having a window appear 
import matplotlib.pyplot as plt
from pathlib import Path
import argparse



parser = argparse.ArgumentParser(description='Multiseries plot of number of DMIs average in time of multiple experiments')
parser.add_argument('--relative', '-r', 
    action='store_true',
    help='X axis relative instead of absolute, use max generation as 100',
    default=False)
parser.add_argument('--style', '-s', 
    help='Define colors using MatPlotLib styles',
    default="default")
parser.add_argument('--xmin', '-xm', type=int, default=0, help='Set xmin')
parser.add_argument('--xmax', '-xM', type=int, help='Set xmax')
parser.add_argument('--name', '-n', type=str, help='Set plot name', default="")
parser.add_argument('--stddev', '-sd', 
    action='store_true',
    help='Uses standard deviation instad of Standard error of the mean',
    default=False)


parser.add_argument('files', nargs='+')

args = parser.parse_args()





allfiles = '' #define empty string for all names in comparison

plt.style.use(args.style)

for file in args.files: # 

    print("Plotting " + file + "...")

    filename = Path(file).name.split("numHybsByNumofDMIs")[0] #extract the name of the experiment
    if (args.name):
        plotname = args.name
    else:
        plotname = filename
    allfiles += filename + '_' 
    df = pd.read_json(file) #read the file and put in the df

    if args.stddev:
        df_grouped = df.groupby(['generation'])['avgDMIsperHyb'].agg(['mean','std']) 
    else:
        df_grouped = df.groupby(['generation'])['avgDMIsperHyb'].agg(['mean','sem']) #group df using generation column but using only avgDMIsperHyb column. The aggregators are passed in tuple to add new columns.

    #add zero to the dataframe of the zero generation
    # new_row = pd.DataFrame({'mean':0, 'sem':0},index =[0])
    # df_grouped = pd.concat([new_row, df_grouped])


    
    if args.relative:
        divisor = df_grouped.index.max()
        relativeString = "Relative"
    else:
        divisor = 1
        relativeString = ""
    if args.stddev:
        plt.errorbar(df_grouped.index / divisor, df_grouped['mean'], yerr='std', data=df_grouped, label=plotname) #add a line with error standar for each file, notice that I am assuming here that the index is the generation, check df before using this and be sure that generation is the index 
    else:
        plt.errorbar(df_grouped.index / divisor, df_grouped['mean'], yerr='sem', data=df_grouped, label=plotname) #add a line with error standar for each file, notice that I am assuming here that the index is the generation, check df before using this and be sure that generation is the index 

# # #Save in png file and adjust some plot parameters




if args.relative:
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    plt.xlabel('Time')
else:
    plt.xlabel('Generations')


if args.xmax:
    finalxmax=args.xmax
    relativeString += "_[{}-{}]".format(args.xmin,args.xmax)
else:
    finalxmax=plt.axis()[1]

plt.xlim(args.xmin,finalxmax)
plt.ylim(0) #adjust y to start in x-intecerpt

plt.legend(loc='best')
plt.ylabel('# of DMIs')
plt.tight_layout() #fit everything in the canvas
plt.savefig(allfiles + "timeComparison" + relativeString + ".png")


if len(args.files) <= 1:
    # df_grouped.to_excel(allfiles + "timeComparison" + relativeString + ".xlsx", engine='xlsxwriter')
    df_grouped.to_csv(allfiles + "timeComparison" + relativeString + ".csv", index=True, index_label="generation_" + plotname, header=["mean_" + plotname, "+-"])
    # df_grouped.to_csv(allfiles + "timeComparison" + relativeString + ".csv", index=True, index_label="generation_" + plotname, header=["mean_" + plotname, "+-","nn"])
    #df.to_csv(allfiles + "pureDataFrame" + relativeString + ".csv", index=True, index_label="index")
    df.to_csv(allfiles + "pureDataFrame" + relativeString + ".csv", index=False, index_label="index", columns=["avgDMIsperHyb", "generation", "s"], header=["avgDMIsperHyb_" + plotname, "generations_" + plotname, "s_" + plotname])