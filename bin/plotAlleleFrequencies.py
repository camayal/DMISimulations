#!/usr/bin/env python3
# -*- coding: utf-8 -*- 

import sys
import json
import pandas as pd
import matplotlib 
matplotlib.use('agg') #avoid tkinter problem, that generate images without having a window appear 
import matplotlib.pyplot as plt
import os.path
import argparse



parser = argparse.ArgumentParser(description='Multiseries plot of number of DMIs average in time of multiple experiments')
parser.add_argument('--nopop2',
    action='store_true',
    help='Plot population 2 (parent)',
    default=False)
parser.add_argument('--nopop3', 
    action='store_true',
    help='Plot population 3 (parent)',
    default=False)
parser.add_argument('--nopop4', 
    action='store_true',
    help='Plot population 4 (hybrids)',
    default=False)
parser.add_argument('--format', '-f', 
    help='Set the format for output, "svg", "png"',
    default="png")
parser.add_argument('--style', '-s', 
    help='Define colors using MatPlotLib styles',
    default="default")
parser.add_argument('--dmis', '-d', 
    help='Define list of DMIs manually to report',
    default=" ")

# parser.add_argument('--xmin', '-xm', type=int, default=0, help='Set xmin')
# parser.add_argument('--xmax', '-xM', type=int, help='Set xmax')
# parser.add_argument('--name', '-n', type=str, help='Set plot name', default="")

parser.add_argument('files', nargs='+')

args = parser.parse_args()


print("Plotting " + args.files[0] + "...")


seed = args.files[0].split("alleleFreq")[0]

df = pd.read_json(args.files[0])


fig, ax = plt.subplots(3, figsize=(4.5, 5))
plt.style.use(args.style)


#print(sum([args.nopop2, args.nopop3, args.nopop4]))


# Get DMI pairs from infoHybs.json file (suggested or passed in the second argument)
def getDMIpairs(infoHybsFile):
  print("\tGetting DMIs from " + infoHybsFile + "...")
  with open(infoHybsFile) as details_json_file: #open file
    details = json.load(details_json_file)
    for sim in details:
            if sim['s'] == str(seed):
                if sim['dmi']:
                    dmis = sim['dmi']
    return dmis


if len(args.files) < 2:
    # if os.path.isfile(seed + 'infoHybs.json'):
    #     dmiPairs =  getDMIpairs(seed + 'infoHybs.json')
    # else:
    if args.dmis:
        dmiPairs = args.dmis
else:
    dmiPairs = getDMIpairs(args.files[1])




data2 = [s for s in df if "_2" in s]
data3 = [s for s in df if "_3" in s]
data4 = [s for s in df if "_4" in s]

if data2 and not args.nopop2: 
    plotPop2 = df[data2].plot(x=df["G"], ax=ax[0], stacked=False, linewidth=0.5)
    plotPop2.set_ylabel('Allele frequency')
    plotPop2.set_xlabel('Parental pop. 1')
    plotPop2.legend(fontsize = 'xx-small',loc='upper left')
    plotPop2.set_ylim(-0.1,1.1)
    plotPop2.text(plotPop2.get_xlim()[0], plotPop2.get_ylim()[1] + 0.05, str(dmiPairs), size = "medium") 


if data4 and not args.nopop4: 
    plotPop4 = df[data4].plot(x=df["G"], ax=ax[1], stacked=False, marker='+', linewidth=0.5)
    plotPop4.set_ylabel('Allele frequency')
    plotPop4.set_xlabel('Hybrids')
    plotPop4.legend(fontsize = 'xx-small',loc='upper left')
    plotPop4.set_ylim(0,1)
    plotPop4.text(plotPop4.get_xlim()[0], plotPop4.get_ylim()[1] + 0.05, str(dmiPairs), size = "medium") 


if data3 and not args.nopop3: 
    plotPop3 = df[data3].plot(x=df["G"], ax=ax[2], stacked=False, linewidth=0.5)
    plotPop3.set_ylabel('Allele frequency')
    plotPop3.set_xlabel('Parental pop. 2')
    plotPop3.legend(fontsize = 'xx-small',loc='upper left')
    plotPop3.set_ylim(-0.1,1.1)
    plotPop3.text(plotPop3.get_xlim()[0], plotPop3.get_ylim()[1] + 0.05, str(dmiPairs), size = "medium")
    






plt.tight_layout()


plt.savefig(args.files[0].split(".")[0] + "_plot." + args.format)



#df.to_csv(sys.argv[1].split(".")[0] + ".csv", index=False)
