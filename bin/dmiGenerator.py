#!/usr/bin/env python3
# -*- coding: utf-8 -*- 

"""
Generate from a given secuence of genes (or a random sequence) all potential DMIs
"""

import argparse
import random

#get some arguments
parser = argparse.ArgumentParser(description="Script find the potential DMIs given a sequence of loci")
parser.add_argument("-v", "--verbose", dest="verbose", help="Show all information", action="store_true")
parser.add_argument("-g", "--graphical", dest="graphical", help="Show a graphical output", action="store_true")
parser.add_argument("-n", "--number",  metavar='N', type=int, default=3, dest="number", help="Generate a random sequence of N genes")
parser.add_argument("input", metavar="AbCDe", nargs='?', type=str, help="Genotype of lineage 1")
parser.add_argument("-s", "--string", metavar="\",\"", dest="string", help="Return as a string using a separator")
parser.add_argument("-t", "--stats", dest="stats", help="Output only stats: genes, dmis, ddis, adis, orrEstimate", action="store_true")

#put them into args
args = parser.parse_args()


##Experimental
#alphabets (pay attention A != Α != А, in latin, greek and russian)
#Latin (26): ABCDEFGHIJKLMNOPQRSTUVWXYZ // abcdefghijklmnopqrstuvwxyz
#Greek (24): ΑΒΓΔΕΖΗΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩϴ // αβγδεζηθικλμνξοπρστυφχψωθ
#Russian (31): АБВГДЕЖЗИКЛМНОПРСТУФХЦЧШЩЪЫЬЭЮЯ // абвгдежзиклмнопрстуфхцчшщъыьэюя
##Danger zone (some versions of Unicode cn lack of those)
#Georgian (38): ԱԲԳԴԵԶԷԸԹԺԻԼԽԾԿՀՁՂՃՄՅՆՇՈՉՊՋՌՍՎՏՐՑՒՓՔՕՖ // աբգդեզէըթժիլխծկհձղճմյնշոչպջռսվտրցւփքօֆ
#Georgian (40): ႠႡႢႣႤႥႦჁႧႨႩႪႫႬჂႭႮႯႰႱႲჃႳჇႴႵႶႷႸႹႺႻႼႽႾჄႿჀჅჍ // ⴀⴁⴂⴃⴄⴅⴆⴡⴇⴈⴉⴊⴋⴌⴢⴍⴎⴏⴐⴑⴒⴣⴓⴧⴔⴕⴖⴗⴘⴙⴚⴛⴜⴝⴞⴤⴟⴠⴥⴭ
#Coptic (33): ⲀⲂⲄⲆⲈⲊⲌⲎⲐⲒⲔⲖⲘⲚⲜⲞⲠϤⲢⲤⲦⲨⲪⲬⲮⲰⳀϢϦϨϪϬϮ // ⲁⲃⲅⲇⲉⲋⲍⲏⲑⲓⲕⲗⲙⲛⲝⲟⲡϥⲣⲥⲧⲩⲫⲭⲯⲱⳁϣϧϩϫϭϯ


#in case use -n generate a random sequence
if args.input is None:
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZΑΒΓΔΕΖΗΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩϴАБВГДЕЖЗИКЛМНОПРСТУФХЦЧШЩЪЫЬЭЮЯԱԲԳԴԵԶԷԸԹԺԻԼԽԾԿՀՁՂՃՄՅՆՇՈՉՊՋՌՍՎՏՐՑՒՓՔՕՖႠႡႢႣႤႥႦჁႧႨႩႪႫႬჂႭႮႯႰႱႲჃႳჇႴႵႶႷႸႹႺႻႼႽႾჄႿჀჅჍⲀⲂⲄⲆⲈⲊⲌⲎⲐⲒⲔⲖⲘⲚⲜⲞⲠϤⲢⲤⲦⲨⲪⲬⲮⲰⳀϢϦϨϪϬϮ"
    lineage1 = ""
    for i in range(args.number):
        if random.choice([True, False]):
            lineage1 += alphabet[i].lower()
        else:
            lineage1 += alphabet[i].upper()
        lineage1 = ''.join(random.sample(lineage1,len(lineage1)))
else:
    lineage1 = args.input



lineage2 = lineage1.swapcase()  #create the lineage 2 secuence (complement of lineage 1)

dmis = []
dmiInfo = {}
ancestraInDMIInfo = []
miCount = 0
diCount = 0

#for both lineages reconstruct the history, slicing from 0 to current index and from index + 1 to end, this last
#part in lowercase to create the non-mutate time
#ABc > 1st slice: A, 2nd: AB, 3rd: ABc.  In the other slicing process, 1stB: Bc (lowercase)> bc, 2ndB: c, 3rdB: ""
#After that join with "" 1st, 1stB, 2st, 2stB, ...  Same for both lineages
lineage1_history = ["".join([lineage1[:index + 1], lineage1[index + 1:].lower()]) for index, gene in enumerate(lineage1)]
lineage2_history = ["".join([lineage2[:index + 1], lineage2[index + 1:].lower()]) for index, gene in enumerate(lineage2)]


#do a loop in the first lineage history, saving the index.
#do a loop in the first genotype, in the first letter, check if have or not a mutation
    #if has the mutation iterate with all genes in lineage 2 in that time (or index)
    #if does not have a mutation move to the lineage 2 and does the iteration with all genes in lineage 1
#while is getting all possible combination, do the evauation, if Aa noDMI, if pair is in the genotype is not a DMI.
#But when the pair is not in the genotype is a potential DMI
for index,genotype_l1 in enumerate(lineage1_history):
    if genotype_l1[index].isupper():
        for gene_l1 in genotype_l1[index]:
            for comparison in lineage2_history[index]:
                pair = gene_l1 + comparison
                #check if the pair is in the genotype of lineage 1 (using a loop to avoid non continuous false negative, Ac is in ABc)
                #also avoid same gene comparing both in lowercase
                if sum([c in pair for c in genotype_l1]) < len(pair) and gene_l1.lower() != comparison.lower():
                    dmis.append(pair) #add dmi to list
                    #store some information about this dmi
                    dmiInfo.setdefault(gene_l1, {})
                    dmiInfo[gene_l1].setdefault("index", index)
                    dmiInfo[gene_l1].setdefault("comparison", []).append(comparison)
                    ancestraInDMIInfo.append(comparison + str(index))


    else:
        for gene_l2 in lineage2_history[index][index]:
            if gene_l2.isupper():
                for comparison in genotype_l1:
                    pair = gene_l2 + comparison
                    #check if the pair is in the genotype of lineage 2 (using a loop to avoid non continuous false negative, Ac is in ABc)
                    #also avoid same gene comparing both in lowercase
                    if sum([c in pair for c in lineage2_history[index]]) < len(pair) and gene_l2.lower() != comparison.lower():
                        dmis.append(pair)#add dmi to list
                        #store some information about this dmi
                        dmiInfo.setdefault(gene_l2, {})
                        dmiInfo[gene_l2].setdefault("index", index)
                        dmiInfo[gene_l2].setdefault("comparison", []).append(comparison)
                        ancestraInDMIInfo.append(comparison + str(index))

                    
#print some extra information if verbose is activated        
if args.verbose:
    K = len(lineage1) 
    print ("Number of genes: " + str(K))
    print ("Original input (final genotype for lineage 1): " + lineage1)
    print ("Final genotype for lineage 2): " + lineage2)
    print ("Mutation per time in lineage 1: " + str(lineage1_history))
    print ("Mutation per time in lineage 2: " + str(lineage2_history))
    for index in range(len(dmis)):
        if str(dmis[index][0]).islower() or str(dmis[index][1]).islower():
            miCount += 1
        else:
            diCount += 1
    print ("Number of DMIs pairs: " + str(len(dmis)) + " (DI: "+ str(diCount) + ", MI: " + str(miCount) + ")" )
    print ("Orr`s DMIs expection (K(K-1)/2): " + str(int(K * (K - 1) / 2)))
    print ("\nPair of genes with potential DMI: " if len(dmis) > 0 else "No DMIs found!")


#print final result
if not args.stats:
    if len(dmis) > 0:
        if args.string:
            print (args.string.join(dmis))
        else:
            print (dmis)






# print dmiInfo
# print ancestraInDMIInfo


#graphical view
if args.graphical:
    # print ("\nGraphical output")
    for i in range(len(lineage1_history)):
        for j in range(len(lineage1_history)):
            # print "|" + lineage1_history[j][i] + "    " + lineage2_history[j][i] + "| ",
           
            #print in colors the table
            if lineage1_history[j][i] in dmiInfo and dmiInfo[lineage1_history[j][i]]["index"] == j or lineage1_history[j][i] + str(j) in ancestraInDMIInfo:
                l1ini = '\x1b[91m'
                l1end = '\x1b[0m'
            else:
                l1ini = ''
                l1end = '' 
            if lineage2_history[j][i] in dmiInfo and dmiInfo[lineage2_history[j][i]]["index"] == j or lineage2_history[j][i] + str(j) in ancestraInDMIInfo:
                l2ini = '\x1b[91m'
                l2end = '\x1b[0m'
            else:
                l2ini = ''
                l2end = ''
            print ("|" + l1ini + lineage1_history[j][i] + l1end + "    " +l2ini+ lineage2_history[j][i] + l2end + "| ", end='')           



        print ("")


#Only print tab delimited output with results 
if args.stats:
    K = len(lineage1) 
    for index in range(len(dmis)):
        if str(dmis[index][0]).islower() or str(dmis[index][1]).islower():
            miCount += 1
        else:
            diCount += 1
    print(K, len(dmis), diCount, miCount, K * (K - 1) / 2, sep='\t') 






