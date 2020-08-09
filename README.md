# *dmiSim.slim* reproducibility documentation <!-- omit in toc -->

- [1. Requirements](#1-requirements)
- [2. Options implemented in *dmiSim.slim* v.2.14](#2-options-implemented-in-dmisimslim-v214)
- [3. Basic demographic simulation as example](#3-basic-demographic-simulation-as-example)
- [4. Simulation of models included in main document](#4-simulation-of-models-included-in-main-document)
  - [4.1. Orr (1995) models](#41-orr-1995-models)
    - [4.1.1. Orr model](#411-orr-model)
    - [4.1.2. Orr-Inf model](#412-orr-inf-model)
    - [4.1.3. Orr-Inf+Rev model](#413-orr-infrev-model)
  - [4.2. Demographic models](#42-demographic-models)
    - [4.2.1. Dem model](#421-dem-model)
    - [4.2.2. Dem+Sel](#422-demsel)
    - [4.2.3. Dem+Rev](#423-demrev)
    - [4.2.4. Dem+Rev+Sel](#424-demrevsel)
- [5. Extraction of DMI information](#5-extraction-of-dmi-information)
- [6. Allele frequencies](#6-allele-frequencies)
- [7. Fitting curves](#7-fitting-curves)
- [8. Simulation of models with selection to show the effect of population size](#8-simulation-of-models-with-selection-to-show-the-effect-of-population-size)
- [9. Mutational-Order analysis](#9-mutational-order-analysis)

## 1. Requirements

Following script and programs can be found in `bin` folder. Parameters for all of these could be consulted using `scriptname -h`

```bash
dmiGenerator.py #Generates lists of theoretical DMIs base on a list of genes
multipleTimeGraphs.py #Uses numHybsByNumofDMIs.json file to extract the number of DMIs in average by hybrid
extractAverageDnMis.py #Uses infoHybs.json file to extract the number of Dobzhansky and Muller incompatibilities in average by hybrid
jsonObjectstoArray.sh #Transform JSON files with objects into one only JSON array
plotAlleleFrequencies.py #Plots allele frequencies obtained from dmiSim when -d freq flag is True
```

Following script are found in `scripts_and_notebooks` folder, these are not developed to run from the terminal, so they do not accept any parameter in command line.

```bash
AICw.R #Fits DMI accumulation curves using a given set of functions
AICw_DDDA.R #Fits DMI accumulation curves discriminating DDI and ADI using a given set of functions
jsonToOrderOfFixtion.ipynb #Quantifies homoplasy in fixation orders
```

## 2. Options implemented in *dmiSim.slim* v.2.14

`dmiSim.slim` can be found in `slim_script` folder.

```text
Usage of dmiSim.slim:
    slim [parameters] dmiSim_v2.14.slim > output.full

Example:
    slim -d N=100000 -d K1=10000 -d K2=10000 -d K3=10000 -d K4=10000 -d freq=T -d "dmi='BA-CA-Cb'" -d "outputPath='freq/'" -s 26737249 dmiSim_v2.14.slim > output.full

Flag                    Default value Explaination
-s int                  [auto]        Seed
-d N=int                [10000]       Number of generations (hybrid generation is N+1)
-d "dmi='str'"          [Ba-CA-CB]    DMIs pairs, must be a string following the format: geneAgeneB-geneXgeneZ, being lowercase non mutated and upper mutated (e.g. 'Ca-Ba-Bc'), single quotes must be included
-d probSup=float        [1.0]         Probability total of been supressed, only use in case DMI pairs do not cause epistatic interaction in populations
-d "outputPath='dir'"   [slimOut/]    Output path relative to the script path or absolute, only for tables, not stadin
-d gl=int               [200]         Number of positions (length) for each gene
-d sl=int               [1000]        Number of positions (length) for each spacer
-d a=int                [3]           Number of genes to activate automatic generation of DMIs. dmiGenerator.py must be in bin and be executable
-d K1=int               [1000]        Population 1 size
-d K2=int               [500]         Population 2 size  
-d K3=int               [500]         Population 3 size
-d K4=int               [500]         Population 4 size (based on migration rate from pop2 and pop3)
-d "note='str'"         [empty]       String to append to outputs file names e.g. if note is 'run1_' the table output will be 'run1_resultTable.json'
-d "uid='str'"          [empty]       String to append to table plots outputs file names
-d verbose=Bool(T/F)    [F]           Show information about each individual
-d history=Bool(T/F)    [F]           Show information about each position mutated by individual
-d hybrid=Bool(T/F)     [F]           Show information about each individual in the hybrid generation
-d totalHybCohort=int   [10]          Create int of cohorts of hybrids during the entire simulation
-d multifiles=Bool(T/F) [F]           Produce multiple files per hybrid population, if F, only produce summary files
-d spacer=Bool(T/F)     [T]           Simulate spacers
-d zoom=Bool(T/F)       [F]           Generate more hybrids at the start of simulation to track fast changes
-d orr=Bool(T/F)        [F]           Modify some variables to emulate Orr simulations
-d rev=Bool(T/F)        [F]           Treat each gene as in a reversible model, if A mutate once A is Derived, if mutate again, become ancestral (even mutations = ancestral, odd mutations = derived)
-d orrNoInf(T/F)        [F]           Emulate a orr simulation but allowing mutate or not any gene many times. Limited by number of generations.
-d orrNoInfCicles=int   [10]          Factor to multiply the number of genes to get the total number of hybrid cohorts (and number of mutation events) Remember in Orr emulations each hybrid cohort is a mutation event
-d sampledHybCohort=int [24]          Only works with orrNoInf models, stablishes number of hybrid cohorts that  will.
-d freq=Bool(T/F)       [F]           Print a JSON file with allele frequencies by subpop and by gene by generation.
```

## 3. Basic demographic simulation as example

The next script is a simple instruction to run `dmiSim.slim` in a local machine using 48 cores. Running a simulation for 10000 generations with Ne = 1000 for all populations, creating 10 hybrid cohorts during the process. Every individual has a genome with 5 genes (of 100 bases lenght) and 6 spacers (of 2000 bases lenght).

```bash
for x in {1..48}
do

experimentName="exampleRun"
inputFolder=~/slimThings
outputFolder=~/slimThings/results/${experimentName}
slimScript=dmiSim_v2.10.slim

S=${RANDOM}${RANDOM}
KT=1000
totalHybCohort=10
parametersToSlim="-d N=10000 -d K1=${KT} -d K2=${KT} -d K3=${KT} -d a=5 -d totalHybCohort=${totalHybCohort} -d probSup=0 -d gl=100 -d sl=2000"
note="note='${experimentName}'"
uid="uid='${S}'"
outfile=${outputFolder}/${S}.full


mkdir -p ${outputFolder}/
slim -s ${S} -d "outputPath='${outputFolder}/'" -d ${note} -d ${uid} ${parametersToSlim} ${inputFolder}/${slimScript} > ${outfile} &
done
```

## 4. Simulation of models included in main document

Following are the contents of bash files to run the code in HPC systems using SLURM. Copy of those files can be found in `sbatchs` folder.

### 4.1. Orr (1995) models

Three different models were simulated under Orr category (see table 1 or suppplementary table 1).

#### 4.1.1. Orr model

No demographic dinamics [instant fixation] with infinite sites.

```bash
#!/bin/bash
#SBATCH --job-name=orr
#SBATCH --output=%x_%A_%a.out
#SBATCH --array=1-100
#SBATCH --time=11:00:00
#SBATCH --account=account
#SBATCH --mem-per-cpu=1GB
#SBATCH --cpus-per-task=1
#SBATCH --workdir=/moto/home/user/logsSlurm/
#SBATCH --mail-user=email@example.com
#SBATCH --mail-type=ALL

#load Slim module
module load SLiM
#load anaconda for python3
module load anaconda/3-2018.12

#constants
inputFolder=~/slimThings
currentArray=${SLURM_ARRAY_TASK_ID}
outputFolder=~/data/resultsForPaper/${SLURM_JOB_NAME}

#define variables and toggles for dmiSim execution
experimentName="orr"
slimScript=dmiSim_v2.10.slim
genes=24
S=${RANDOM}${RANDOM}
parametersToSlim="-d a=${genes} -d probSup=0 -d orr=T"
note="note='${experimentName}'"
uid="uid='${S}'"
outfile=${outputFolder}/${S}.full

mkdir -p ${outputFolder}/
slim -s ${S} -d "outputPath='${outputFolder}/'" -d ${note} -d ${uid} ${parametersToSlim} ${inputFolder}/${slimScript} > ${outfile}
```

#### 4.1.2. Orr-Inf model

No demographic dinamics [instant fixation] with finite sites.

```bash
#!/bin/bash
#SBATCH --job-name=orr-inf
#SBATCH --output=%x_%A_%a.out
#SBATCH --array=1-100
#SBATCH --time=11:00:00
#SBATCH --account=account
#SBATCH --mem-per-cpu=1GB
#SBATCH --cpus-per-task=1
#SBATCH --workdir=/moto/home/user/logsSlurm/
#SBATCH --mail-user=email@example.com
#SBATCH --mail-type=ALL

#load Slim module
module load SLiM
#load anaconda for python3
module load anaconda/3-2018.12

#constants
inputFolder=~/slimThings
currentArray=${SLURM_ARRAY_TASK_ID}
outputFolder=~/data/resultsForPaper/${SLURM_JOB_NAME}

#define variables and toggles for dmiSim execution
experimentName="orr-inf"
slimScript=dmiSim_v2.10.slim
genes=24
S=${RANDOM}${RANDOM}
totalHybCohort=10
parametersToSlim="-d a=${genes} -d probSup=0 -d orrNoInf=T -d orrNoInfCicles=20 -d sampledHybCohort=24"
note="note='${experimentName}'"
uid="uid='${S}'"
outfile=${outputFolder}/${S}.full

mkdir -p ${outputFolder}/
slim -s ${S} -d "outputPath='${outputFolder}/'" -d ${note} -d ${uid} ${parametersToSlim} ${inputFolder}/${slimScript} > ${outfile}
```

#### 4.1.3. Orr-Inf+Rev model

No demographic dinamics [instant fixation] with finite sites and with reversibility.

```bash
#!/bin/bash
#SBATCH --job-name=orr-inf+rev
#SBATCH --output=%x_%A_%a.out
#SBATCH --array=1-100
#SBATCH --time=11:00:00
#SBATCH --account=account
#SBATCH --mem-per-cpu=1GB
#SBATCH --cpus-per-task=1
#SBATCH --workdir=/moto/home/user/logsSlurm/
#SBATCH --mail-user=email@example.com
#SBATCH --mail-type=ALL

#load Slim module
module load SLiM
#load anaconda for python3
module load anaconda/3-2018.12

#constants
inputFolder=~/slimThings
currentArray=${SLURM_ARRAY_TASK_ID}
outputFolder=~/data/resultsForPaper/${SLURM_JOB_NAME}

#define variables and toggles for dmiSim execution
experimentName="orr-inf+rev"
slimScript=dmiSim_v2.10.slim
genes=24
S=${RANDOM}${RANDOM}
totalHybCohort=10
parametersToSlim="-d a=${genes} -d probSup=0 -d orrNoInf=T -d orrNoInfCicles=20 -d rev=T -d sampledHybCohort=24"
note="note='${experimentName}'"
uid="uid='${S}'"
outfile=${outputFolder}/${S}.full

mkdir -p ${outputFolder}/
slim -s ${S} -d "outputPath='${outputFolder}/'" -d ${note} -d ${uid} ${parametersToSlim} ${inputFolder}/${slimScript} > ${outfile}
```

### 4.2. Demographic models

#### 4.2.1. Dem model

Demographic dinamics.

```bash
#!/bin/bash
#SBATCH --job-name=dem
#SBATCH --output=%x_%A_%a.out
#SBATCH --array=1-100
#SBATCH --time=11:00:00
#SBATCH --account=account
#SBATCH --mem-per-cpu=2GB
#SBATCH --cpus-per-task=1
#SBATCH --workdir=/moto/home/user/logsSlurm/
#SBATCH --mail-user=email@example.com
#SBATCH --mail-type=ALL

#load Slim module
module load SLiM
#load anaconda for python3
module load anaconda/3-2018.12

#constants
inputFolder=~/slimThings
currentArray=${SLURM_ARRAY_TASK_ID}
outputFolder=~/data/resultsForPaper/${SLURM_JOB_NAME}

#define variables and toggles for dmiSim execution

experimentName="dem"
slimScript=dmiSim_v2.10.slim
S=${RANDOM}${RANDOM}
KT=100
totalHybCohort=10
parametersToSlim="-d N=250000 -d K1=${KT} -d K2=${KT} -d K3=${KT} -d a=24 -d totalHybCohort=${totalHybCohort} -d probSup=0 -d gl=300 -d sl=2000 -d zoom=T"
note="note='${experimentName}'"
uid="uid='${S}'"
outfile=${outputFolder}/${S}.full

mkdir -p ${outputFolder}/
slim -s ${S} -d "outputPath='${outputFolder}/'" -d ${note} -d ${uid} ${parametersToSlim} ${inputFolder}/${slimScript} > ${outfile}
```

#### 4.2.2. Dem+Sel

Demographic dinamics with multiple incompabitilities (sensu Orr & Turelli, 2001) like selection.

```bash
#!/bin/bash
#SBATCH --job-name=dem+sel
#SBATCH --output=%x_%A_%a.out
#SBATCH --array=1-100
#SBATCH --time=119:00:00
#SBATCH --account=account
#SBATCH --mem-per-cpu=2GB
#SBATCH --cpus-per-task=1
#SBATCH --workdir=/moto/home/user/logsSlurm/
#SBATCH --mail-user=email@example.com
#SBATCH --mail-type=ALL

#load Slim module
module load SLiM
#load anaconda for python3
module load anaconda/3-2018.12

#constants
inputFolder=~/slimThings
currentArray=${SLURM_ARRAY_TASK_ID}
outputFolder=~/data/resultsForPaper/${SLURM_JOB_NAME}

#define variables and toggles for dmiSim execution
experimentName="dem+sel"
slimScript=dmiSim_v2.10.slim
S=${RANDOM}${RANDOM}
KT=100
totalHybCohort=10
parametersToSlim="-d N=500000 -d K1=${KT} -d K2=${KT} -d K3=${KT} -d a=24 -d totalHybCohort=${totalHybCohort} -d probSup=1 -d gl=300 -d sl=2000 -d zoom=T"
note="note='${experimentName}'"
uid="uid='${S}'"
outfile=${outputFolder}/${S}.full

mkdir -p ${outputFolder}/
slim -s ${S} -d "outputPath='${outputFolder}/'" -d ${note} -d ${uid} ${parametersToSlim} ${inputFolder}/${slimScript} > ${outfile}
```

#### 4.2.3. Dem+Rev

Demographic dinamics with reversible states.

```bash
#!/bin/bash
#SBATCH --job-name=dem+rev
#SBATCH --output=%x_%A_%a.out
#SBATCH --array=1-100
#SBATCH --time=11:00:00
#SBATCH --account=account
#SBATCH --mem-per-cpu=2GB
#SBATCH --cpus-per-task=1
#SBATCH --workdir=/moto/home/user/logsSlurm/
#SBATCH --mail-user=email@example.com
#SBATCH --mail-type=ALL

#load Slim module
module load SLiM
#load anaconda for python3
module load anaconda/3-2018.12

#constants
inputFolder=~/slimThings
currentArray=${SLURM_ARRAY_TASK_ID}
outputFolder=~/data/resultsForPaper/${SLURM_JOB_NAME}

#define variables and toggles for dmiSim execution
experimentName="dem+rev"
slimScript=dmiSim_v2.10.slim
S=${RANDOM}${RANDOM}
KT=100
totalHybCohort=10
parametersToSlim="-d N=250000 -d K1=${KT} -d K2=${KT} -d K3=${KT} -d a=24 -d totalHybCohort=${totalHybCohort} -d probSup=0 -d gl=300 -d sl=2000 -d zoom=T -d rev=T"
note="note='${experimentName}'"
uid="uid='${S}'"
outfile=${outputFolder}/${S}.full

mkdir -p ${outputFolder}/
slim -s ${S} -d "outputPath='${outputFolder}/'" -d ${note} -d ${uid} ${parametersToSlim} ${inputFolder}/${slimScript} > ${outfile}
```

#### 4.2.4. Dem+Rev+Sel

Demographic dinamics with reversible states and multiple incompabitilities (sensu Orr & Turelli, 2001) like selection.

```bash
#!/bin/bash
#SBATCH --job-name=dem+rev+sel
#SBATCH --output=%x_%A_%a.out
#SBATCH --array=1-100
#SBATCH --time=119:00:00
#SBATCH --account=account
#SBATCH --mem-per-cpu=2GB
#SBATCH --cpus-per-task=1
#SBATCH --workdir=/moto/home/user/logsSlurm/
#SBATCH --mail-user=email@example.com
#SBATCH --mail-type=ALL

#load Slim module
module load SLiM
#load anaconda for python3
module load anaconda/3-2018.12

#constants
inputFolder=~/slimThings
currentArray=${SLURM_ARRAY_TASK_ID}
outputFolder=~/data/resultsForPaper/${SLURM_JOB_NAME}

#define variables and toggles for dmiSim execution
experimentName="dem+rev+sel"
slimScript=dmiSim_v2.10.slim
S=${RANDOM}${RANDOM}
KT=100
totalHybCohort=10
parametersToSlim="-d N=500000 -d K1=${KT} -d K2=${KT} -d K3=${KT} -d a=24 -d totalHybCohort=${totalHybCohort} -d probSup=1 -d gl=300 -d sl=2000 -d zoom=T -d rev=T"
note="note='${experimentName}'"
uid="uid='${S}'"
outfile=${outputFolder}/${S}.full

mkdir -p ${outputFolder}/
slim -s ${S} -d "outputPath='${outputFolder}/'" -d ${note} -d ${uid} ${parametersToSlim} ${inputFolder}/${slimScript} > ${outfile}
```

## 5. Extraction of DMI information

With JSON files produced by dmiSim, we used the follow instructions in Bash executing python programs (`multipleTimeGraphs.py` and `extractAverageDnMis.py`) and extracting the number of DMIs per hybrid. This step was performed for all previous models. Here is a example of how to run for a particular model (e.g. orr model).

```console
jsonObjectstoArray.sh *.json; multipleTimeGraphs.py *numHybsByNumofDMIs.json; extractAverageDnMis.py *infoHybs.json *resultTable.json
```

Or in batch for a set of folders.  

```console
for d in */; do cd $d; jsonObjectstoArray.sh *.json; multipleTimeGraphs.py *numHybsByNumofDMIs.json; extractAverageDnMis.py *infoHybs.json *resultTable.json; cd ..; done
```

## 6. Allele frequencies

In order to create allele frequencies plots the parameter `-d freq=T` must be enable. `dmiSim.slim` will generate a JSON file (`alleleFreq.json`). With this file the script `plotAlleleFrequencies.py` may be used to create allele frequency plots. For example:

```console
plotAlleleFrequencies.py 3258914890alleleFreq.json -f svg
```

## 7. Fitting curves

We tested what pattern follow each accumulation curve obtained for each model. First visual examination for four different expressions was performed using Past3. To corroborate this estimations R was used, nslML function included minpack.lm library was used to calculate AIC, weights of AIC and final parameters. The script to do it is `AICw.R`. As explained before, this script is not intended to be use as a command line program; in order to change files names or folder paths, it must be done within the script itself.

Also we test the fit for DDI and ADI with the following R script `AICw_DDDA.R`.

## 8. Simulation of models with selection to show the effect of population size

A predefine list of 100 random numbers were used to be sure the conditions of all experiments are the same, with just the predefine DMIs as variable, following bash command was used.

```bash
for i in {1..100}; do echo ${RANDOM}${RANDOM} >> 100randoms.txt; done
```

We ran a set of three experiments for different Ne (100, 1000, and 10000), with 100 independent runs each one. This run uses 75 cores!
In this experiment mutation rate and recombination rate were scaled by a factor of 100 to increase the computational speed. Following bash script was used.

```bash
for y in {0..3}
do

for x in {1..25}
do

#Define input folder
inputFolder=~/slimThings

#Use a list of 100 random numbers, in this way, I use same random numbers for every set of simulations
line=$((${x}+(${y}*25)))
S=$(sed "${line}q;d" 100randoms.txt)

slimScript=dmiSim_v2.13.slim
totalHybCohort=100

experimentName="dem+sel_f100_m5r6_100i_24g_final"
KT=100
parametersToSlim="-d mutRate=1e-5 -d recRate=1e-6 -d N=3500 -d K1=${KT} -d K2=${KT} -d K3=${KT} -d a=24 -d totalHybCohort=${totalHybCohort} -d probSup=1 -d gl=300 -d sl=2000 -d zoom=F"
note="note='${experimentName}'"
uid="uid='${S}'"
outputFolder=~/slimThings/results/${experimentName}
outfile=${outputFolder}/${S}.full
mkdir -p ${outputFolder}/
slim -s ${S} -d "outputPath='${outputFolder}/'" -d ${note} -d ${uid} ${parametersToSlim} ${inputFolder}/${slimScript} > ${outfile} &


experimentName="dem+sel_f100_m5r6_1000i_24g_final"
KT=1000
parametersToSlim="-d mutRate=1e-5 -d recRate=1e-6 -d N=3500 -d K1=${KT} -d K2=${KT} -d K3=${KT} -d a=24 -d totalHybCohort=${totalHybCohort} -d probSup=1 -d gl=300 -d sl=2000 -d zoom=F"
note="note='${experimentName}'"
uid="uid='${S}'"
outputFolder=~/slimThings/results/${experimentName}
outfile=${outputFolder}/${S}.full
mkdir -p ${outputFolder}/
slim -s ${S} -d "outputPath='${outputFolder}/'" -d ${note} -d ${uid} ${parametersToSlim} ${inputFolder}/${slimScript} > ${outfile} &


experimentName="dem+sel_f100_m5r6_10000i_24g_final"
KT=10000
parametersToSlim="-d mutRate=1e-5 -d recRate=1e-6 -d N=3500 -d K1=${KT} -d K2=${KT} -d K3=${KT} -d a=24 -d totalHybCohort=${totalHybCohort} -d probSup=1 -d gl=300 -d sl=2000 -d zoom=F"
note="note='${experimentName}'"
uid="uid='${S}'"
outputFolder=~/slimThings/results/${experimentName}
outfile=${outputFolder}/${S}.full
mkdir -p ${outputFolder}/
slim -s ${S} -d "outputPath='${outputFolder}/'" -d ${note} -d ${uid} ${parametersToSlim} ${inputFolder}/${slimScript} > ${outfile} &

done

wait
echo "Round ${y} finished" >> status.log

done

```

Extraction of DMI information followed this analysis (explained in section 5).

```console
for d in dem+sel_f100_m5r6*/; do cd $d; jsonObjectstoArray.sh *.json; multipleTimeGraphs.py *numHybsByNumofDMIs.json; extractAverageDnMis.py *infoHybs.json *resultTable.json; cd ..; done
```

All CSV files were used as input for the python code wrote in the Jupyter Notebook called `JSON-to-DMIs-2.ipynb`

## 9. Mutational-Order analysis

For quantifying homoplasy in fixation orders we ran three different simulations using the same set of randome seed thatn in section 8. Here is the script to run it using a predefine non-variable set of DMIs for each simulation.

```bash
for y in {0..3}
do

for x in {1..25}
do

#constants
inputFolder=~/slimThings

#Use a list of 100 random numbers, in this way, I use same random numbers for every set of simulations
line=$((${x}+(${y}*25)))
S=$(sed "${line}q;d" 100randoms.txt)


slimScript=dmiSim_v2.14.slim
totalHybCohort=100


experimentName="test_Mutation-Order_Ne100_24genes"
KT=100
parametersToSlim="-d mutRate=1e-5 -d recRate=1e-6 -d N=3500 -d K1=${KT} -d K2=${KT} -d K3=${KT} -d totalHybCohort=${totalHybCohort} -d probSup=1 -d gl=300 -d sl=2000 -d zoom=F -d freq=T"
dmis="dmi='BA-Ca-CB-Da-DB-Dc-EA-Eb-EC-ED-FA-Fb-FC-FD-Fe-GA-Gb-GC-GD-Ge-Gf-Ha-HB-Hc-Hd-HE-HF-HG-IA-Ib-IC-ID-Ie-If-Ig-IH-Ja-JB-Jc-Jd-JE-JF-JG-Jh-JI-Ka-KB-Kc-Kd-KE-KF-KG-Kh-KI-Kj-La-LB-Lc-Ld-LE-LF-LG-Lh-LI-Lj-Lk-MA-Mb-MC-MD-Me-Mf-Mg-MH-Mi-MJ-MK-ML-NA-Nb-NC-ND-Ne-Nf-Ng-NH-Ni-NJ-NK-NL-Nm-OA-Ob-OC-OD-Oe-Of-Og-OH-Oi-OJ-OK-OL-Om-On-Pa-PB-Pc-Pd-PE-PF-PG-Ph-PI-Pj-Pk-Pl-PM-PN-PO-Qa-QB-Qc-Qd-QE-QF-QG-Qh-QI-Qj-Qk-Ql-QM-QN-QO-Qp-Ra-RB-Rc-Rd-RE-RF-RG-Rh-RI-Rj-Rk-Rl-RM-RN-RO-Rp-Rq-Sa-SB-Sc-Sd-SE-SF-SG-Sh-SI-Sj-Sk-Sl-SM-SN-SO-Sp-Sq-Sr-TA-Tb-TC-TD-Te-Tf-Tg-TH-Ti-TJ-TK-TL-Tm-Tn-To-TP-TQ-TR-TS-Ua-UB-Uc-Ud-UE-UF-UG-Uh-UI-Uj-Uk-Ul-UM-UN-UO-Up-Uq-Ur-Us-UT-VA-Vb-VC-VD-Ve-Vf-Vg-VH-Vi-VJ-VK-VL-Vm-Vn-Vo-VP-VQ-VR-VS-Vt-VU-WA-Wb-WC-WD-We-Wf-Wg-WH-Wi-WJ-WK-WL-Wm-Wn-Wo-WP-WQ-WR-WS-Wt-WU-Wv-XA-Xb-XC-XD-Xe-Xf-Xg-XH-Xi-XJ-XK-XL-Xm-Xn-Xo-XP-XQ-XR-XS-Xt-XU-Xv-Xw'"
note="note='${experimentName}'"
uid="uid='${S}'"
outputFolder=~/slimThings/results/${experimentName}
outfile=${outputFolder}/${S}.full
mkdir -p ${outputFolder}/
slim -s ${S} -d "outputPath='${outputFolder}/'" -d ${dmis} -d ${note} -d ${uid} ${parametersToSlim} ${inputFolder}/${slimScript} > ${outfile} &



experimentName="test_Mutation-Order_Ne1000_24genes"
KT=1000
parametersToSlim="-d mutRate=1e-5 -d recRate=1e-6 -d N=3500 -d K1=${KT} -d K2=${KT} -d K3=${KT} -d totalHybCohort=${totalHybCohort} -d probSup=1 -d gl=300 -d sl=2000 -d zoom=F -d freq=T"
dmis="dmi='BA-Ca-CB-Da-DB-Dc-EA-Eb-EC-ED-FA-Fb-FC-FD-Fe-GA-Gb-GC-GD-Ge-Gf-Ha-HB-Hc-Hd-HE-HF-HG-IA-Ib-IC-ID-Ie-If-Ig-IH-Ja-JB-Jc-Jd-JE-JF-JG-Jh-JI-Ka-KB-Kc-Kd-KE-KF-KG-Kh-KI-Kj-La-LB-Lc-Ld-LE-LF-LG-Lh-LI-Lj-Lk-MA-Mb-MC-MD-Me-Mf-Mg-MH-Mi-MJ-MK-ML-NA-Nb-NC-ND-Ne-Nf-Ng-NH-Ni-NJ-NK-NL-Nm-OA-Ob-OC-OD-Oe-Of-Og-OH-Oi-OJ-OK-OL-Om-On-Pa-PB-Pc-Pd-PE-PF-PG-Ph-PI-Pj-Pk-Pl-PM-PN-PO-Qa-QB-Qc-Qd-QE-QF-QG-Qh-QI-Qj-Qk-Ql-QM-QN-QO-Qp-Ra-RB-Rc-Rd-RE-RF-RG-Rh-RI-Rj-Rk-Rl-RM-RN-RO-Rp-Rq-Sa-SB-Sc-Sd-SE-SF-SG-Sh-SI-Sj-Sk-Sl-SM-SN-SO-Sp-Sq-Sr-TA-Tb-TC-TD-Te-Tf-Tg-TH-Ti-TJ-TK-TL-Tm-Tn-To-TP-TQ-TR-TS-Ua-UB-Uc-Ud-UE-UF-UG-Uh-UI-Uj-Uk-Ul-UM-UN-UO-Up-Uq-Ur-Us-UT-VA-Vb-VC-VD-Ve-Vf-Vg-VH-Vi-VJ-VK-VL-Vm-Vn-Vo-VP-VQ-VR-VS-Vt-VU-WA-Wb-WC-WD-We-Wf-Wg-WH-Wi-WJ-WK-WL-Wm-Wn-Wo-WP-WQ-WR-WS-Wt-WU-Wv-XA-Xb-XC-XD-Xe-Xf-Xg-XH-Xi-XJ-XK-XL-Xm-Xn-Xo-XP-XQ-XR-XS-Xt-XU-Xv-Xw'"
note="note='${experimentName}'"
uid="uid='${S}'"
outputFolder=~/slimThings/results/${experimentName}
outfile=${outputFolder}/${S}.full
mkdir -p ${outputFolder}/
slim -s ${S} -d "outputPath='${outputFolder}/'" -d ${dmis} -d ${note} -d ${uid} ${parametersToSlim} ${inputFolder}/${slimScript} > ${outfile} &


done

wait
echo "Round ${y} finished" >> status.log

done

```

Allele frequency files were used as input for the python code wrote in the Jupyter Notebook called `jsonToOrderOfFixtion.ipynb` in order to plot the order of fixation.
