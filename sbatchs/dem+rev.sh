#!/bin/bash
#SBATCH --job-name=dem+rev
#SBATCH --output=%x_%A_%a.out
#SBATCH --array=1-100
#SBATCH --time=11:00:00
#SBATCH --account=eaton
#SBATCH --mem-per-cpu=2GB
#SBATCH --cpus-per-task=1
#SBATCH --workdir=/moto/home/cam2384/logsSlurm/
#SBATCH --mail-user=cam2384@columbia.edu
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