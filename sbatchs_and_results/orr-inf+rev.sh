#!/bin/bash
#SBATCH --job-name=orr-inf+rev
#SBATCH --output=%x_%A_%a.out
#SBATCH --array=1-100
#SBATCH --time=11:00:00
#SBATCH --account=eaton
#SBATCH --mem-per-cpu=1GB
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