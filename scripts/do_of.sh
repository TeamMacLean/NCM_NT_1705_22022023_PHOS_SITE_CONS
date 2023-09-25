#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

#source orthofinder-2.3.7

#source package /tsl/software/testing/bin/mafft-7.490
source diamond-2.0.14
source orthofinder-2.3.7

#orthofinder -f $1 -o $2 -og -t 32
orthofinder -f /tsl/scratch/macleand/phos_cons/proteomes/  -o /tsl/scratch/macleand/phos_cons/orthof/ -og -t 32