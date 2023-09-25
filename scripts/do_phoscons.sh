#!/bin/bash

if [ -z "$1" ]
then
   sbatch -J phoscons \
    -o phoscons.log \
    --wrap="source snakemake-5.5.3; snakemake -s Snakefile all --cluster \"sbatch  --output=/dev/null --error=/dev/null --partition={params.queue} -c {threads} --mem={params.mem}\" --jobs 80 --latency-wait 60"
elif [ $1 = 'unlock' ]
then
    sbatch -J unlock \
        -o phoscons.log \
        --wrap="source snakemake-6.10.0_gk; snakemake -s Snakefile --unlock" \
        --partition="tsl-short" \
        --mem="16G"
elif [ $1 = "dryrun" ]
then
    sbatch -J dryrun \
    -o phoscons.log \
    --wrap="source snakemake-6.10.0_gk; snakemake -s Snakefile -n" \
    --partition="tsl-short" \
    --mem="16G"
elif [ $1 = "rule" ]
then
    sbatch -J phoscons \
    -o phoscons.log \
    --wrap="source snakemake-5.5.3; snakemake -s Snakefile ${2} --cluster \"sbatch --partition={params.queue} -c {threads} --mem={params.mem}\" --jobs 80 --latency-wait 60"
elif [ $1 = "dag" ]
then
    sbatch -J dag \
    -o phoscons.log \
    --wrap="source snakemake-6.10.0_gk; snakemake -s Snakefile  > phoscons.dot" \
    --partition="tsl-short" \
    --mem="16G"
elif [ $1 = '-h' ]
then
  usage
fi
