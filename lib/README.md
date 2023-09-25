# Pipeline


## Prerequisites

This pipeline needs sonicparanoid. It is not available on the cluster. A singularity def file is available in `scripts/sonicparanoid.def`.
Build that and edit the `scripts/do_sonicparanoid.sh` file to refer to the new image.

## Notes

24-Feb-2023: The sonicparanoid process mmseqs broke on the initial data but is fine on the supplied test data.
May need to clean up fasta files for it.

Done using

```shell
[macleand@TSL-HPC proteomes]$ awk -F'|' '/^>/ {print ">"$2; next} 1' UP000002489_660025.fasta > fo5176.fa
[macleand@TSL-HPC proteomes]$ awk '{print $1}' Saccharomyces_cerevisiae.R64-1-1.pep.all.fa > ScR64.fa
[macleand@TSL-HPC proteomes]$ awk '{print $1}' Magnaporthe_oryzae.MG8.pep.all.fa > Mo8.fa
```

Cleaning fasta header and file names seems to allow sonicparanoid to work