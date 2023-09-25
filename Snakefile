import os


'''
Running this pipeline:
1. Pregenerate the json files describing the phossites from the excel file with
`mo_sites_to_json.py`. Put the json files in the jsons scratch folder

2. Prepare the fasta proteome files: 
    * rename to a very short key (see the lib/all_proteomes.csv file for ideas (this is important, some file names are hard
    coded in scripts. 
    * reformat each sequence name to its briefest identifier
    * put the files into the folder `scratch/proteomes`. 


3. Run the orthofinder step manually:
    Though there is an orthofinder rule it will not run because of problems in how orthofinder deals with its output 
    directories. To get round it you need to run the sbatch scripts/do_of.sh before running this pipeline. Then update the 
    variable ortho_info with the new output file.
'''



scratch = "/tsl/scratch/macleand/phos_cons/"
#sonic_out = scratch + "sonicrun/"
#sonic_bin = "/hpc-home/macleand/singularity_images/"
proteomes = scratch + "proteomes/"
tmp = scratch + "tmp/"
jsons = scratch + "jsons/"
ortho_group_hits = scratch + "ortho_group_hits/"
#ortho_info = sonic_out + "runs/snakerun2/ortholog_groups/flat.ortholog_groups.tsv"

orthof_out = scratch + "orthof/"
#ortho_info = orthof_out + "Results_Mar29/Orthogroups/Orthogroups.tsv"
ortho_info = orthof_out + "Results_May23/Orthogroups/Orthogroups.tsv"

info_file = "lib/all_proteomes.csv"
phos_sites = [os.path.splitext(os.path.basename(p))[0] for p in os.listdir(jsons)]

rule all:
    input:
        test="results/merged_test_sites.csv",
        control="results/merged_control_sites.csv",
        rooted_species_tree="results/SpeciesTreeRooted.txt",
        node_labels="results/SpeciesTree_rooted_node_labels.txt"

rule trees:
    input: orthof_out
    output:
        rooted_species_tree="results/Species_Tree/SpeciesTree_rooted.txt",
        node_labels="results/SpeciesTree_rooted_node_labels.txt"
    params:
        mem='128G',
        queue="tsl-medium",
        outdir="ortho_info"

    shell:
        "source diamond-2.0.14; source orthofinder-2.3.7; orthofinder -fg /tsl/scratch/macleand/phos_cons/orthof/{params.outdir} -t 32 && cp {params.outdir}/Species_Tree/SpeciesTree_rooted.txt {output.rooted_species_tree} && cp {params.outdir}/Species_Tree/SpeciesTree_rooted_node_labels.txt {output.node_labels}"
rule make_orthogroups:
    input: proteomes
    output: ortho_info
    threads: 32
    params:
        mem='128G',
        queue="tsl-medium",
        outdir=orthof_out
    shell:
        "bash scripts/do_of.sh {input} {params.outdir}"

rule find_site_in_orthogroup:
    input:
        phos_site= jsons + "{phos_site}.json",
        ortho_data=ortho_info
    output: ortho_group_hits + "{phos_site}_test_hit_info.csv"
    params:
        mem='8G',
        queue='tsl-short',
        proteomes=proteomes,
        tmpdir=tmp
    shell:
        "bash scripts/do_find_site.sh {params.proteomes} {input.phos_site} {input.ortho_data} {params.tmpdir} TEST {output}"

rule merge_test_sites:
    input:
        ortho_group_hits=expand(ortho_group_hits + "{ps}_test_hit_info.csv", ps=phos_sites)
    output: "results/merged_test_sites.csv"
    params:
        mem='64G',
        queue='tsl-short'
    shell:
        "awk 'FNR>1' {ortho_group_hits}*_test_hit_info.csv > results/tmpout && cat lib/header.csv results/tmpout > results/merged_test_sites.csv && rm results/tmpout"

rule rename_fasta:
    '''uses the proteomes.csv file to change the name of the proteome fasta to 
    a shorter version'''
    input: info_file
    params:
        mem = '8G',
        queue = 'tsl-short'
    shell:
        "python scripts/rename_files.py {input}"

rule control_by_reverse:
    '''assess the false positive against a random background by reversing the protein sequence in the alignment'''
    input:
        phos_site= jsons + "{phos_site}.json",
        ortho_data=ortho_info
    output: ortho_group_hits + "{phos_site}_control_hit_info.csv"
    params:
        mem='8G',
        queue='tsl-short',
        proteomes=proteomes,
        tmpdir=tmp
    shell:
        "bash scripts/do_find_site.sh {params.proteomes} {input.phos_site} {input.ortho_data} {params.tmpdir} CONTROL {output}"

rule merge_control_sites:
    input:
        ortho_group_hits=expand(ortho_group_hits + "{ps}_control_hit_info.csv", ps=phos_sites)
    output: "results/merged_control_sites.csv"
    params:
        mem='64G',
        queue='tsl-short'
    shell:
        "awk 'FNR>1' {ortho_group_hits}*_control_hit_info.csv > results/ctrl_tmpout && cat lib/header.csv results/ctrl_tmpout > results/merged_control_sites.csv && rm results/ctrl_tmpout"