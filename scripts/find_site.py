#!/usr/bin/env python
import os
import sys
from collections import defaultdict
import json
import subprocess
import tempfile
from Bio import SeqIO, Seq

'''
find_site.py <proteome fasta directory> <phos_site_info.json> <orthogroups flatfile> <tempdir> 

This script depends on blastp.



The script loads in proteome sequences for Mo and related organisms, and a json file describing a phos site in a known 
Mo peptide. The file names should correspond to the key in the file lib/proteomes.csv
 
The script loads in a list of precomputed orthogroups from orthofinder. These must have been computed on the filenames
in lib/all_proteomes.csv 

The orthologues of the source protein of the phos site are loaded and the phospho peptide sequence is then aligned 
against each in turn using blastp . If any matches are found according to BLAST defaults, 
the best hsp (by bitscore) is retained. The hsp in the orthologue is extended to match the full range of the length of the
phosphopeptide sequence and if the phosphosite lies in the range the site is retained, The hsp of the orthologue is then checked
to see whether the corresponding residue in the orthologue has an exact match to the phos site residue in the 
Mo protein.

'''



pfiles = [os.path.join(sys.argv[1], f) for f in os.listdir(sys.argv[1])]
json_file = sys.argv[2]
ortho_file = sys.argv[3]
tmpdir = sys.argv[4]

CONTROL_RUN = False
if sys.argv[5] == 'CONTROL':
    CONTROL_RUN = True

MAGNAPORTHE_PROTEOME_TAG = "Mo8"

def make_sequence_list(fas):
    '''
    loads fasta sequences to dict, {'id' : Bio::SeqRecord}
    :param fas: fasta file
    :return: dict
    '''
    r = {}
    for f in fas:
        for record in SeqIO.parse(f, "fasta"):
            r[record.id] = record
    return r


def make_ortho_list(orthogroups):
    '''
    loads orthologue ids for Mo proteins
    :param orthogroups: sonicparanoid flat orthogroups file - latterly orthofinder orthogroups file
    :return: dict of dicts {'mo_protein_id : {'other_species_key': [other_species_protein_id1, other_species_protein_id2 ] }
    '''
    r = defaultdict(defaultdict)
    with open(orthogroups, "r") as f:
        col_headers = None
        mg_col = None
        other_cols = None
        for line in f:
            line = line.rstrip()
            a = line.split("\t")
  #          print(len(a), file=sys.stderr)
  #          print(a, file=sys.stderr)
            if line.startswith("Orthogroup"):
                col_headers = a
                mg_col = a.index(MAGNAPORTHE_PROTEOME_TAG)
                other_cols = [i for i in range(len(a)) if not i in [0, mg_col] ]
            else:
                if len(a) == len(col_headers) and len(a[mg_col]):
                    for m in a[mg_col].split(","): # if
                        for g in other_cols:
                            #print(other_cols, file=sys.stderr)
                            #print(f"g is {g}", file=sys.stderr)
                            try:
                                if len(a[g]):
                                    r[m][col_headers[g]] = [i.strip(" ") for i in a[g].split(",")]
                                else:
                                    r[m][col_headers[g]] = None  # a[g].split(",")
                            except IndexError:
                                print("Some sort of out of range error occurred", file=sys.stderr)
                                print(f"col_headers is : {col_headers}", file=sys.stderr)
                                print(f"g is : {g}", file=sys.stderr)
                                print(f"len(a) is : {len(a)}", file=sys.stderr)
    #print(r, file=sys.stderr)
    #sys.exit()
    return r


def get_mg_subseq(ps, seq):
    '''
    get a subseq SeqIO record of the phospho site
    :param ps: the loaded phosphosite json entry
    :param seq: the full sequence of the Mo protein
    :return:
    '''
    return seq[ps['seq_region_start']:ps['seq_region_end']]

def parse_aln(file):
    '''
    gets hsp information from a BLASTP of two sequences
    :param file: the tempfile from a local BLAST subprocess
    :return: a dict of BLAST stats
    '''
    #print(file)
    with open(file, "r") as f:
        b = f.read().rstrip().split("\t")
        if len(b) == 12:
            #print(b)
            nm = ['query', 'subject', 'ident', 'aln_length', 'mismatches', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bits']
            return {x[0]:x[1] for x in zip(nm,b)}
        else:
            return None


def do_aln(mg_sub, o, dir=".", control_run=False):
    '''
    do a BLASTP of two sequences, creates temp files for sequences and spawns a subprocess for the BLAST
    :param mg_sub: the SeqRecord of the mo protein containing the phos_site
    :param o: the SeqRecord of the orthologue
    :param dir: temp directory to write files to
    :return: dict of BLAST stats or None in case of no hits.
    '''
    lib, query, out = [ tempfile.NamedTemporaryFile(mode="w", delete=False, dir=dir).name for i in range(3) ]
    try:
        query_fh = open(query, "w")
        lib_fh = open(lib, "w")
        if control_run:
            o = Seq(str(o)[::-1])
        SeqIO.write(mg_sub, query_fh, "fasta")
        SeqIO.write(o, lib_fh, "fasta")
    #finally:
        query_fh.flush()
        query_fh.close
        lib_fh.flush()
        lib_fh.close

    #try:
        cmd = "blastp -query {} -subject {} -out {} -outfmt 6".format(query, lib, out)
        subprocess.call(cmd, shell=True)

        return parse_aln(out)
    except:
        return None
    finally:
        for f in [lib, query, out]:
            os.remove(f)

def get_residue_match(ps, aln, mg_sub, ortho):
    '''
    works out whether the phos site residue has identity in the Mo and target sequence
    extends the hit in the orthologue sequence to the full length of the Mo region that
    initially contained the phos site
    :param ps: phos site data from json
    :param aln: alignment result stats
    :param mg_sub: SeqRecord for Mo phosphosite region
    :param ortho: SeqRecord for orthologue
    :return:
    '''


    full_start = int(aln['sstart']) - int(aln['qstart'])
    full_end = int(aln['send']) + (len(mg_sub.seq) - int(aln['qend']))

    orth_seq = ortho[full_start:full_end]


    result = {
        'num_sites_to_find': str(len(ps['mods'])),
        'found_count':0,
        'found_sites': list(),
        'mg_seq_in_aln': str(mg_sub[int(aln['qstart']):int(aln['qend'])].seq),
        'orth_seq_in_aln': str(orth_seq.seq),
        'orth_hit_id': ortho.id,
        'annotated_seq': ps['annotated_seq'],
        'mod_pattern': ps['mod_pat'],
        'ps_site_in_seq_region': list()
    }
    for ps_res, ps_site in ps['mods']:
        #print("looking for: {} {}".format(ps_res, ps_site))

        ps_site_in_sub = (ps_site - ps['seq_region_start'] ) - 1
        result['ps_site_in_seq_region'].append(ps_site_in_sub)
        #print("looking for residue:{}".format(ps_site_in_sub))
        #print("identity: {}".format(mg_sub[ps_site_in_sub]))
        #print("orth_seq residue at pos: {}".format(orth_seq[ps_site_in_sub]))
        #print(mg_sub)
        #print(orth_seq)
    try:
        if mg_sub[ps_site_in_sub] == orth_seq[ps_site_in_sub]:
            result['found_count'] += 1
            result['found_sites'].append([ps_res, str(ps_site)])
    except IndexError: #something went wrong with counting of site. Let data pass to output for checking
        pass

    if len(result['found_sites']) > 0:
        result['found_sites'] = ";".join([":".join(a) for a in result['found_sites']])
    else:
        result['found_sites'] = "NA"

    result['found_count'] = str(result['found_count'])
    result['ps_site_in_seq_region'] = ";".join(str(result['ps_site_in_seq_region']))
    return result

def has_site_in_aln(ps, aln):

    real_start_q = int(aln['qstart']) + int(ps['seq_region_start'])
    real_end_q = int(aln['qend']) + int(ps['seq_region_end'])
    for r, s in ps['mods']:
        if s >= real_start_q and s <= real_end_q:
            #print("has site in aln {} {} {}".format(real_start_q, real_end_q, s))
            return True
    return False

sequences = make_sequence_list(pfiles)
ortho_list = make_ortho_list(ortho_file)
phos_site = json.load(open(json_file))



header = ",".join(["phos_site_id", "mo_protein_id", "species_compared", "has_ortholog", "best_hit_ortholog",
                   "has_hit_in_ortholog", "best_hit_ortholog_score", "best_hit_ortholog_identity", "num_p_sites_in_mo",
                   "num_sites_matched_in_ortholog", "sites_found_in_ortho","mo_peptide","mo_seq_in_hit", "orth_seq_in_hit", "annotated_seq", "mod_pattern", "ps_site_in_seq_region"
                   ])
print(header)
for species in ortho_list[phos_site['protein_id']]:

    #print("whole mg seq is", sequences[phos_site['protein_id']].seq)
    orthos = ortho_list[phos_site['protein_id']][species]  #
    mg_sub = get_mg_subseq(phos_site, sequences[phos_site['protein_id']])
    #print("expected region is: {} {}".format(phos_site['seq_region_start'], phos_site['seq_region_end']))
    #print("extracted subseq is: {}".format(mg_sub.seq))

#    for i in orthos:
#        i = i.lstrip(" ")
#        print(f"i is : {i}", file=sys.stderr)
#        if not i in sequences:
#            print(f"not found {i}", file=sys.stderr)
#            sys.exit()
    if ortho_list[phos_site['protein_id']][species] is not None:
        #print("has ortho in {}".format(species))
        aln_results = None
        for i in orthos:
            if i in sequences:
                aln_pairs = [[mg_sub, sequences[i]] for i in orthos] #_list[phos_site['protein_id']][species]]  # seqrecord objects
                aln_results = [do_aln(ap[0], ap[1], dir=tmpdir, control_run=CONTROL_RUN) for ap in aln_pairs] #pairwise alignments
                aln_results = [a for a in aln_results if a is not None] #remove None alignments
                aln_results = [a for a in aln_results if has_site_in_aln(phos_site, a)]
            else:
                print(f"Could not find: {i} in sequences... ", file=sys.stderr)

        #print(aln_results)
        if aln_results:

            aln_scores = [a['bits'] for a in aln_results] #scores
            best_aln = aln_results[aln_scores.index(max(aln_scores))]
            best_score = best_aln['bits']
            identity = best_aln['ident']
            best_ortho = sequences[best_aln['subject']]
            #print("best_ortho seq is: {}".format(best_ortho.seq))

            residue_match = get_residue_match(phos_site, best_aln, mg_sub, best_ortho)
            out_str = ",".join([
                  phos_site['phos_site_id'], phos_site['protein_id'], species, "TRUE", residue_match['orth_hit_id'],
                  "TRUE", best_score, identity,
                  residue_match['num_sites_to_find'], residue_match['found_count'],
                  residue_match['found_sites'], str(mg_sub.seq), residue_match['mg_seq_in_aln'],
                  residue_match['orth_seq_in_aln'], phos_site['annotated_seq'], phos_site['mod_pat'], str(residue_match['ps_site_in_seq_region'])
                 ])
            print(out_str)

        else: #no hit within found ortholog
            ps_site_in_seq_region = ";".join([ str(s - phos_site['seq_region_start']) for i,s in phos_site['mods']])
            out_str = ",".join(
                [phos_site['phos_site_id'], phos_site['protein_id'], species, "TRUE", "NA", "FALSE", 'NA', 'NA', 'NA', 'NA', 'NA', "NA", 'NA', 'NA', phos_site['annotated_seq'], phos_site['mod_pat'], ps_site_in_seq_region ]#.extend(["NA" for i in range(7)])
            )
            print(out_str)
    else: #no ortholog found
        ps_site_in_seq_region = ";".join([str(s - phos_site['seq_region_start']) for i, s in phos_site['mods']])
        out_str = ",".join(
            [phos_site['phos_site_id'], phos_site['protein_id'], species, "FALSE",'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', "NA", 'NA', "NA", phos_site['annotated_seq'], phos_site['mod_pat'], ps_site_in_seq_region] #.extend(["NA" for i in range(9)])
        )
        print(out_str)
