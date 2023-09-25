#!/usr/env python
'''
script to parse and convert each phosphorylation site from the Excel file to a JSON
representation. Each json can then be the input in a single CLUSTALO alignment to assess
which genomes the phos site it can be found in.


Usage:
python mo_sites_to_json.py <scratch dir for JSON>


'''

import sys
import json
import re
import os
import pandas as pd





outdir = sys.argv[1]
xl = pd.read_excel("lib/Nrc210121_rgm2r.xlsx", sheet_name="PeptideIsoforms")
xl = xl[xl['Modifications in Master Proteins'].str.contains("Phospho") == True]
protein_ids = xl['Master Protein Accessions'] .tolist()
#protein_ids  = ['MGG_01819T0']
positions = xl['Positions in Master Proteins'].tolist()
#positions = ['MGG_01819T0 [15-38]']
mod_info =  xl['Modifications in Master Proteins'].tolist()
#mod_info = ['MGG_01819T0 1xPhospho [S18(90.8)]']
row_nums = xl.index.tolist()
#row_nums = [33]

annot_seq = xl['Annotated Sequence']
mod_pattern = xl['Modification Pattern']


ps_raw = list(zip(row_nums, protein_ids, positions, mod_info,annot_seq,mod_pattern))


class PhosSite:

    def __init__(self, row_num, protein_id, seq_region, mods,an_seq, mod_patt):
        self.row_num = int(row_num)
        self.protein_id = protein_id
        self._get_seq_region(seq_region)
        self._get_mods(mods)
        self._phos_site_id(mods)
        self.annotated_seq = an_seq
        self.mod_pat = mod_patt

    def _get_seq_region(self, seq_region):
        m = re.search("\[(\d+)\-(\d+)\]", seq_region)
        self.seq_region_start = int(m.group(1))
        self.seq_region_end = int(m.group(2))

    def _get_mods(self, mods):
        m = re.search(".*(\d)xPhospho \[(.*)\]", mods)
        self.mod_count = int(m.group(1))
        mods = [re.sub("\(.*\)?","", i).strip() for i in m.group(2).split(";")]
        self.mods = [ [i[0], int(i[1:]) ] for i in mods if i[1:].isnumeric() ] ##removes any mods that are ambiguous (only keeps ones with a number after S/T/Y).

    def _phos_site_id(self,mods):
        self.phos_site_id =  self.protein_id + "-" + "-".join(["".join(map(str, i)) for i in self.mods])

    def to_json(self):
        r = {
            "row_num": self.row_num,
            "phos_site_id": self.phos_site_id,
            "protein_id": self.protein_id,
            "seq_region_start": self.seq_region_start,
            "seq_region_end": self.seq_region_end,
            "mods":self.mods,
            "annotated_seq":self.annotated_seq,
            "mod_pat":self.mod_pat
        }
        return json.dumps(r)


phos_sites = [PhosSite(*i) for i in ps_raw]
for p in phos_sites:
    fname = "{}.json".format(p.phos_site_id)
    out = os.path.join(outdir, fname)
    with open( out, "w") as f:
        f.write(p.to_json())

