#!/bin/bash
# -*- coding: utf-8 -*-
"""
  ============================================================================
  
File Description:
    Obtain gene modules associated with Kegg species ID
    
Input: 
    Panda's dataframe in .csv format. Filename is defined by NCBI_GENES_FROM_MODULES in 
    find_modules_in_species.py. columns=['k_gene_id', 'module', 'k_species_id','k_species_name', 'aaseq'].
    
    taxid_2_subtree_species_only.txt 
    Generated Aug 12 2019 from https://www.ncbi.nlm.nih.gov/taxonomy/advanced
    using search '(((txid2[Subtree]) AND species[Rank])) '
Output: 

    
Author: Clinton Elg 
Contact: elg6752@vandals.uidaho.edu or clintelg@gmail.com

  ===========================================================================
"""
import numpy as np
import pandas as pd
from Bio.Blast.Applications import NcbiblastpCommandline
import os
import subprocess
import re
import time
from Bio.KEGG import REST


GENE_SEQS_FROM_MODULES="~/correlogy/distance_formula_evaluations/test/genes_modules_aaseq_test.csv"

df = pd.read_csv(GENE_SEQS_FROM_MODULES, dtype=object)

for counter, seq in df.aaseq.items():
    #Open file with taxid inclusion list, remove the taxid of current species being blasted
    with open ('/home/clint/correlogy/distance_formula_evaluations/taxid_2_subtree_species_only.txt', 'r' ) as f:
        original_taxid_list = f.read()
        current_taxid_list= re.sub(r'^{}\n'.format(df.ncbi_tax_id[counter]), r'', original_taxid_list, flags = re.M)
        
    with open ('/home/clint/correlogy/distance_formula_evaluations/current_taxid_list.txt', 'w+' ) as f:
        f.write(current_taxid_list)
    
    with open ('/home/clint/correlogy/distance_formula_evaluations/current_aaseq.fasta', 'w+' ) as f:
        f.write(df.aaseq[counter])
    
    output_filename=re.sub(r'[^\w\s]','',df.k_gene_id[counter]+"_"+ df.module[counter])
    
    print("Current aa seq: {}".format(df.aaseq[counter]))
    print("Current output filename: {}".format(output_filename))
    print("Currently BLASTing '{}' with taxid '{}' excluded from search...".format(df.k_gene_id[counter], df.ncbi_tax_id[counter]))
    
    bash_cmd="/bin/bash -c \"cd ~/correlogy/distance_formula_evaluations && blastp -db nr_v5 -num_threads 8 -taxidlist current_taxid_list.txt \
    -query current_aaseq.fasta  -outfmt 5 -evalue 0.0001 -out ./test/01_BLAST_results/{} && rm aa_seq.fasta\" ".format(output_filename+".xml")
#    print(bash_cmd)
    os.system(bash_cmd)
#    print("{} / {} completed".format((counter+1, len(df.k_gene_id))))
    