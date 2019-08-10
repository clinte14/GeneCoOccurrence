#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
  ============================================================================
  
File Description:
    Obtain gene modules associated with Kegg species ID
    
Input: 
    Panda's dataframe in .csv format. Filename is defined by NCBI_GENES_FROM_MODULES in 
    find_modules_in_species.py. columns=['k_gene_id', 'module', 'k_species_id','k_species_name', 'aaseq'].
    
Output: 

    
Author: Clinton Elg 
Contact: elg6752@vandals.uidaho.edu or clintelg@gmail.com

  ===========================================================================
"""

import numpy as np
import pandas as pd
from Bio.Blast.Applications import NcbiblastpCommandline
import os.path
from Bio.KEGG import REST

GENE_SEQS_FROM_MODULES="~/correlogy/distance_formula_evaluations/varied_taxa/genes_modules_aaseq_varied_short.csv"

df = pd.read_csv(GENE_SEQS_FROM_MODULES, dtype=object)

for counter, seq in df.aaseq.items():
#    bash_cmd="awk '/start/{getline;print;}' FILENAME"
    print("Currently BLASTing {}...".format(df.k_gene_id.at[counter]))
    bash_cmd="time blastp -db nr_v5 -num_threads 6 -taxidlist LIST.txt -query INPUT.fasta -evalue 0.0001 -out {}.xml".format(df.k_gene_id.at[counter])
    print(bash_cmd)
    blast_handle = 
    out_handle = open(flag_values['output'] + "/" + "01_BLAST_results" + "/" + i + "_BLAST_results.xml", "w")
    out_handle.write(result_handle.read())
    result_handle.close()  