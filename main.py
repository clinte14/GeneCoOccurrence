#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""   
See README.md @ https://github.com/clinte14/correlogy
Contact: elg6752@vandals.uidaho.edu or clintelg@gmail.com
"""
import pandas as pd
from housekeeping import create_folders, parse_flags, save_flags, convert_protID_to_common_names
from blast_parse import parse_merge_BLAST, create_pa, clean_BLAST_xml
from correlation_calcs import correlation_calcs, calc_mrs
from create_visuals import create_heatmap, create_network_map

def main():
    
    # Obtain flag values by parsing keyboard input at BASH command
    flag_values = parse_flags()

    # Create folders for project in path defined by -o flag. Default is current directory.
    create_folders(flag_values['output'])
    
    # Save/record flags options and paramaters to "command.txt" file in 00_settings folder
    save_flags(flag_values)

    '''
    Bug workaround for NCBIWWW. Occasionally NCBI BLAST online writes 'CREATE_VIEW' into 
    the BLAST .xml output with no associated '<>' tags, and this causes xml parsing in 
    Biopython to fail. See: https://github.com/biopython/biopython/issues/3342
    '''
    clean_BLAST_xml(flag_values)
    
    # Parse XML files located in 01_BLAST_results folder.
    hits = parse_merge_BLAST(flag_values)

    
    # Converts BLAST query protein ID's to common gene names using .csv if -c flag was invoked.
    convert_protID_to_common_names(flag_values, hits)

    #write hits dictionary to Pandas DF, query gene is column and hits are rows. Write dataframe to 02_PA_matrix      
    #folder as PA_matrix.csv.
    create_pa(hits, flag_values)
    

    #calculate Cij (# species common to gene i & j), Rij (Pearson Correlation)
    #and Wij (Partial Correlation)
    Wij_df= correlation_calcs(flag_values)
    network_list = calc_mrs(Wij_df, flag_values)

    create_network_map(Wij_df, flag_values, network_list)
    create_heatmap(Wij_df, flag_values)

if __name__ == "__main__":
    main()

'''
To do:
In command.txt add full keyboard input command
In command.txt clean up writing of input, output, -c file
Copy .xml input file to 01 folder
For two nodes that point to each other, only show one correlog score
In MRS visuals: Add color that corresponds with heatmap colors/ improve edge thickness
Clean up notes to PEP-8 standard
Add help info for each function
Standardize output text
Print statements relating to graphical MRS creation/generation
'''