#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""   
See README.md @ https://github.com/clinte14/GeneCoOccurrence
Contact: elg6752@vandals.uidaho.edu or clintelg@gmail.com
"""
# Fixing text so it's not cutoff on the heatmap is on my to do list as well
# Fix log to show actual command line entry, not a list
import sys
import os
import shutil
import pandas as pd
from GeneCoOccurrence.housekeeping import create_folders, parse_flags, save_flags, convert_protID_to_common_names
from GeneCoOccurrence.blast_parse import parse_merge_BLAST, create_pa, clean_BLAST_xml
from GeneCoOccurrence.correlation_calcs import correlation_calcs, calc_mrs
from GeneCoOccurrence.create_visuals import create_heatmap, create_network_map

def main():
    
    # Obtain flag values by parsing keyboard input at BASH command
    flag_values = parse_flags()

    # Create folders for project in path defined by -o flag. Default is current directory.
    create_folders(flag_values['output'])
    # Save/record flags options and parameters to "command.txt" file in 00_settings folder
    save_flags(flag_values)

    if flag_values['input'].split('.')[-1]=='xml':    
        # Parse XML files located in 01_BLAST_results folder.
        hits = parse_merge_BLAST(flag_values)

        # Converts BLAST query protein ID's to common gene names using .csv if -c flag was invoked.
        convert_protID_to_common_names(flag_values, hits)

        #write hits dictionary to Pandas DF, query gene is column and hits are rows. Write dataframe to 02_PA_matrix      
        #folder as PA_matrix.csv.
        create_pa(hits, flag_values)
    
    elif flag_values['input'].split('.')[-1]=='csv':

        pa_dir = os.path.join(flag_values['output'] ,'01_PA_matrix','pa_matrix.csv')
        shutil.copyfile(flag_values['input'],pa_dir)
        print("  --->Copied custom matrix '{}' to '{}'".format(flag_values['input'],pa_dir + '\n'))

    else:
        print("FATAL ERROR. Input must be either BLAST formatted '.xml' file or comma-seperated presence/absence '.csv' file")
        sys.exit(0)

    #calculate Cij (# species common to gene i & j), Rij (Pearson Correlation)
    #and Wij (Partial Correlation)
    Wij_df = correlation_calcs(flag_values)
    network_list = calc_mrs(Wij_df, flag_values)

    create_network_map(flag_values, network_list)
    create_heatmap(Wij_df, flag_values)

if __name__ == "__main__":
    main()