"""   
See README.md @ https://github.com/clinte14/correlogy
See coding documentation & TO DO list: @ https://docs.google.com/document/d/1jad90RfWmCfFbAY6HheRkbyQ1qc92bQbZlQKdDB_9KE/edit?usp=sharing
Contact @ elg6752@vandals.uidaho.edu or clintelg@gmail.com
"""
from housekeeping import create_folders, parse_input, save_flags
from blast_parse import parse_merge_BLAST, create_pa
from correlation_calcs import correlation_calcs

def main():
    
    # obtain flag values in BASH by parsing input.
    flag_values = parse_input()

    # create folders for project in path defined by -o "/some/place" flag. Default is the directory program ran from.
    create_folders(flag_values['output'])
    
    save_flags(flag_values)
    # save/record flags options and paramaters to "command.txt" file in 00_settings folder
    
    # performs online NCBIWWW.qblast if 'skipblast' flag is set to False (default). Output is XML formatted BLAST 
    # results (one per search query/line in input file) in 01_BLAST_results folder.  
    
    # parse XML files located in 01_BLAST_results folder. Discard MULTISPECIES hits. Strip remaining BLAST hits to  
    # only include specices name [hit_id] and add to hits dictionary, query gene is key and associated value is list   
    # of strings formatted as species name [hit_id]
    hits = parse_merge_BLAST(flag_values)

    # write hits dictionary to Pandas DF, query gene is column and hits are rows. Write dataframe to 02_PA_matrix      
    # folder as PA_matrix.csv.
    create_pa(hits, flag_values)
    
    # calculate Cij (# species common to gene i & j), Rij (Pearson Correlation), and Wij (Partial Correlation, ie Correlogy)
    correlation_calcs(flag_values)


if __name__ == "__main__":
    main()