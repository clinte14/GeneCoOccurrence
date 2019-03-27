from housekeeping import create_folders, parse_input, save_flags
from blast_parse import query_BLAST, parse_BLAST, create_pa

def main():
    # Spyder/IDE: set use_IDE to True
    # BASH: set use_IDE to False
    use_IDE = False
    
    # obtain flag values. For BASH, parse input. For Spyder/IDE, use a pre-defined dictionary.
    flag_values = parse_input(use_IDE)

    # create folders for project in path defined by -o "/some/place" flag. Default is the directory program ran from.
    create_folders(flag_values['output'])
    
    # save/record flags options and paramaters to "command.txt" file in 00_settings folder
    save_flags(flag_values)
    
    # performs online NCBIWWW.qblast if 'skipBLAST' flag is set to False (default). Output is XML formatted BLAST 
    # results (one per search query/line in input file) in 01_BLAST_results folder.  
    if flag_values['skipBLAST'] is False:
        query_BLAST(flag_values)
    else:
        print('-----Skipping BLAST search...\n-> Skipped\n')
    
    # parse XML files located in 01_BLAST_results folder. Discard MULTISPECIES hits. Strip remaining BLAST hits to  
    # only include specices name [hit_id] and add to hits dictionary, query gene is key and associated value is list   
    # of strings formatted as species name [hit_id]
    hits = parse_BLAST(flag_values)

    # write hits dictionary to Pandas DF, query gene is column and hits are rows. Write dataframe to 02_PA_matrix      
    # folder as PA_matrix.csv.
    create_pa(hits, flag_values)

    print('debug stop point')


if __name__ == "__main__":
    main()
    
