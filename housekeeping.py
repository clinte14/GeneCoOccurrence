import os.path, sys
import argparse, datetime, json

def create_folders(dir_path):
    folders = ['00_Settings','01_BLAST_results', '02_PA_matrix', '03_Correlation_calcs', '04_Correlog_values']
    for i in folders:
        
        if not os.path.exists(dir_path + "/" + i):
            os.mkdir(dir_path + "/" + i)
            print("Checking directory '{}'... ".format(dir_path + "/" + i))
            print('  --->Directory Created')
            
        else:
            print("Checking directory '{}'... ".format(dir_path + "/" + i))
            print('  --->Directory Already Exists')
    print('')
    
# obtain flag values. For BASH, parse input. For Spyder/IDE, use a pre-defined dictionary.    
def parse_input(use_IDE):
    dir_path = os.path.dirname(os.path.realpath('__file__'))   
    
    # use flag values pre-defined in dictionary, allows program to run in Spyder/IDE 
    if use_IDE is True:
        arg_dict= {'input':"input.txt", 'output':dir_path, 'evalue':"0.0001", 'entrez':'(txid2 [ORGN]) NOT (txid666 [ORGN])', 'skipblast':True}    
    
    # capture flag options from command line for BASH usage    
    elif use_IDE is False:
        parser = argparse.ArgumentParser(description="Correlation of Gene Homologues", prog='correlog')
        
        parser.add_argument("-i", "--input", required=True, help="List of genes of interest (one gene per line, UTF-8 encoded, Unix LF Newlines). REQUIRED.")
        parser.add_argument("-o", "--output", default=dir_path, help="Desired directory of output file. Default=Current Directory.")
        parser.add_argument("--evalue", default = 0.0001, help="Desired E-Value Cutoff for BLAST Search of Homologous Genes. Default=0.0001.")
        parser.add_argument("--entrez", default = '', help="ENTRTEZ filtering string. Default=empty string.")
    #    parser.add_argument("-e", "--exlude", help="File with tax ID to EXCLUDE from BLAST search (one ID per line, UTF-8 encoded")
    #    parser.add_argument("-l", "--local", default="false", help="Local or Online NCBI BLAST search?")
        parser.add_argument("-s", "--skipblast", default = False, help="Do Not BLAST: Bring Your Own Blast Files (BYOB). Default=False.")
    
        args = parser.parse_args()
        
        arg_dict=vars(args)
    
    # throw error & quit if use_IDE != True or False
    else:
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print(" 'use_IDE value' must be set to True or False; quitting now")
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        sys.exit(0)

    return arg_dict

# save flags options and paramaters to "command.txt" file in 00_Settings folder
def save_flags(flags):
    with open(flags['output'] + "/" + "00_Settings" + "/"'command.txt', 'w') as file:
        currentDT = datetime.datetime.now()
        file.write("Executed the following in " + flags['output'] + " at " + str(currentDT) + ":\n" + json.dumps(flags))
    print("Writing flags to 'command.txt' in '{}'... ".format(flags['output'] + "/" + "00_settings"))
    print("  --->Flags written\n")