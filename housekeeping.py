import os.path
import argparse, datetime, json
import pandas as pd
import sys
def create_folders(dir_path):
    dirs = ['00_pre_process',
            '01_PA_matrix', 
            '02_correlation', 
            '03_visual_output']

    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
        print("Checking directory '{}'... ".format(dir_path))
        print('  --->Directory Created')
    else:
        print("Checking directory '{}'... ".format(dir_path))
        print('  --->Directory Already Exists')

    for current_dir in dirs:
        current_path = os.path.join(dir_path,current_dir)
        if not os.path.exists(current_path):
            os.mkdir(current_path)
            print("Checking directory '{}'... ".format(current_path))
            print('  --->Directory Created')
            
        else:
            print("Checking directory '{}'... ".format(current_path))
            print('  --->Directory Already Exists')
    print('')
    
def parse_flags():
    dir_path = os.path.dirname(os.path.realpath('__file__'))  
    
    parser = argparse.ArgumentParser(description="GeneCoOccurence v0.1. Please see https://github.com/clinte14/GeneCoOccurrence", prog='gco')

    parser.add_argument("-i", "--input", 
        required=True, 
        help="REQUIRED. Input file with either BLAST results in .xml format OR a comma-seperated presence/absence .csv file.")

    parser.add_argument("-o", "--output", 
        default=dir_path,
        help="Desired directory of output file. Default=Current Directory.")

    parser.add_argument("-c", "--common_name",
        help="A .csv file that allows conversion of BLAST query protein ID's \
        to common gene names.")

    args = parser.parse_args()
    arg_dict=vars(args)

    # Convert relative paths (eg './' in BASH) to absolute paths
    if arg_dict['input']!=None:
        arg_dict['input'] = os.path.abspath(arg_dict['input'])

    if arg_dict['output']!=None:
        arg_dict['output'] = os.path.abspath(arg_dict['output'])  

    if arg_dict['common_name']!=None:
        arg_dict['common_name'] = os.path.abspath(arg_dict['common_name'])

    print("dir_path: " + str(dir_path))
    arg_dict['command'] = " ".join(sys.argv)

    return arg_dict


# save flags options and paramaters to "log.txt" file in 00_Settings folder
def save_flags(flags):
    with open(os.path.join(flags['output'],"00_pre_process",'log.txt'), 'w') as file:

        currentDT = datetime.datetime.now()

        file.write("Executed the following in " + flags['output'] + " at " \
            + str(currentDT) + ":\n\n" + flags['command'])

    print("Writing flags to 'log.txt' in '{}'... ".\
        format(os.path.join(flags['output'],"00_settings")) )

    print("  --->Flags written\n")

def convert_protID_to_common_names(flag_values, hits):
    '''
    Converts BLAST query protein ID's to common gene names using
    .csv if -c flag was invoked.

    Parameters
    ----------
    flag_values: dict 
        Arg values capatured from keyboard in housekeeping.py
 
    Returns
    -------
    '''
    if flag_values['common_name']:
        common_name_df = pd.read_csv(flag_values['common_name'], names=['protID','common'])

        # Create a subset dataframe of every row with duplicated common names
        subset_df=common_name_df.duplicated(subset='common', keep=False)

        # Grab index/row numbers of duplicated common names
        duplicated_rows=common_name_df[subset_df].index.tolist()
        
        # Duplicated common names are reassigned unique names (common name + protein id)
        for index in duplicated_rows:
            common_name_df.iloc[index,1] = \
            str(common_name_df.iloc[index,1]) + ' ' + str(common_name_df.iloc[index,0])
        
        # Convert dataframe with prot ID + unique common names to a dict for further use
        common_name = \
        pd.Series(common_name_df.common.values,index=common_name_df.protID).to_dict()

        '''
        We need to loop through and update key names in the 'hits' dictionary, but
        we can't loop through keys that are being updated. So, here we create 
        'unmutable' lists of the dictionary keys to loop through
        '''
        common_name_list=list(common_name.keys())
        hits_list=list(hits.keys())
        
        # Rename keys that are protein IDs with common gene names. Using lists from comment above.
        for k in common_name_list:
            for hit in hits_list:
                if k == hit:
                    print("Common name:'{}' is in hit name: '{}'".format(k,hit))
                    hits[common_name[k]] = hits[hit]
                    del hits[hit]