import os.path
import argparse, datetime, json
import pandas as pd
def create_folders(dir_path):
    dirs = ['00_log',
            '01_blast_input', 
            '02_PA_matrix', 
            '03_correlation_calcs', 
            '04_correlog_values', 
            '05_final_outputs']

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
    
    parser = argparse.ArgumentParser(description="Correlogy", prog='correlog')

    parser.add_argument("-i", "--input", 
        required=True, 
        help="File with BLAST results containing homologues in .xml format. REQUIRED.")

    parser.add_argument("-o", "--output", 
        default=dir_path,
        help="Desired directory of output file. Default=Current Directory.")

    parser.add_argument("-c", "--common_name",
        help="A .csv file that allows conversion of BLAST query protein ID's \
        to common gene names.")

    args = parser.parse_args()
    arg_dict=vars(args)

    # Convert relative paths (eg './' in BASH) to absolute paths
    arg_dict['input'] = os.path.abspath(arg_dict['input'])
    arg_dict['output'] = os.path.abspath(arg_dict['output'])
    arg_dict['common_name'] = os.path.abspath(arg_dict['common_name'])

    print("dir_path: " + str(dir_path))
    return arg_dict


# save flags options and paramaters to "command.txt" file in 00_Settings folder
def save_flags(flags):
    with open(os.path.join(flags['output'],"00_log",'command.txt'), 'w') as file:

        currentDT = datetime.datetime.now()

        file.write("Executed the following in " + flags['output'] + " at " \
            + str(currentDT) + ":\n" + json.dumps(flags))

    print("Writing flags to 'command.txt' in '{}'... ".\
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
                if k in hit:
                    print("Common name:'{}' is in hit name: '{}'".format(k,hit))
                    hits[common_name[k]] = hits[hit]
                    del hits[hit]