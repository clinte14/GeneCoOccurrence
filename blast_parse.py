import re, os, sys
from Bio.Blast import NCBIXML
import pandas as pd

def parse_merge_BLAST(flag_values):
    output_hits = {}
    file = []
#    all_files = os.listdir(flag_values['output'])
#    all_files = flag_values['output']
    input_dir = flag_values['input']
    print("Input Dir: " + input_dir)
    
    # take list of files from input directory ('all_files') and place only files ending with ".xml" in 'file' variable
    for index, value in enumerate(os.listdir(input_dir)):
        if value[-4:] == ".xml":
            file.append(value)
    
    # throw error and quit if no .xml BLAST input file(s)
    if len(file) == 0:
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("No '.xml' formatted output detected in /01_BLAST_results/; quitting  now...")
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        sys.exit(0)
    
    file.sort()
    print(file)
    # loop through .xml files in /01_BLAST_results
    for i in file:
        print("Parsing BLAST file {}...".format(i))
        print()
        xml_handle = open(os.path.join(input_dir,i))
        xml_record = NCBIXML.read(xml_handle)
        file_name = i[:-4]
        output_hits[file_name] = []
        
        # loop through single .xml file
        # writes these to 'hits' dictionary with gene as key and value is list containing associated hits
        for hit_title in xml_record.alignments:
            blast_hits=hit_title.title.split('>')

            for h in blast_hits:
                if h.startswith('ref') == False:
                    db=h.split('|')[0]
                    species=re.findall('\[(.*?)\]',h)[0]
                    protein_ID=re.findall('\|(.*?)\|',h)[0].split('.')[0][:-4]
                    UID=species+ '[' + db + '|' + protein_ID + 'xxxx]'
                    output_hits[file_name].append(UID)
    print('  --->Finished Parsing')
    print()
    return output_hits 

# write hits dictionary to Pandas DF, query gene is column and hits are rows. Write dataframe to 02_PA_matrix      
# folder as PA_matrix.csv.
def create_pa(hits, flag_values):
    print('Merging BLAST XML files...')
    # create 'merged_BLAST.csv' file in /02_PA_Matrix/ containing all values in 'hits' dictionary from parse_merge_BLAST()
    merged_df = pd.DataFrame({ key:pd.Series(value) for key, value in hits.items() })
    merged_df = merged_df.reindex(sorted(merged_df.columns), axis=1)
    merged_df.to_csv(flag_values['output'] + '/02_PA_matrix/' + 'merged_BLAST.csv')
    
    print("  --->Wrote 'merged_BLAST.csv' to '{}".format(flag_values['output']) + "/02_PA_matrix/'..." + '\n') 
    
#    # read merged_df back from file, in case user made manual changes
#    merged_csv = pd.read_csv(flag_values['output'] + '/02_PA_matrix/' + 'merged_BLAST.csv', index_col=0)

    # create 'pa_df', an empty binary presence/absence Pandas DF. make sure rows and columns are neatly sorted for us humans    
    pa_df_columns = merged_df.columns.values.astype(str).tolist()
    pa_df_columns.sort()
    
    pa_df_rows = pd.unique(merged_df[pa_df_columns].values.ravel()).astype(str).tolist()
    pa_df_rows.sort()
    
    pa_df = pd.DataFrame(index=pa_df_rows, columns=pa_df_columns)

    
    # populate values into presence/absence pa_df from merged_df
    print("Creating presence/absence matrix...") 
    for column_num, column_name in enumerate(merged_df):
#        total_rows = merged_df[column].count()
        for row, value in enumerate(merged_df[column_name]):
            pa_df.at[value,column_name] = 1      
    # Replace nan with zero
    pa_df = pa_df.fillna(value=0)
    pa_df.to_csv(flag_values['output'] + '/02_PA_matrix/' + 'pa_matrix.csv')
    print("  --->Wrote 'pa_matrix.csv' to '{}".format(flag_values['output']) + "/02_PA_matrix/'..." + '\n') 