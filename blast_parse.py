import re, os, sys
from Bio.Blast import NCBIXML
import pandas as pd

def clean_BLAST_xml(flag_values):
    '''
    Search and remove 'CREATE_VIEW from BLAST .xml input file.

    Bug workaround for NCBIWWW. Occasionally NCBI BLAST online writes 'CREATE_VIEW' into 
    the BLAST .xml output with no associated '<>' tags, and this causes xml parsing in 
    Biopython to fail. See: https://github.com/biopython/biopython/issues/3342
    '''

    print('Cleaning up BLAST results: {}'.format(flag_values['input']))
    with open(flag_values['input'], 'r') as file:
        xml_data = file.read()
        xml_data = xml_data.replace('CREATE_VIEW','')

    with open(flag_values['input'], 'w') as file:
        file.write(xml_data)

def parse_merge_BLAST(flag_values):
    '''
    Converts BLAST query protein ID's to common gene names using \
    .csv if -c flag was invoked.

    Parameters
    ----------
    flag_values: dict 
        Arg values capatured from keyboard in housekeeping.py
 
    Returns
    -------
    '''
    output_hits = {}
    file = []

    xml_input = flag_values['input']
    print("Parsing BLAST file '{}'...".format(xml_input))
    xml_handle = open(xml_input)
    xml_records = NCBIXML.parse(xml_handle)

    for xml_record in xml_records:
        gene_name=xml_record.query_id
        output_hits[gene_name]=[]

        for hit_title in xml_record.alignments:
            blast_hits=hit_title.title.split('>')

            for h in blast_hits:
                if h.startswith('ref') == False:

                    if bool(h.split('|')):
                        db=h.split('|')[0]
                    else:
                        db='NODB'
                        UID += (db + '|')

                    if bool(re.findall('\[(.*?)\]',h)):
                        species=re.findall('\[(.*?)\]',h)[0]
                    else:
                        species='NOSPECIES'

                    if bool(re.findall('\|(.*?)\|',h)):
                        protein_ID=re.findall('\|(.*?)\|',h)[0].split('.')[0][:-4]
                    else:
                        protein_ID='NOPROTEINID'

                    UID=species+ '[' + db + '|' + protein_ID + 'xxxx]'
                    #print("UID: '{}".format(UID))
                    output_hits[gene_name].append(UID)
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
    merged_df_dir=os.path.join(flag_values['output'],'02_PA_matrix','merged_BLAST.csv')
    merged_df.to_csv(merged_df_dir)
    
    print("  --->Wrote 'merged_BLAST.csv' to '{}'".format(merged_df_dir + '\n'))
    
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
    pa_dir = os.path.join(flag_values['output'] ,'02_PA_matrix','pa_matrix.csv')
    pa_df.to_csv(pa_dir)
    print("  --->Wrote 'pa_matrix.csv' to '{}'".format(pa_dir + '\n'))