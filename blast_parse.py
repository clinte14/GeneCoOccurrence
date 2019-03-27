import time, re
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd
    
# performs online NCBIWWW.qblast if 'skipBLAST' flag is set to False (default). Output is XML formatted BLAST 
# results (one per search query/line in input file) in 01_BLAST_results folder.  
def query_BLAST(flag_values):
    file = open(flag_values['input'], 'r').read().split('\n') #rU, was r
    for i in file:
        print("Currently BLASTing {}...".format(i))
        result_handle = NCBIWWW.qblast("blastp", "refseq_protein", i, expect=flag_values['evalue'], entrez_query = flag_values['entrez'])
        time.sleep(3)
    
    out_handle = open(flag_values['output'] + "/" + "01_BLAST_results" + "/" + i + "_BLAST_results.xml", "w")
    print(i)
    out_handle.write(result_handle.read())
    result_handle.close()  

# parse XML files located in 01_BLAST_results folder. Discard MULTISPECIES hits. Strip remaining BLAST hits to  
# only include specices name [hit_id] and add to hits dictionary, query gene is key and associated value is list   
# of strings formatted as species name [hit_id]
def parse_BLAST(flag_values):
    hits={}
    file = open(flag_values['input'], 'r').read().split('\n')
    for i in file:
        print("-----Parsing XML BLAST output {}...".format(i))
        xml_handle = open (flag_values['output'] + "/" + "01_BLAST_results" + "/" + i + "_BLAST_results.xml", "r")
        xml_record = NCBIXML.read(xml_handle)
#        hits[i] = ["test", "some", "strings"]
        hits[i] = []
        
        for hit_title in xml_record.alignments:
#            print("hit:",hit_title.title)
            if "MULTISPECIES" not in hit_title.title:
                stripped = re.findall('\[(.*?)\]', hit_title.title)
#                stripped.replace('[', '')
#                stripped2 = re.findall(r"\w+", stripped)
                #\[(.*?)\]
#                print(stripped[0] + ' ' + hit_title.hit_id)
                hits[i].append(stripped[0] + ' ' + hit_title.hit_id)
        print('->Finished Parsing')
    print()
    return hits

# write hits dictionary to Pandas DF, query gene is column and hits are rows. Write dataframe to 02_PA_matrix      
# folder as PA_matrix.csv.
def create_pa(hits, flag_values):
    print('-----Merging BLAST XML files...')
    merged_df = pd.DataFrame({ key:pd.Series(value) for key, value in hits.items() })
    merged_df.to_csv(flag_values['output'] + '/02_PA_matrix/' + 'merged_BLAST.csv')
    print("->Wrote 'merged_BLAST.csv' to '{}".format(flag_values['output']) + "/02_PA_matrix/'..." + '\n') 
    
#    # read merged_df back from file, in case user made manual changes
#    merged_csv = pd.read_csv(flag_values['output'] + '/02_PA_matrix/' + 'merged_BLAST.csv', index_col=0)

    # create binary presence/absence Pandas DF    
    pa_df_columns = merged_df.columns.values.tolist()
    pa_df_rows = pd.unique(merged_df[pa_df_columns].values.ravel()).tolist()
    pa_df = pd.DataFrame(index=pa_df_rows, columns=pa_df_columns)
    
    # copy values from  from merged_df to binary presence/absence pa_df
    print("-----Creating presence/absence matrix...") 
    for column_num, column_name in enumerate(merged_df):
#        total_rows = merged_df[column].count()
        for row, value in enumerate(merged_df[column_name]):
            pa_df.at[value,column_name] = 1      
        # Replace nan with zero
    pa_df = pa_df.fillna(value=0)
    pa_df.to_csv(flag_values['output'] + '/02_PA_matrix/' + 'pa_matrix.csv')
    print("->Wrote 'pa_matrix.csv' to '{}".format(flag_values['output']) + "/02_PA_matrix/'..." + '\n') 