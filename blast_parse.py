import time, re, os, sys
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd
    
# performs online NCBIWWW.qblast if 'skipblast' flag is set to False (default). Output is XML formatted BLAST 
# results (one per search query/line in input file) in 01_BLAST_results folder.  
def query_BLAST(flag_values):
    file = open(flag_values['input'], 'r').read().split('\n') #rU, was r
    print("BLASTing using evalue: '{}' and entrez query '{}'".format(flag_values['evalue'], flag_values['entrez']))
    for i in file:
        print("Currently BLASTing {}...".format(i))
        result_handle = NCBIWWW.qblast("blastp", "refseq_protein", i, expect=flag_values['evalue'], entrez_query = flag_values['entrez'])
        time.sleep(3)
        
        out_handle = open(flag_values['output'] + "/" + "01_BLAST_results" + "/" + i + "_BLAST_results.xml", "w")
        print('  --->Finished {}'.format(i))
        out_handle.write(result_handle.read())
    result_handle.close()  

# parse XML files located in 01_BLAST_results folder. Discard MULTISPECIES hits. Strip remaining BLAST hits to  
# only include specices name [hit_id] and add to hits dictionary, query gene is key and associated value is list   
# of strings formatted as species name [hit_id]
def parse_merge_BLAST(flag_values):
    hits = {}
    file = []
    all_files = os.listdir(flag_values['output'] + "/" + "01_BLAST_results" + "/")
    
    # take list of files from input directory ('all_files') and place only files ending with ".xml" in 'file' variable
    for index, value in enumerate(all_files):
        if all_files[index][-4:] == ".xml":
            file.append(all_files[index])
    
    # throw error and quit if no .xml BLAST input file(s)
    if len(file) == 0:
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("No '.xml' formatted output detected in /01_BLAST_results/; quitting  now...")
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        sys.exit(0)
    
    file.sort()
    
    # loop through .xml files in /01_BLAST_results
    for i in file:
        print("Parsing XML BLAST output {}...".format(i))
        xml_handle = open(flag_values['output'] + "/" + "01_BLAST_results" + "/" + i, "r")
        xml_record = NCBIXML.read(xml_handle)
# debug        hits[i] = ["test", "some", "strings"]
        hits[i] = []
        
        # loop through single .xml file, remove 'MULTISPECIES' hits. Such hits are always non-species taxonomic levels (eg genus, family, order)
        # parses remaing hits, strips off unwanted data
        # writes these to 'hits' dictionary with gene as key and value is list containing associated hits
        for hit_title in xml_record.alignments:
# debug            print("hit:",hit_title.title)
            if "MULTISPECIES" not in hit_title.title:
                stripped = re.findall('\[(.*?)\]', hit_title.title)
                hits[i].append(stripped[0])
#08:08                hits[i].append(stripped[0] + ' ' + hit_title.hit_id)
        print('  --->Finished Parsing')
    print()
    return hits

# write hits dictionary to Pandas DF, query gene is column and hits are rows. Write dataframe to 02_PA_matrix      
# folder as PA_matrix.csv.
def create_pa(hits, flag_values):
    print('Merging BLAST XML files...')
    # create 'merged_BLAST.csv' file in /02_PA_Matrix/ containing all values in 'hits' dictionary from parse_merge_BLAST()
    merged_df = pd.DataFrame({ key:pd.Series(value) for key, value in hits.items() })
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