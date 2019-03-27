# Correlogy v0.1a
Gene Co-Occurrence Across Phylogeny

### How to run
##### Choose BASH or IDE (eg Spyder)
1. **BASH**: Must set 'use_IDE' variable to 'False' in main.py
	- minimal example: 
`python3 main.py -i input_gene_list`
	- full example: 
`some text`
	- usage: 
`correlog [-h] -i INPUT [-o OUTPUT] [--evalue EVALUE] [--entrez ENTREZ] [-s SKIPBLAST]`
	- get help, list all flag defaults: 
`python3 main.py -h`

2. **IDE**: Must set 'use_IDE' variable to 'True' in main.py
	- Hardcode flag settings in `housekeeping.py --> parse_input() --> arg_dict dictionary`

### Known bugs

1. NCBI BLASTing can sometimes time out, resulting in an empty XML file in 01_BLAST_results. Workaround is to use '-s' flag to skip BLAST query and put your own BLAST results in XML format (one file per query) into 01_BLAST-results.

2. blast_parse.py --> create_pa() does 'flattening' of pa_df via 'pa_df_rows = pd.unique(merged_df[pa_df_columns].values.ravel()).tolist()'. This returns some hits as "NaN", not sure why probably related to flattening.

### To do 
in `housekeeping.py --> create_folders()`
  1) validate input file exists (otherwise throw error)
  2) validate output directory exists (otherwise throw error)
