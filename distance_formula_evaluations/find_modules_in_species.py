#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
  ============================================================================
  
File Description:
    Obtain gene modules associated with Kegg species ID
    
Input: 
    .csv file with column for species ID's (header = 'k_species_id') 
    and species name ('k_species_name')
    
Output: 
    .csv of pandas data frame:
            k_id species_name module
    gene1
    gene2
    
Author: Clinton Elg 
Contact: elg6752@vandals.uidaho.edu or clintelg@gmail.com

  ===========================================================================
"""
import numpy as np
import pandas as pd
import os.path, sys
from Bio.KEGG import REST
from urllib.error import URLError, HTTPError

SPECIES_LIST_FILENAME = "/home/clint/Documents/correlogy/v2/species_list_varied_taxa.csv"
GENE_MODULES_OUTPUT_FILENAME = "/home/clint/Documents/correlogy/v2/genes_modules_pathogens_df.csv"
NCBI_GENES_FROM_MODULES = "/home/clint/Documents/correlogy/v2/genes_modules_aaseq_pathogens.csv"

#Obtain gene modules associated with Kegg species ID
def find_modules_in_species(df):
    #module_dict contains {Kegg bacteria ID:list[Kegg module ID's]}
    modules_in_species_dict = {}
    k_id_REST_failures = []
    all_modules = []
    counter = 0
    
    for k_id in df.k_id:
        current_module_list = []
        counter+=1
        print("{}/{}...Retrieving modules in '{}'".format(counter, len(df.k_id), k_id))

        
        try:
            bacteria_modules = REST.kegg_("module", k_id).read()
        except (URLError, HTTPError, TimeoutError):
            print("WARNING: URLError, HTTPError, or TimeoutError with species: '{}'".format(k_id))
            k_id_REST_failures.append(k_id)
        else:
            #kegg returns modules as string, so parse line by line and split
            for line in bacteria_modules.rstrip().split('\n'):
                module, description = line.split('\t')
                #module = module.strip('md:') MIGHT NEED THIS LATER?
                #add modules stripped from kegg to current_module_list[]
                current_module_list.append(module)
                all_modules.append(module)
                
            #write current_module_list to modules_in_species_dict
            #module_dict contains {Kegg bacteria ID:list[Kegg module ID's]}    
            modules_in_species_dict[k_id] = current_module_list
            
#    new_df = pd.DataFrame.from_dict(module_dict, orient ="index")
#    new_df = new_df.transpose()di
    if len(k_id_REST_failures)==0:
        print("{}/{} All species had associated modules successfully returned".format((counter-1), len(df.k_id)))
    else:
        print("WARNING: {}/{} Species had URLError, HTTPError, or Timeout Error".format(len(k_id_REST_failures),len(df.k_id)))
        print("The following failed: {}".format(k_id_REST_failures))
    return all_modules, modules_in_species_dict

#
def find_genes_in_modules(species_df):
    counter = 0
    REST_failures = []
    temp_gene_list = []
    temp_module_list = []
    temp_k_id = []
    temp_k_species_name = []

    
    for k_id in species_df.k_species_id:
        df = pd.DataFrame(columns=['k_gene_id', 'module', 'k_species_id','k_species_name', 'aaseq'])
        counter+=1
        print("{}/{}...Retrieving modules and associated genes in species '{}'".format(counter, len(species_df.k_species_id), k_id))
        
        try:
            kegg_handle = REST.kegg_link("module", k_id).read()
        except (URLError, HTTPError, TimeoutError):
            print("***WARNING***: URLError, HTTPError, or TimeoutError with species: '{}'".format(k_id))
            print("\n")

            REST_failures.append(k_id)
        else:
            for line in kegg_handle.rstrip().split('\n'):
                gene, module = line.split('\t')
                print(gene, module, end=',')
                temp_gene_list.append(gene)
                temp_module_list.append(module)
                temp_k_id.append(k_id)
                temp_k_species_name.append(species_df.loc[species_df['k_species_id']==k_id, 'k_species_name'].iloc[0])
            print("\n")
            df['k_gene_id'] = temp_gene_list
            df['module'] = temp_module_list
            df['k_species_id'] = temp_k_id
            df['k_species_name'] = temp_k_species_name
    if len(REST_failures)==0:
        print("{}/{} All species had associated modules successfully returned".format((counter), len(species_df.k_species_id)))
    else:
        print("{}/{} Species had associated modules successfully returned".format(((counter)-len(REST_failures)), len(species_df.k_species_id)))
        print("***WARNING***: {}/{} Species had URLError, HTTPError, or Timeout Error".format(len(REST_failures),len(species_df.k_species_id)))
        print("The following failed: {}".format(REST_failures))       

            
#    df.to_csv('find_genes_in_modules_df.csv', encoding='utf-8')
#    print("DEBUG")
    df.to_csv(GENE_MODULES_OUTPUT_FILENAME,index=False, header=True)
    return df
            
            
            

"""    
    for k_id in species_df.k_species_id:
        df = pd.DataFrame(columns=['k_gene_id', 'module'])
        counter+=1
        print("{}/{}...Retrieving modules and associated genes in species '{}'".format(counter, len(k_id), k_id), end=',')
        try:
            k_gene_module = REST.kegg_link("module", k_id).read()
        except (URLError, HTTPError, TimeoutError):
            print("WARNING: URLError, HTTPError, or TimeoutError with species: '{}'".format(k_id))
            REST_failures.append(module)
        else:
            for line in k_gene_module.rstrip().split('\n'):
                gene, module = line.split('\t')
                temp_gene_list.append(gene)
                temp_module_list.append(module)
        df['k_gene_id'] = temp_gene_list
        df['module'] = temp_module_list
        df['k_species_id'] = k_id
        df['k_species_name'] = species_df.loc[k_id,'k_species_name']
        df.append(df)
    
    if len(module_REST_failures)==0:
        print("{}/{} All modules had associated genes successfully returned".format((counter-1), len(all_modules)))
    else:
        print("WARNING: {}/{} Species had URLError, HTTPError, or Timeout Error".format(len(module_REST_failures),len(all_modules)))
        print("The following failed: {}".format(module_REST_failures))
        
    return genes_in_module_dict
        
    print("DEBUG")
"""
def convert_k_gene_id_to_aaseq(df,species_df):
    counter = 0
    REST_failures = []
    no_ncbi_gene_id_returned = []
    
    for k_gene_id in df.k_gene_id:
        counter+=1
        
        #print("{}/{}...Converting KEGG gene ID '{}'to NCBI gene ID".format(counter, len(df.k_gene_id), k_gene_id))
        print("{}/{}...Retreiving aaseq associated with KEGG gene ID '{}'".format(counter, len(df.k_gene_id), k_gene_id))

        
        try:
            my_handle = REST.kegg_get(k_gene_id, "aaseq").read()
            df.at[(counter-1),'aaseq'] = my_handle
            #kegg_handle = REST.kegg_conv("ncbi-geneid", k_gene_id).read()
            #print(kegg_handle)
        except (URLError, HTTPError, TimeoutError):
            print("***WARNING***: URLError, HTTPError, or TimeoutError with KEGG gene ID: '{}'".format(k_gene_id))
            print("\n")

            REST_failures.append(k_gene_id)
#        else:
#            if kegg_handle !="\n":
#                ncbi_gene_id = kegg_handle.split(':')
#                df.at[(counter-1),'ncbi_gene_id'] = ncbi_gene_id[2]
#            else:
#                print("KEGG returned no NCBI ID for '{}', obtaining ntseq from KEGG again".format(k_gene_id))
#                df.at[(counter-1),'ntseq'] = REST.kegg_get(k_gene_id, "ntseq").read()
#                no_ncbi_gene_id_returned.append(k_gene_id)


#    if len(REST_failures)==0:
#        print("{}/{}...All KEGG gene ID's converted to NCBI gene ID's".format(counter, len(df.k_gene_id)))
#    else:
#        print("{}/{} Species had associated modules successfully returned".format(((counter)-len(REST_failures)), len(df.k_gene_id)))
    print("*Note*:{}/{} Kegg ID's had URLError, HTTPError, or Timeout Error".format(len(REST_failures),len(df.k_gene_id)))
    print("The following failed: {}".format(no_ncbi_gene_id_returned))
#        

    return df
    
def main():

#    gi_output_df = pd.Series()
    species_df = pd.read_csv(SPECIES_LIST_FILENAME, header=0, sep=",")
#    all_modules, modules_in_species_dict = find_modules_in_species(species_df)
        
    if os.path.isfile(GENE_MODULES_OUTPUT_FILENAME):
        print("***KEGG MODULES FOUND: {}***".format(GENE_MODULES_OUTPUT_FILENAME))
        df = pd.read_csv(GENE_MODULES_OUTPUT_FILENAME, dtype=object)

        
    else:
        print("***QUERING KEGG FOR MODULES***")
        df = find_genes_in_modules(species_df)


    #df = find_genes_in_modules(species_df)
    #df = pd.read_csv('genes_modules_df.csv')
    print("***QUERING KEGG TO OBTAIN NCBI NUMBER OF KEGG GENE ID***")
    df = convert_k_gene_id_to_aaseq(df, species_df)
    df.to_csv(NCBI_GENES_FROM_MODULES,index=False, header=True)
    print("***DEBUG: End of program***")
    
#    print(df)



if __name__ == "__main__":
    main()