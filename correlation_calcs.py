import numpy as np
import pandas as pd
import math as math
import csv
import graphviz
import os


def correlation_calcs(flag_values):
    '''
    '''
    # Open prescence/abscence matrix from /02_PA_matrix/ directory 
    with open(os.path.join(flag_values['output'],'02_PA_matrix','pa_matrix.csv'), 'r') as infile:
            # Create dataframe from BLAST results, header = row# w/ column names
            pa_df = pd.read_csv(infile, sep=',', header=0)
    
    #clean up indexing and names in our pa_df (prescence/abscence dataframe)
    pa_df = pa_df.rename(columns={"Unnamed: 0" : "species"})
    
    #an index of gene_names. Do not include indexer row '0' or 'species' column '1' of pa_df (prescence/abscence dataframe)
    gene_names =  pa_df.dtypes.index[1:]
    
    #this function calculates all the terms needed for pearson correlation
    #FUNCTION: calc_pearson_terms(pa_df, gene_names), returns E_i (dict gene:# species containing gene) N is total # of species
    #Cij_df contains number of species containing both genes i and j
    E_i, N, Cij_df = calc_pearson_terms(flag_values, pa_df, gene_names)

    #this function calculates pearson correlation values
    #FUNCTION: pearson_corr_calc(flag_values, pa_df, gene_names, E_i, N, Cij_df), returns Rij_inverse_df (inverse matrix of pearson correlations)
    Rij_inverse_df = pearson_corr_calc(flag_values, pa_df, gene_names, E_i, N, Cij_df)
    
    # Wij = Partial correlation
    Wij_df = pd.DataFrame(index=gene_names, columns=gene_names)
    
    max_rel_df = pd.DataFrame(columns=['source', 'target', 'Wij_value', 'correlogy_group', 'correlogy_dir', 'directed'])
    
    # Create empty maximum relationship dataframe 
    max_rel_df.source = gene_names
    for i, name in enumerate(gene_names):
        max_rel_df.loc[i, 'directed']='TRUE'
    

    # Calculate all Wij permutations.
    for column_index, column_name in enumerate(gene_names):
        # Loop through all columns to the right of current column (column_index), so column_index + 1 = column to right
        for column_to_right_index, column_to_right_name in enumerate(gene_names[column_index+1:]):
            #print(column_name)
            Pij = Rij_inverse_df.loc[column_name, column_to_right_name]
            Pii = Rij_inverse_df.loc[column_name, column_name]
            Pjj = Rij_inverse_df.loc[column_to_right_name, column_to_right_name]
            #print("Pij=", Pij, "Pii=", Pii, "Pjj=", Pjj)
            Wij_df.at[column_name, column_to_right_name] = -((Pij) / (math.sqrt(Pii * Pjj)))
    Wij_df.to_csv(os.path.join(flag_values['output'],'03_correlation_calcs','Wij_df_triangle.csv'))
    print("->Wrote 'Wij_df_triangle.csv' to '{}".format(os.path.join(flag_values['output']),"04_correlog_values..") + '\n') 

           
    # Convert Wij_df from a triangle matrix into a symmetrical matrix
    for column in gene_names:
        for row in gene_names:
            # Create main diagonal of 1
            if column == row:
                Wij_df.loc[column, row]=1
            # Transpose rows/columns with columns/rows to fill in bottom triangles
            else:
                Wij_df.loc[row, column] = Wij_df.loc[column, row]

    Wij_df.to_csv(os.path.join(flag_values['output'],'03_correlation_calcs','Wij_df.csv'))

    #need print statement here to console
    return Wij_df

def calc_mrs(Wij_df, flag_values):
    Wij_outerDict = {}
    gene_names =  Wij_df.dtypes.index
    
    f = open(os.path.join(flag_values['output'],'04_correlog_values',"Wij_matrix.csv"), "w")
    writer = csv.writer(f)
    
    network_list=[]
    for row in gene_names:
        #this portion of loop writes Wij_matrix.csv file
        #sort columns by a given row highest to lowest
        Wij_df = Wij_df.sort_values(by=row, ascending=False, axis=1, na_position='last')

        current_columns = Wij_df.columns.tolist()
        current_columns.remove(row)
        current_columns.insert(0,row)
        writer.writerow(current_columns)
        
        current_row = Wij_df.loc[row,:].tolist()
        current_row.remove(1)
        current_row.insert(0,"")
        writer.writerow(current_row)
        Wij_outerDict[row] = dict(zip(Wij_df.columns, Wij_df.loc[row,:]))

        #this data is for creating visual graphics for maximum related networks
        temp_list=[]
        temp_list.append(current_columns[0])
        temp_list.append(current_columns[1])
        temp_list.append(current_row[1])
        network_list.append(temp_list)
    f.close()
    return network_list

########################################################################################################################
#this function calculates all the terms needed for pearson correlation
#FUNCTION: calc_pearson_terms(pa_df, gene_names), returns E_i (dict gene:# species containing gene) N is total # of species
#Cij_df contains number of species containing both genes i and j
def calc_pearson_terms(flag_values, pa_df, gene_names):
 

    # E_i is dictionary(gene:number species containing gene)
    E_i = {}
    for i in gene_names: #do not include indexer row '0' or 'species' column '1'
        E_i.update({i:len(pa_df.loc[pa_df[i] == 1, ['species']])})
        
    N = len(pa_df.species)
    
    # Cij = Number of species containing both i and j. 
    Cij_df = pd.DataFrame(index=gene_names, columns=gene_names)
    # Calculates all Cij permutations.Move through columns L to R, stopping at each column to compare 
    # current columns gene presence in species (rows) to the following columns gene presence.
        # Loop through all columns except last
    for column_index, column_name in enumerate(gene_names[0:-1]):
        print(column_name)
        # Loop through all columns to the right of current column (column_index), so column_index + 1 = column to right
        for column_to_right_index, column_to_right_name in enumerate(gene_names[column_index+1:]):
            # Cij = # species containing both genes i and j, zeroes out everytime you compare new columns
            Cij = 0
            # Loop through current column and column to right (provided by above 'for' loop) row by row 
            # to check if they both have '1' values ie gene is present in both organisms
            for current_index, current_item in enumerate(pa_df.loc[:, column_name]):
                #print(item1) #testing purposes
                if current_item == 1 and pa_df.loc[current_index, column_to_right_name] == 1:
                        Cij += 1
            print("  --->", column_to_right_name, Cij)
            Cij_df.at[column_name, column_to_right_name] = Cij
    # write Cij_df to disk
    cij_dir=os.path.join(flag_values['output'],'03_correlation_calcs','Cij_df.csv')
    Cij_df.to_csv(cij_dir)
    print("Wrote 'Cij_df.csv' aka # species with both gene 'i' & gene 'j' to '{}".format(cij_dir + '\n')) 
    # Debug print(Cij, index, df.loc[index, "species"]) 
    
    return E_i, N, Cij_df

########################################################################################################################
#this function calculates pearson correlation values
#FUNCTION: pearson_corr_calc(flag_values, pa_df, gene_names, E_i, N, Cij_df), returns Cij_df (matrix of pearson correlations)
def pearson_corr_calc(flag_values, pa_df, gene_names, E_i, N, Cij_df):
    #Rij = Pearson Correlation Values of gene i to j. stored in Rij_df
    Rij_df = pd.DataFrame(index=gene_names, columns=gene_names, dtype = float)

    # Calculate all Rij permutations. Move through Cij _df columns L to R, stopping at each column to 
    # calculate current column in relationship to following columns.
        # Loop through all columns except last
    for column_index, column_name in enumerate(gene_names[0:-1]):
        print(column_name)
        # Loop through all columns to the right of current column (column_index), so column_index + 1 = column to right
        for column_to_right_index, column_to_right_name in enumerate(gene_names[column_index+1:]):
            current_Cij = Cij_df.loc[column_name, column_to_right_name]
            current_E_i = E_i[column_name]
            current_E_j = E_i[column_to_right_name]
            Rij_df.at[column_name, column_to_right_name] = ((current_Cij * N) - (current_E_i * current_E_j)) / (
                    math.sqrt(current_E_i) * math.sqrt(current_E_j) * math.sqrt(N-current_E_i) * math.sqrt(N-current_E_j))
    # write Rij_df to disk
    rij_triangle_dir=os.path.join(flag_values['output'],'03_correlation_calcs','Rij_df_triangle.csv')
    Rij_df.to_csv(rij_triangle_dir)
    print("->Wrote 'Rij_df_triangle.csv' aka Pearson Correlations to '{}".format(rij_triangle_dir + '\n'))
            
    # Convert Rij_df from a triangle matrix into a symmetrical matrix. We need this to create our inverse matrix later.
    for column in gene_names:
        for row in gene_names:
            # Create main diagonal of 0's
            if column == row:
                Rij_df.loc[column, row]=1
            # Transpose rows/columns with columns/rows to fill in bottom triangles
            else:
                Rij_df.loc[row, column] = Rij_df.loc[column, row]
    # write Rij_df symmetrical to disk
    rij_sym_dir=os.path.join(flag_values['output'], '03_correlation_calcs', 'Rij_df_symmetrical.csv')
    Rij_df.to_csv(rij_sym_dir)
    print("->Wrote 'Rij_df_symmetrical.csv' to '{}".format(rij_sym_dir + '\n'))
    
    # Create inverse matrix from symmetrical Rij matrix
    Rij_inverse_df = pd.DataFrame(np.linalg.pinv(Rij_df.values), Rij_df.columns, Rij_df.index)

    # write Rij_inverse_df to disk
    Rij_inverse_dir=os.path.join(flag_values['output'],'03_correlation_calcs','Rij_inverse_df.csv')
    Rij_inverse_df.to_csv(Rij_inverse_dir)
    print("->Wrote 'Rij_inverse_df.csv' to '{}".format(Rij_inverse_dir + '\n'))
    
    # Verify matrix multiplication results in the identity matrix.
#    identify_matrix_test = Rij_inverse_df.dot(Rij_df)
    return Rij_inverse_df
