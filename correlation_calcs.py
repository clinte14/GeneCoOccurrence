import numpy as np
#import matplotlib.pyplot as  plt
import pandas as pd
import math as math

def correlation_calcs(flag_values):
    with open(flag_values['output'] + '/02_PA_matrix/' + 'pa_matrix.csv', 'r') as infile:
            # Create dataframe from BLAST results, header = row# w/ column names
            pa_df = pd.read_csv(infile, sep=',', header=0)
            
    pa_df = pa_df.rename(columns={"Unnamed: 0" : "species"})
    
    # Dict w/ gene as key with corresponding # of species that gene is found in as. Term 'Ei' in Pearson correlation.
    gene_species_count = {}
    
    # Populate 'gene_species_count' ie 'Ei'. See previous comment.
    gene_column = pa_df.dtypes.index[1:] #do not include indexer row '0' or 'species' column '1'
    for i in gene_column:
        counter = 0
        for x in pa_df.loc[:, i]:
            if x == 1:
                counter += 1
            gene_species_count.update({i:counter})
            
    N = len(pa_df.species)

    # Cij = Number of species containing both i and j. Used to calculate Rij (see next comment)
    Cij_df = pd.DataFrame(index=gene_column, columns=gene_column)
    
    # Rij = Pearson Correlation Values of gene i to j. Used to calculate an inverse matrix of itself, which is
    # then used to calculate Wij(see next comment)
    Rij_df = pd.DataFrame(index=gene_column, columns=gene_column, dtype = float)
    
    Rij_inverse_df = pd.DataFrame(index=gene_column, columns=gene_column, dtype = float)

    # Wij = Partial correlation
    Wij_df = pd.DataFrame(index=gene_column, columns=gene_column)
    
    max_rel_df = pd.DataFrame(columns=['source', 'target', 'Wij_value', 'correlogy_group', 'correlogy_dir', 'directed'])
    
    # Create empty maximum relationship dataframe 
    max_rel_df.source = gene_column
    for i, name in enumerate(gene_column):
        max_rel_df.loc[i, 'directed']='TRUE'
    
    # Calculates all Cij permutations.Move through columns L to R, stopping at each column to compare 
    # current columns gene presence in species (rows) to the following columns gene presence.
        # Loop through all columns except last
    for column_index, column_name in enumerate(gene_column[0:-1]):
        print(column_name)
        # Loop through all columns to the right of current column (column_index), so column_index + 1 = column to right
        for column_to_right_index, column_to_right_name in enumerate(gene_column[column_index+1:]):
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
    Cij_df.to_csv(flag_values['output'] + '/03_Correlation_calcs/' + 'Cij_df.csv')
    print("Wrote 'Cij_df.csv' aka # species with gene 'i' & gene 'j' to '{}".format(flag_values['output']) + "/03_Correlation_calcs/..." + '\n') 
    # Debug print(Cij, index, df.loc[index, "species"]) 
    
    # Calculate all Rij permutations. Move through Cij _df columns L to R, stopping at each column to 
    # calculate current column in relationship to following columns.
        # Loop through all columns except last
    for column_index, column_name in enumerate(gene_column[0:-1]):
        print(column_name)
        # Loop through all columns to the right of current column (column_index), so column_index + 1 = column to right
        for column_to_right_index, column_to_right_name in enumerate(gene_column[column_index+1:]):
            current_Cij = Cij_df.loc[column_name, column_to_right_name]
            Ei = gene_species_count[column_name]
            Ej = gene_species_count[column_to_right_name]
            Rij_df.at[column_name, column_to_right_name] = ((current_Cij * N) - (Ei * Ej)) / (
                    math.sqrt(Ei) * math.sqrt(Ej) * math.sqrt(N-Ei) * math.sqrt(N-Ej))
    # write Rij_df to disk
    Rij_df.to_csv(flag_values['output'] + '/03_Correlation_calcs/' + 'Rij_df_triangle.csv')
    print("->Wrote 'Rij_df_triangle.csv' aka Pearson Correlations to '{}".format(flag_values['output']) + "/03_Correlation_calcs/..." + '\n')
            
    # Convert Rij_df from a triangle matrix into a symmetrical matrix. We need this to create our inverse matrix later.
    for column in gene_column:
        for row in gene_column:
            # Create main diagonal of 0's
            if column == row:
                Rij_df.loc[column, row]=1
            # Transpose rows/columns with columns/rows to fill in bottom triangles
            else:
                Rij_df.loc[row, column] = Rij_df.loc[column, row]
    # write Rij_df symmetrical to disk
    Rij_df.to_csv(flag_values['output'] + '/03_Correlation_calcs/' + 'Rij_df_symmetrical.csv')
    print("->Wrote 'Rij_df_symmetrical.csv' to '{}".format(flag_values['output']) + "/03_Correlation_calcs/..." + '\n')
    
    # Create inverse matrix from symmetrical Rij matrix
    Rij_inverse_df = pd.DataFrame(np.linalg.pinv(Rij_df.values), Rij_df.columns, Rij_df.index)
    # write Rij_inverse_df to disk
    Rij_inverse_df.to_csv(flag_values['output'] + '/03_Correlation_calcs/' + 'Rij_inverse_df.csv')
    print("->Wrote 'Rij_inverse_df.csv' to '{}".format(flag_values['output']) + "/03_Correlation_calcs/..." + '\n') 
    
    # Verify matrix multiplication results in the identity matrix.
#    identify_matrix_test = Rij_inverse_df.dot(Rij_df)
    
    # Calculate all Wij permutations.
    for column_index, column_name in enumerate(gene_column):
        # Loop through all columns to the right of current column (column_index), so column_index + 1 = column to right
        for column_to_right_index, column_to_right_name in enumerate(gene_column[column_index+1:]):
            #print(column_name)
            Pij = Rij_inverse_df.loc[column_name, column_to_right_name]
            Pii = Rij_inverse_df.loc[column_name, column_name]
            Pjj = Rij_inverse_df.loc[column_to_right_name, column_to_right_name]
            #print("Pij=", Pij, "Pii=", Pii, "Pjj=", Pjj)
            Wij_df.at[column_name, column_to_right_name] = -((Pij) / (math.sqrt(Pii * Pjj)))
    Wij_df.to_csv(flag_values['output'] + '/03_Correlation_calcs/' + 'Wij_df_triangle.csv')
    print("->Wrote 'Wij_df_triangle.csv' to '{}".format(flag_values['output']) + "/04_Correlog_values/..." + '\n') 

           
    # Convert Wij_df from a triangle matrix into a symmetrical matrix
    for column in gene_column:
        for row in gene_column:
            # Create main diagonal of 0's
            if column == row:
                Wij_df.loc[column, row]=0
            # Transpose rows/columns with columns/rows to fill in bottom triangles
            else:
                Wij_df.loc[row, column] = Wij_df.loc[column, row]
    Wij_df.to_csv(flag_values['output'] + '/03_Correlation_calcs/' + 'Wij_df.csv')
    print("->Wrote 'Wij_df.csv' to '{}".format(flag_values['output']) + "/04_Correlog_values/..." + '\n') 

"""

    max_rel = 0
    max_rel_dict = {}
    max_rel_df = pd.DataFrame(columns=['gene_i', 'gene_j','Wij_max'])
    test = gene_column.to_series()
    another = max_rel_df['gene_i'].append(test)
    okay = another = max_rel_df['gene_i'].append(test, ignore_index = True)
    
        
    for i, index in enumerate(gene_column):
        print(index,"...")
        max_rel = Wij_df.at[index,'VC0175'] #VC0175 in VSP-I, VC0490
        print("Max_rel starts at:",max_rel)
        for column in gene_column:
           #print(column) 
           if max_rel < Wij_df.at[index, column]:
                max_rel = Wij_df.at[index, column]
                print("-->New max_rel =", max_rel, "@", column)
                max_rel_df.at[i, 'target'] = column
                max_rel_df.at[i, 'Wij_value'] = max_rel
                Wij_dict[index] = [max_rel, column]
                
    max_rel = 0
    for loop_index, row in enumerate(gene_column):
        print(row,"...")
        max_rel = Wij_df.at[index,'VC0175'] #VC0175 in VSP-I, VC0490
        print("Max_rel starts at:",max_rel)
        for column in gene_column:
           #print(column) 
           if max_rel < Wij_df.at[row, column]:
                max_rel = Wij_df.at[row, column]
                print("-->New max_rel =", max_rel, "@", column)
                max_rel_df.at[loop_index, 'target'] = column
                max_rel_df.at[loop_index, 'Wij_value'] = max_rel
        min

#                Wij_dict[index] = [max_rel, column]
    
    
    
    #for i in Wij_dict:
    #   Wij_dict[i]
"""