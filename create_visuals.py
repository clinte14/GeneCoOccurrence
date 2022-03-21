import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os
import numpy as np
import graphviz

#https://towardsdatascience.com/visualising-stocks-correlations-with-networkx-88f2ee25362e
#https://www.cl.cam.ac.uk/teaching/1314/L109/tutorial.pdf
#https://valdanchev.github.io/reproducible-data-science-python/network_analysis.html
#see above! Meso-scale network diagnostics...NETWORK COMMUNITIES
#https://vamshij.com/blog/network-science/basics-network-science-1/
#https://stats.stackexchange.com/questions/138325/clustering-a-correlation-matrix
#https://stackoverflow.com/questions/32303217/how-to-use-pearson-correlation-as-distance-metric-in-scikit-learn-agglomerative
#https://stats.stackexchange.com/questions/275720/does-any-other-clustering-algorithms-take-correlation-as-distance-metric-apart
#https://online.stat.psu.edu/stat555/node/85/
def create_network_map(Wij_df, flag_values, network_list):
    # Create 'dot' engine style graphviz network visual
    network = graphviz.Digraph('Maximum_Related_Networks_dot', engine='dot', comment='')

    for i in network_list:
        edge_thickness = str(i[2]*5)
        correlog_score = '{:.4f}'.format(i[2])
        network.edge(i[0], i[1], penwidth=edge_thickness, xlabel=correlog_score)

    network.render(directory=os.path.join(flag_values['output'], '05_final_outputs'))

    # Create 'circl' engine style graphviz network visual
    network = graphviz.Digraph('Maximum_Related_Networks_circo', engine='circo', comment='') 

    for i in network_list:
        edge_thickness = str(i[2]*5)
        correlog_score = '{:.4f}'.format(i[2])
        network.edge(i[0], i[1], penwidth=edge_thickness, xlabel=correlog_score)

    network.render(directory=os.path.join(flag_values['output'], '05_final_outputs'))

    #need print statement here to console

def create_heatmap(Wij_df, flag_values):

    # Seaborn likes all data as floats
    Wij_df = Wij_df.astype(float)

    # Generate a mask for the upper triangle; True = do NOT show
    mask = np.zeros_like(Wij_df, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    # Set up the matplotlib figure
    fig = plt.figure(num=None, figsize=(
        len(Wij_df.index)+2, 
        len(Wij_df.columns)),
        dpi=300
        )
#    f, ax = plt.subplots(figsize=(11, 9))
#    f, ax = plt.subplots(figsize=(len(Wij_df.index)+2, len(Wij_df.columns)))
#    breakpoint()
    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    # More details at https://seaborn.pydata.org/generated/seaborn.heatmap.html

    res = sns.heatmap(
        Wij_df,          # The data to plot
        mask=mask,     # Mask some cells
        cmap=cmap,     # What colors to plot the heatmap as
        annot=True,    # Should the values be plotted in the cells?
        vmax=1,       # The maximum value of the legend. All higher vals will be same color
        vmin=-1,      # The minimum value of the legend. All lower vals will be same color
        center=0,      # The center value of the legend. With divergent cmap, where white is
        square=True,   # Force cells to be square
        linewidths=.5, # Width of lines that divide cells
        cbar_kws={'label' : 'Correlog Value','shrink': .75}  # Extra kwargs for the legend; 
        # in this case, shrink by 50%
    )

    res.set_xticklabels(res.get_xmajorticklabels(), fontsize = 12)
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize = 12)
    res.yaxis.label.set_size(20)

    heatmap_dir = os.path.join(flag_values['output'], '05_final_outputs','heatmap.png')
    fig.savefig(heatmap_dir)
    print('Heatmap written to: {}'.format(heatmap_dir) + '\n')