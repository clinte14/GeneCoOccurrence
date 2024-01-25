import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os
import numpy as np
import graphviz
import networkx as nx
from pyvis.network import Network

#https://towardsdatascience.com/visualising-stocks-correlations-with-networkx-88f2ee25362e
#https://www.cl.cam.ac.uk/teaching/1314/L109/tutorial.pdf
#https://valdanchev.github.io/reproducible-data-science-python/network_analysis.html
#see above! Meso-scale network diagnostics...NETWORK COMMUNITIES
#https://vamshij.com/blog/network-science/basics-network-science-1/
#https://stats.stackexchange.com/questions/138325/clustering-a-correlation-matrix
#https://stackoverflow.com/questions/32303217/how-to-use-pearson-correlation-as-distance-metric-in-scikit-learn-agglomerative
#https://stats.stackexchange.com/questions/275720/does-any-other-clustering-algorithms-take-correlation-as-distance-metric-apart
#https://online.stat.psu.edu/stat555/node/85/

def create_network_map(flag_values, network_list):
    #breakpoint()
    for count, value in enumerate(network_list):
        network_list[count][2] = '{:.3f}'.format(value[2])

    #DG = nx.DiGraph()
    #DG.add_weighted_edges_from(network_list)
    nx_graph = nx.DiGraph()
    for item in network_list:
        nx_graph.add_node(item[0])
    for item in network_list:
        #add weighted line
        if float(item[2]) > 0:
            nx_graph.add_edge(item[0],item[1],label=item[2],weight=float(item[2])*10)
        else:
            nx_graph.add_edge(item[0],item[1],label=item[2],weight=1)

    nt = Network('800px', '1080px')
    nt.from_nx(nx_graph)
    nt.toggle_physics(True)
    nt.show_buttons(True)
    os.chdir(os.path.join(flag_values['output'], '03_visual_output'))
    nt.write_html("MSR_network.html",notebook=False)

    # Create 'dot' engine style graphviz network visual
    network = graphviz.Digraph('Maximum_Related_Networks_dot', engine='dot', comment='')

    for i in network_list:
        co_score = i[2]
        if float(i[2]) <= 0:
            edge_thickness = 1      
        else:
            edge_thickness = int(i[2].split('.')[-1][0])
        network.edge(i[0], i[1], penwidth=str(edge_thickness), xlabel=co_score, dir='none')
    network.render(directory=os.path.join(flag_values['output'], '03_visual_output'))

    # Create 'circo' engine style graphviz network visual
    network = graphviz.Digraph('Maximum_Related_Networks_circo', engine='circo', comment='') 

    for i in network_list:
        co_score = i[2]
        if float(i[2]) <= 0:
            edge_thickness = 1      
        else:
            edge_thickness = int(i[2].split('.')[-1][0])
        network.edge(i[0], i[1], penwidth=str(edge_thickness), xlabel=co_score, dir='none')

    network.render(directory=os.path.join(flag_values['output'], '03_visual_output'))


    #need print statement here to console

def create_heatmap(Wij_df, flag_values):

    # Seaborn likes all data as floats
    Wij_df = Wij_df.astype(float)

    # Generate a mask for the upper triangle; True = do NOT show
    mask = np.zeros_like(Wij_df, dtype=bool)
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
        square=False,   # Force cells to be square
        linewidths=.5, # Width of lines that divide cells
        cbar_kws={'label' : 'Co-Occurrence Value','shrink': .75}  # Extra kwargs for the legend; 
        # in this case, shrink
    )

    res.set_xticklabels(res.get_xmajorticklabels(), fontsize = 12)
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize = 12)
    res.yaxis.label.set_size(20)

    heatmap_dir = os.path.join(flag_values['output'], '03_visual_output','heatmap.png')
    fig.tight_layout()
    fig.savefig(heatmap_dir)
    print('Heatmap written to: {}'.format(heatmap_dir) + '\n')