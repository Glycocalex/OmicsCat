import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.spatial.distance import cdist
import ipywidgets as widgets
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
pip install r'C:\Users\Alex Bennett\Desktop\rust\OmicsCat\distance_calc\target\wheels\distance_calc-0.1.0-cp311-none-win_amd64.whl'
# Todo: Setup rust lib for first time users


def oc_import(z, label): # Ensure data is structured so that each row represents a species and each column represents a time point
    omics = pd.read_csv(z)
    label = str(label)
    omics.insert(1, "type", label)
    omics.iloc[:, 2:] = omics.iloc[:, 2:].astype(np.float64)
    x = omics
    return x


def oc_norm(z):  # Data normalisation
    # Currently just a z - score normalisation
    # Expects a pd.DF with molecules by row, timepoints/sample by column and numerical expression values in cells
    # The first columns should contain the molecule name, and second should flag what type of molceules (protein, transcript, ect.)

    x = z.copy()
    x.iloc[:, 2:] = stats.zscore(x.iloc[:, 2:], axis=1)
    x = pd.DataFrame(x)
    return x

def oc_cat(*args):  # Data aggregation Todo: introduce a tag for each dataframe indicating contents (genes, proteins, ect.)
    # Not necessary to normalise different datasets first (but you really should)
    # momics = multi-omics. Not creative, but fun to say!

    for df in args:
        df.columns.values[0] = 'Species name'
        #level = str(df)
        #df[level] = level

    molecules = pd.concat([df.iloc[:, 0:2] for df in args])
    values = pd.concat([df.iloc[:, 2:] for df in args])
    momics = pd.concat([molecules, values], axis=1)
    return momics

# Todo: can we deal with replicate data here? Necessary for statistical test
    # Is the best way to deal with replicates is to assume no sig. differences between induviduals at each timepoint?
    # With this assumption, we can treat replicates as 'matched' and compare each replicate only to the corresponding
    # replicates at each timepoint
# Todo: For stats test, read up on permutation test and mantel test


def oc_graph(x): # Todo: fix for new labelling and add ego graph
    #  Produce network graph:
    g = nx.Graph()

    # Create network: nodes= molecules, edges lengths = distances
    node_colours = {'t': 'skyblue', 'p': 'red'}
    nodes = set()

    for _, row in x.iterrows():
        species1, level1 = map(str, row['Species1'])
        species2, level2 = map(str, row['Species2'])
        distance = row['Distance']

        if species1 not in nodes:
            g.add_node(species1, level=level1)
            nodes.add(species1)
            #print(f"Node {species1} has level {level1}")
            #nx.draw_networkx_nodes(g, pos={species1:(0,0)}, nodelist=[species1], node_size=300,
            #                       node_color=node_colours.get(level1))

        if species2 not in nodes:
            g.add_node(species2, level=level2)
            nodes.add(species2)
            #print(f"Node {species2} has level {level2}")
            #nx.draw_networkx_nodes(g, pos={species2: (0, 0)}, nodelist=[species2], node_size=300,
            #                       node_color=node_colours.get(level2))

        g.add_edge(species1, species2, weight=distance)

    #pos = nx.spring_layout(g)

   # # Colour nodes based upon level, by iterating through nodes and setting colour
   # for node, data in g.nodes(data=True):
   #     node_type = data.get('level', None)
   #     colour = node_colours.get(node_type)
   #     nx.draw_networkx_nodes(g, pos, nodelist=[node], node_size=300, node_color=colour)
#
   # edge_labels = {(species1, species2): g.edges[species1, species2]['weight'] for species1, species2 in g.edges()}
   # #nx.draw_networkx_nodes(g, pos, node_size=300)
   # nx.draw_networkx_labels(g, pos, font_size=10, font_family='sans-serif')
   # nx.draw_networkx_edges(g, pos)
#
    #  Visualise with pyvis #Todo: Visualisation not finished. By constructing the graph in pyvis we can colour node by level but it looks shit. Is there a way to colour node by level in netowrkx?
    net = Network()

    for node, data in g.nodes(data=True):
        node_type = data.get('level', None)
        colour = node_colours.get(node_type, 'gray')
        net.add_node(node, color=colour, borderWidth=2, size=15)

    for edge in g.edges():
        net.add_edge(edge[0], edge[1], value=((g.edges[edge]['weight'])/4), width=1,  color='gray')
        net.set_edge_smooth('dynamic')
    #net.add_edges(g.edges())
    #for edge_id in net.edges:
    #    net.edges[edge_id]['color'] = 'gray'
    #net.set_edge_smooth('dynamic')
#

    #net.from_nx(g)
    net.toggle_physics(True)  # Change to True for small graphs, provides dynamic/more fun graph interactions
    # net.show("file\\path.html")  # Note: a html file will be created which contains the network - open in your web browser
    #net.show_buttons(filter_=['physics'])
    net.show("C:\\Users\\Alex Bennett\\Desktop\\Python\\Omics integration\\OmicsCat\\OmicsCatDev\\outputs\\OCTestegograph.html", notebook=False)



tomics = oc_import("C:\\Users\\Alex Bennett\\Desktop\\Python\\Omics integration\\Omics integration\\Distance calculations\\data\\renal testing\\renal tomics small.csv", 't')
pomics = oc_import("C:\\Users\\Alex Bennett\\Desktop\\Python\\Omics integration\\Omics integration\\Distance calculations\\data\\renal testing\\renal pomics small.csv", 'p')
gomics = oc_import("C:\\Users\\Alex Bennett\\Desktop\\Python\\Omics integration\\Omics integration\\Distance calculations\\data\\renal testing\\renal gomics small.csv", 't')

tomics = oc_norm(tomics)
pomics = oc_norm(pomics)
gomics = oc_norm(gomics)
momics = oc_cat(tomics, pomics)
dists = calculate_distance_matrix(momics, 0.05)
oc_graph(dists)


## Todo: make an egograph function


#def oc_egograpgh(node, radius, distance):
#
g = nx.Graph()

start_node = "Leo1"
radius = 2

ego_graph = nx.ego_graph(g, start_node, radius=radius)

pos = nx.spring_layout(ego_graph)
edge_labels = {(gene1, gene2): ego_graph.edges[gene1, gene2]['weight'] for gene1, gene2 in ego_graph.edges()}
nx.draw_networkx_nodes(ego_graph, pos, node_size=300, node_color='skyblue')
nx.draw_networkx_labels(ego_graph, pos, font_size=10, font_family='sans-serif')
nx.draw_networkx_edges(ego_graph, pos)


# Following two lines are for simple visualisations, recommend to use pyvis in the block below
# plt.axis('off')
# plt.show()

#  Attempting to visualise with pyvis
net = Network()
net.from_nx(ego_graph)
net.toggle_physics(True)  # Change to True for small graphs, provides dynamic/more fun graph interactions
net.show("C:\\Users\\Alex Bennett\\Desktop\\Python\\Omics integration\\OmicsCat\\OmicsCatDev\\outputs\\OCTestEgograph.html", notebook=False)
