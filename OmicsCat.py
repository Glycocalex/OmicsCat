# Part 1 of my similarity calculation script.
# Import modules

import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.spatial.distance import cdist
import ipywidgets as widgets
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
# pip install C:\Users\Alex Bennett\Desktop\rust\OmicsCat\distance_calc\target\wheels\distance_calc-0.1.0-cp311-none-win_amd64.whl


# Data import and preparation
# Ensure data is structured so that each row represents a species and each column represents a time point

def oc_import(z, label):
    omics = pd.read_csv(z)
    label = str(label)
    omics.insert(1, "type", label)
    omics.iloc[:, 2:] = omics.iloc[:, 2:].astype(np.float64)
    x = omics
    return x

tomics = oc_import("C:\\Users\\Alex Bennett\\Desktop\\Python\\Omics integration\\Omics integration\\Distance calculations\\data\\renal testing\\renal tomics small.csv", 't')
pomics = oc_import("C:\\Users\\Alex Bennett\\Desktop\\Python\\Omics integration\\Omics integration\\Distance calculations\\data\\renal testing\\renal pomics small.csv", 'p')
gomics = oc_import("C:\\Users\\Alex Bennett\\Desktop\\Python\\Omics integration\\Omics integration\\Distance calculations\\data\\renal testing\\renal gomics small.csv", 't')

def oc_norm(z):  # Data normalisation
    # Currently just a z - score normalisation
    # Expects a pd.DF with molecules by row, timepoints/sample by column and numerical expression values in cells
    # The first columns should contain the molecule name, and second should flag what type of molceules (protein, transcript, ect.)

    x = z.copy()
    x.iloc[:, 2:] = stats.zscore(x.iloc[:, 2:], axis=1)
    x = pd.DataFrame(x)
    return x

tomics = oc_norm(tomics)
pomics = oc_norm(pomics)
gomics = oc_norm(gomics)

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
momics = oc_cat(tomics, pomics)

def oc_single(z):

    # Use to calculate distances for a single timepoint/ replicate.
    # Not  sure if this is at all useful, but I suppose it produces an output of molecular elements with
    # a similar magnitude of output.

    for i in range(len(z)):
        for j in range(i + 1, len(z)):
            outputdists = []
            species1 = z[i]
            species2 = z[j]
            distance = np.linalg.norm(species1 - species2, axis=0)
            outputdists.append(([i], [j], distance))
            return outputdists

# Todo: can we deal with replicate data here? Necessary for statistical test
    # Is the best way to deal with replicates is to assume no sig. differences between induviduals at each timepoint?
    # With this assumption, we can treat replicates as 'matched' and compare each replicate only to the corresponding
    # replicates at each timepoint
# Todo: For stats test, read up on permutation test and mantel test


def oc_dist_full(z, tps, reps):  # Calculate distances (z = normalised expression data, tps = timepoints, reps = replicates)
    # Calculate distances with a timeseries for-loop:
    # Todo: fix for new import and labeling
    # Code below is attempt at producing a full distance matrix for DTW that can handle varying tps/reps
    momicsdist = np.zeros((z.shape[0], tps))
    z_np = z.to_numpy()
    #z_col_labels = z_np[0, :]
    z_row_labels = z_np[:, 0]
    z_vals = z_np[:, 1:]

    z_vals.astype(np.float64)
    z_vals = stats.zscore(z_vals, axis=1)

    for i in range(0, tps):
        momicsdist[:, i] = np.mean(z_vals[:, (i*reps):((i*reps)+reps)], axis=1)

    distances = cdist(momicsdist, momicsdist, metric='euclidean')

    return distances


def oc_dist_simple(z, threshold):
    # The code below produces a pair-wise distance matrix for each sample (samples defined by column)
    momicsdist = []
    momicszs = np.array(z.iloc[:, 2:])
    for column in range(momicszs.shape[1]):
        column_values = momicszs[:, column:column + 1]
        distances = cdist(column_values, column_values, metric='euclidean')
        for i in range(len(momicszs)):
            for j in range(i + 1, len(momicszs)):
                species1 = tuple(z.iloc[i, 0:2])
                species2 = tuple(z.iloc[j, 0:2])
                distance = distances[i, j]  # i and j represent the indexes of what molecules are being compared
                momicsdist.append(([column], species1, species2, i, j, distance)) # Todo: are i/j necessary?

    #  Producing a time-averaged distance matrix:
    columns, species1, species2, i, j, distances = zip(*momicsdist)
    i_values = np.array(list(i), dtype=float).flatten()
    j_values = np.array(list(j), dtype=float).flatten()
    columnsf = np.concatenate(columns)
    agdist = pd.DataFrame({'Column': columnsf, 'i': i_values, 'j': j_values, 'Species1': species1, 'Species2': species2, 'Distance': distances})
    agdist['Column'] = pd.Series(agdist['Column'])
    avdist = agdist.groupby(['i', 'j', 'Species1', 'Species2'])['Distance'].mean().reset_index()

    #  Threshold the time-averaged distance matrix:
    y = float(threshold)
    thresholddf = avdist['Distance'].quantile(y)
    coregindex = pd.DataFrame(avdist['Distance'] <= thresholddf)
    coregdists = avdist[coregindex['Distance']]
    coregdists.reset_index(drop=True, inplace=True)

    return coregdists

dists = oc_dist_simple(momics, 0.05)

def oc_warp(x, threshold):
    """Performs dynamic time warping -
    Answers two questions:  1) How similar are two signals?
                            2) Which timepoints correspond to one another (from two signals with different periods,
                                durations or resolutions)
    x = a distance matrix
    threshold = a cutoff determining what alignment costs are reported (higher cost indicates a weaker relationship)"""

# Note: DTW requires a full distance matrix of each waveform. Currently, oc_dist only produces distances between
# each species at corresponding time points.

    cost_map = []

    for i in range(0, len(x[:, ])):
        for j in range(i+1, len(x[:, ])):
            costmatrix = np.zeros((len(x[i, :])+1, len(x[j, :])+1))
            costmatrix[0, :] = np.inf
            costmatrix[:, 0] = np.inf
            v1 = x[i, :]
            v2 = x[j, :]
            cost = np.linalg.norm(v1 - v2)
            costmatrix[i + 1, j + 1] = cost + min(costmatrix[i, j + 1], costmatrix[i + 1, j], costmatrix[i, j])
            costmatrix = costmatrix[1:, 1:]
            i = costmatrix.shape[0] - 1
            j = costmatrix.shape[1] - 1
            matches = []
            mappings_v1 = [list() for v in range(costmatrix.shape[0])]
            mappings_v2 = [list() for v in range(costmatrix.shape[1])]
            while i > 0 or j > 0:
                matches.append((i, j))
                mappings_v1[i].append(j)
                mappings_v2[j].append(i)
                option_diag = costmatrix[i - 1, j - 1] if i > 0 and j > 0 else np.inf
                option_up = costmatrix[i - 1, j] if i > 0 else np.inf
                option_left = costmatrix[i, j - 1] if j > 0 else np.inf
                move = np.argmin([option_diag, option_up, option_left])
                if move == 0:
                    i -= 1
                    j -= 1
                elif move == 1:
                    i -= 1
                else:
                    j -= 1
            matches.append((0, 0))
            mappings_v1[0].append(0)
            mappings_v2[0].append(0)
            matches.reverse()
            for mp in mappings_v1:
                mp.reverse()
            for mp in mappings_v2:
                mp.reverse()
                cost_map.extend([matches, costmatrix[-1, -1], mappings_v1, mappings_v2])

    return cost_map

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
