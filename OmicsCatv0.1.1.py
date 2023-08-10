# Part 1 of my similarity calculation script.
# Import modules

import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.spatial.distance import cdist
from sklearn.metrics.pairwise import euclidean_distances
import ipywidgets as widgets
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network


# Data import and preparation
# Ensure data is structured so that each row represents a species and each column represents a time point

tomics = pd.read_csv("file\\path.csv")
tomics_np = tomics.to_numpy()
tomics_col_labels = tomics_np[0, :]
tomics_row_labels = tomics_np[:, 0]
tomics_np = tomics_np[:, 1:]
tomics_np = tomics_np.astype(np.float64)

pomics = pd.read_csv("file\\path.csv")
pomics_np = pomics.to_numpy()
pomics_col_labels = pomics_np[0, :]
pomics_row_labels = pomics_np[:, 0]
pomics_np = pomics_np[:, 1:]
pomics_np = pomics_np.astype(np.float64)

gomics = pd.read_csv("file\\path.csv")
gomics_np = gomics.to_numpy()
gomics_col_labels = gomics_np[0, :]
gomics_row_labels = gomics_np[:, 0]
gomics_np = gomics_np[:, 1:]
gomics_np = gomics_np.astype(np.float64)


# Data normalisation (z-score)

tomicszs = stats.zscore(tomics_np, axis=1)
pomicszs = stats.zscore(pomics_np, axis=1)
gomicszs = stats.zscore(gomics_np, axis=1)

# Data aggregation
tomicszsdf = pd.DataFrame(tomicszs)
pomicszsdf = pd.DataFrame(pomicszs)
gomicszsdf = pd.DataFrame(gomicszs)
momicszsdf = pd.concat([tomicszsdf, pomicszsdf, gomicszsdf], axis=0)
momicszs = momicszsdf.to_numpy()


## Use code below to calculate distances for a single timepoint
#tomicsdist = []
#for i in range(len(tomicszs)):
#    for j in range(i + 1, len(tomicszs)):
#        species1 = tomicszs[i]
#        species2 = tomicszs[j]
#        distance = np.linalg.norm(species1 - species2, axis=0)
#        tomicsdist.append(([i], [j], distance))


# Calculate distances with a timeseries for-loop:
# The code below produces a pair-wise distance matrix at each time point (t defined by column)
momicsdist = []
for column in range(momicszs.shape[1]):
    column_values = momicszs[:, column:column + 1]
    distances = cdist(column_values, column_values, metric='euclidean')
    for i in range(len(momicszs)):
        for j in range(i + 1, len(momicszs)):
            species1 = column_values[i]
            species2 = column_values[j]
            distance = distances[i, j]
            momicsdist.append(([column], [i], [j], distance))

#  Producing a time-averaged distance matrix:
columns, i_values, j_values, distances = zip(*momicsdist)
i_valuesl = list(i_values)
i_valuesf = np.array(i_valuesl, dtype=float).flatten()
j_valuesl = list(j_values)
j_valuesf = np.array(j_valuesl, dtype=float).flatten()
columnsf = np.concatenate(columns)
AgDist = pd.DataFrame({'Column': columnsf, 'i': i_valuesf, 'j': j_valuesf, 'Distance': distances})
AgDist['Column'] = pd.Series(AgDist['Column'])
AvDist = AgDist.groupby(['i', 'j'])['Distance'].mean().reset_index()

#  Threshold the time-averaged distance matrix:
threshold = AvDist['Distance'].quantile(0.01)
CoRegIndex = pd.DataFrame(AvDist['Distance'] <= threshold)
CoRegDists = AvDist[CoRegIndex['Distance']]
CoRegDists.reset_index(drop=True, inplace=True)

#  Labelling
tomics_row_labels_df = pd.DataFrame(tomics_row_labels)
tomics_row_labels_df.columns = ['Species']
pomics_row_labels_df = pd.DataFrame(pomics_row_labels)
pomics_row_labels_df.columns = ['Species']
gomics_row_labels_df = pd.DataFrame(gomics_row_labels)
gomics_row_labels_df.columns = ['Species']
momics_row_labels_df = pd.concat([tomics_row_labels_df, pomics_row_labels_df, gomics_row_labels_df], axis=0)

momics_row_labels_df['i_left'] = momics_row_labels_df.reset_index(drop=False).index.astype(float)
momics_row_labels_df['j_left'] = momics_row_labels_df.reset_index(drop=False).index.astype(float)
CoRegDistsLab = pd.merge(momics_row_labels_df, CoRegDists, left_on='i_left', right_on='i', how='left')
CoRegDistsLab.rename(columns={'Species': 'Species_Labels_i'}, inplace=True)
CoRegDistsLab.drop(columns=['i_left', 'j_left'], inplace=True)
# Absolute hatchet job on the row below, no idea why NaNs were being added upon merging.
CoRegDistsLab = CoRegDistsLab.dropna()
CoRegDistsLab = pd.merge(CoRegDistsLab, momics_row_labels_df, left_on='j', right_on='j_left', how='left')
CoRegDistsLab.rename(columns={'Species': 'Species_Labels_j'}, inplace=True)
CoRegDistsLab.drop(columns=['i_left', 'j_left'], inplace=True)


#  Produce network graph:
G = nx.Graph()
for _, row in CoRegDistsLab.iterrows():
    gene1 = row['Species_Labels_i']
    gene2 = row['Species_Labels_j']
    Distance = row['Distance']
    G.add_edge(gene1, gene2, weight=Distance)
pos = nx.spring_layout(G)
edge_labels = {(gene1, gene2): G.edges[gene1, gene2]['weight'] for gene1, gene2 in G.edges()}
nx.draw_networkx_nodes(G, pos, node_size=300, node_color='skyblue')
nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif')
nx.draw_networkx_edges(G, pos)
#  nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)

# Following two lines are for simple visualisations, recommend to use pyvis in the block below
# plt.axis('off')
# plt.show()

#  Attempting to visualise with pyvis
net = Network()
net.from_nx(G)
net.toggle_physics(True)  # Change to True for small graphs, provides dynamic/more fun graph interactions
net.show("file\\path.html")  # Note: a html file will be created which contains the network - open in your web browser
