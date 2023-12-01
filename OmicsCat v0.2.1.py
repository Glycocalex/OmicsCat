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

tomics = pd.read_csv("C:\\Users\\Alex Bennett\\Desktop\\Python\\Omics integration\\Omics integration\\Distance calculations\\data\\renal testing\\renal gomics small.csv")
tomics_np = tomics.to_numpy()
tomics_col_labels = tomics_np[0, :]
tomics_row_labels = tomics_np[:, 0]
tomics_np = tomics_np[:, 1:]
tomics_np = tomics_np.astype(np.float64)

pomics = pd.read_csv("C:\\Users\\Alex Bennett\\Desktop\\Python\\Omics integration\\Omics integration\\Distance calculations\\data\\renal testing\\renal pomics small.csv")
pomics_np = pomics.to_numpy()
pomics_col_labels = pomics_np[0, :]
pomics_row_labels = pomics_np[:, 0]
pomics_np = pomics_np[:, 1:]
pomics_np = pomics_np.astype(np.float64)

gomics = pd.read_csv("C:\\Users\\Alex Bennett\\Desktop\\Python\\Omics integration\\Omics integration\\Distance calculations\\data\\renal testing\\renal tomics small.csv")
gomics_np = gomics.to_numpy()
gomics_col_labels = gomics_np[0, :]
gomics_row_labels = gomics_np[:, 0]
gomics_np = gomics_np[:, 1:]
gomics_np = gomics_np.astype(np.float64)


def oc_import(z):
    omics = pd.read_csv(z)
    omics = omics.to_numpy()
    omics_c_labs = omics[0, :]  # Omics column labels
    omics_r_labs = omics[:, 0]  # Omics row labels
    omics = omics[1:, 1:]
    omics = pd.DataFrame(omics, columns=omics_c_labs, index=omics_r_labs)
    omics.astype(np.float64)
    x = omics
    return x


def oc_norm(z):  # Data normalisation
    # Currently just a z - score normalisation
    # Expects a pd.DF with molecules by row, timepoints/sample by column and numerical expression values in cells

    z = stats.zscore(z, axis=1)
    z = pd.DataFrame(z)
    return z


def oc_cat(*args):  # Data aggregation
    # Not necessary to normalise different datasets first (but you really should)
    # momics = multi-omics. Not creative, but fun to say!

    momics = pd.concat([*args], axis=0)
    momics = momics.to_numpy()
    return momics


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

# Todo: work on functions past this point, address labs in oc_dist
# Todo: can we deal with replicate data here? Necessary for statistical tests
    # Is the best way to deal with replicates is to assume no sig. differences between induviduals at each timepoint?
    # With this assumption, we can treat replicates as 'matched' and compare each replicate only to the corresponding
    # replicates at each timepoint
# Todo: For stats test, read up on permutation test and mantel test


def oc_dist_full(z, tps, reps):  # Calculate distances (z = normalised expression data, tps = timepoints, reps = replicates)
    # Calculate distances with a timeseries for-loop:

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


def oc_dist_simple(z):
    # The code below produces a pair-wise distance matrix for each sample (samples defined by column)
    momicsdist = []
    momicszs = np.array(z)
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
    agdist = pd.DataFrame({'Column': columnsf, 'i': i_valuesf, 'j': j_valuesf, 'Distance': distances})
    agdist['Column'] = pd.Series(agdist['Column'])
    avdist = agdist.groupby(['i', 'j'])['Distance'].mean().reset_index()

    #  Threshold the time-averaged distance matrix:
    threshold = avdist['Distance'].quantile(0.01)
    coregindex = pd.DataFrame(avdist['Distance'] <= threshold)
    coregdists = avdist[coregindex['Distance']]
    coregdists.reset_index(drop=True, inplace=True)

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
    coregdistslab = pd.merge(momics_row_labels_df, coregdists, left_on='i_left', right_on='i', how='left')
    coregdistslab.rename(columns={'Species': 'Species_Labels_i'}, inplace=True)
    coregdistslab.drop(columns=['i_left', 'j_left'], inplace=True)
    # Absolute hatchet job on the row below, no idea why NaNs were being added upon merging.
    coregdistslab = coregdistslab.dropna()
    coregdistslab = pd.merge(coregdistslab, momics_row_labels_df, left_on='j', right_on='j_left', how='left')
    coregdistslab.rename(columns={'Species': 'Species_Labels_j'}, inplace=True)
    coregdistslab.drop(columns=['i_left', 'j_left'], inplace=True)

    return coregdistslab

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





def oc_graph(x):
    #  Produce network graph:
    g = nx.Graph()

    for _, row in x.iterrows():
        gene1 = row['Species_Labels_i']
        gene2 = row['Species_Labels_j']
        distance = row['Distance']
        g.add_edge(gene1, gene2, weight=distance)
    pos = nx.spring_layout(g)
    edge_labels = {(gene1, gene2): g.edges[gene1, gene2]['weight'] for gene1, gene2 in g.edges()}
    nx.draw_networkx_nodes(g, pos, node_size=300, node_color='skyblue')
    nx.draw_networkx_labels(g, pos, font_size=10, font_family='sans-serif')
    nx.draw_networkx_edges(g, pos)

    # Following two lines are for simple visualisations, recommend to use pyvis in the block below
    # plt.axis('off')
    # plt.show()

    #  Attempting to visualise with pyvis
    net = Network()
    net.from_nx(g)
    net.toggle_physics(True)  # Change to True for small graphs, provides dynamic/more fun graph interactions
    # net.show("file\\path.html")  # Note: a html file will be created which contains the network - open in your web browser
    net.show("C:\\Users\\Alex Bennett\\Desktop\\Python\\Omics integration\\OmicsCat\\OmicsCatDev\\outputs\\OCTest4.html")

