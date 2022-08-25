# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and extrats a thresholded network.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import ensurePathExists
import argparse


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Defaults to spermatocyte")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    databases = ["GeneMANIA", "IntAct", "CPDB", "Mentha", "StringDB"]
    parser.add_argument("--database", default="StringDB", type=str, choices = databases, help="Database. Defaults to StringDB")
    parser.add_argument("--TPM", default=1, type=float, help="TPM cutoff. Defaults to 1")

    args = parser.parse_args()
    #
    cell_type = args.celltype  # spermatocyte or enterocyte
    network = 'full'
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')
    TPM_cutoff = args.TPM
    database = args.database
    #
    #
    #

    print('Reading Network')
    rGfile_gpickle = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/DGE_network/DGE_{c}_multilayer_network_TPM_cutoff_{tpm}.gpickle'.format(db = database, tpm = TPM_cutoff, c = cell_type)
    G = nx.read_gpickle(rGfile_gpickle)

    #
    # Thresholding
    #
    print('Removing edges (remove edges < threshold)')
    print('Number of edges: {:,d}'.format(G.number_of_edges()))
    edges_to_remove = [(i, j) for i, j, d in G.edges(data=True) if (
        ((d.get('type') == 'intra') and (d.get('weight') < threshold))
    )]
    G.remove_edges_from(edges_to_remove)
    print('Number of edges: {:,d}'.format(G.number_of_edges()))

    # Removal of Nodes with only type=='cross' edges
    print('Removing nodes (with only cross edges)')
    print('Number of nodes: {:,d}'.format(G.number_of_nodes()))
    remove_isolates_nodes = []
    for node in G.nodes():
        if not any([True for i, j, d in G.edges(node, data=True) if d['type'] == 'intra']):
            remove_isolates_nodes.append(node)

    G.remove_nodes_from(remove_isolates_nodes)

    print('Number of nodes: {:,d}'.format(G.number_of_nodes()))

    ##
    # Export
    ##
    print('Exporting')
    wGfile_gpickle = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/DGE_network/Threshold/DGE_{c}_multilayer_network_TPM_cutoff_{tpm}_thr{t}.gpickle'.format(db = database, tpm = TPM_cutoff, c = cell_type,t = threshold)
    ensurePathExists(wGfile_gpickle)
    nx.write_gpickle(G, wGfile_gpickle)
    wGfile_graphml = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/DGE_network/Threshold/DGE_{c}_multilayer_network_TPM_cutoff_{tpm}_thr{t}.graphml'.format(db = database, tpm = TPM_cutoff, c = cell_type,t = threshold)
    ensurePathExists(wGfile_graphml)
    nx.write_graphml(G, wGfile_graphml)
