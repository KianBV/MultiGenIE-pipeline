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
    parser.add_argument("--input_file",type = str, help = "Path to and file name of the .gpickle file for the given inout network, eg. folder1/folder2/folder3/expression_multilayer_network.gpickle")
    parser.add_argument("--output_path", default = "C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/StringDB/expression_network/expression_multilayer_network",type = str, help = "Path to and name for the output file")
    parser.add_argument("--output_type", nargs = '+',default = ["gpickle", "graphml"], type = str, help = "Which file format should the output files be." )
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    args = parser.parse_args()
    #
    input_path = args.input_file
    output_path = args.output_path
    output_files = set(args.output_type)
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')

    #
    #
    #

    print('Reading Network')
    rGfile_gpickle = input_path
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
    #Write in gpickle
    if "gpickle" in output_files:
        path_gpickle = output_path + "thr{thr}_.gpickle".format(thr = threshold_str)
        ensurePathExists(path_gpickle)
        nx.write_gpickle(G,path_gpickle)
        print("Exported the .gpickle")
    #Write in edgelist
    if "edgelist" in output_files:
        path_edgelist = output_path + 'thr{thr}_.edgelist'.format(thr = threshold_str)
        ensurePathExists(path_gpickle)
        nx.write_edgelist(G, path_edgelist)
        print("Exported the .edgelist")
    #Write in graphml
    if "graphml" in output_files:
        path_graphml = output_path + 'thr{thr}_.graphml'.format(thr = threshold_str)
        ensurePathExists(path_graphml)
        nx.write_graphml(G, path_graphml)
        print("Exported the .graphml")
