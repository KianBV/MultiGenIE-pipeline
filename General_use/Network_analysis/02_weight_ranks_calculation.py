#Date: 19.8.2022.
#Written by: Kian Bigović Villi
#Written in: Python 3.9.12
#Takes in multilayer graphs and outputs the weight ranks
import pandas as pd
import networkx as nx
from utils import ensurePathExists, get_network_layer
import argparse
pd.options.mode.chained_assignment = None  # Safe to turn off


if __name__ == '__main__':
	#
	#Arguments
	#
	parser = argparse.ArgumentParser()
	parser.add_argument("--input_file",type = str, help = "Path to and file name of the .gpickle file for the given inout network, eg. folder1/folder2/folder3/expression_multilayer_network.gpickle")
	parser.add_argument("--output_path", default = "C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/StringDB/expression_network/expression_multilayer_network",type = str, help = "Path to and name for the output file")
	args = parser.parse_args()
	#
	#
	#
	input_path = args.input_file
	output_path = args.output_path
	#
	#
	#Stats
	#
	
	print("Loading network")
	#Load the graph gickle
	G = nx.read_gpickle(input_path)
	#Get the edge attributes and put the in a DataFrame
	edge_data = [i[2] for i in G.edges(data = True)]
	weights = pd.DataFrame(edge_data)
	#Filter the intra edges
	weights_intra = weights.loc[weights['type'] == "intra"].drop('type', axis = 1)

	#Gte the nodes and their layer
	nodes = nx.get_node_attributes(G, 'layer')
	
	#List of organisms in the layer
	organism_list = set(nodes.values())
	
	#Loop over it to get the stats
	for organism in organism_list:
		df_weights = weights.loc[weights['layer'] == organism]

		#Export it
		path_weights = output_path + "{org}_weights.txt".format(org = organism)
		ensurePathExists(path_weights)
		df_weights.to_csv(path_weights, index = False)

		#Calculate the rank
		df_weights['rank'] = df_weights['weight'].rank(method = 'min', ascending = False)
		df_weights = df_weights.loc[:, ['weight', 'rank']].drop_duplicates(subset = ['weight','rank'])
		
		#Export it
		path_ranks = output_path + "{org}_ranks_weights.txt".format(org = organism)
		ensurePathExists(path_ranks)
		df_weights.to_csv(path_ranks, index = False)
		print("Done for {org} ".format(org = organism))
