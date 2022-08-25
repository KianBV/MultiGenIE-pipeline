import pandas as pd
import networkx as nx
from utils import ensurePathExists
import argparse


if __name__ == '__main__':
	#
	#Arguments
	#
	parser = argparse.ArgumentParser()
	parser.add_argument("--input_file",type = str, help = "Path to and file name of the .gpickle file for the given inout network, eg. folder1/folder2/folder3/expression_multilayer_network.gpickle")
	parser.add_argument("--output_file", default = "C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/StringDB/expression_network/expression_multilayer_network",type = str, help = "Path to and name for the output file")
	args = parser.parse_args()
	#
	#
	#
	input_path = args.input_file
	output_path = args.output_file
	#
	report = pd.DataFrame(columns = ['organism', 'nodes', 'edges', 'density'])
	print("Loading network")
	#Load the graph
	G = nx.read_gpickle(input_path)
	#Get nodes and edges
	nodes = nx.get_node_attributes(G, 'layer')
	edges = nx.get_edge_attributes(G, 'layer')
	#A list of organisms in the given network
	organism_list = set(nodes.values())
	#Loop over it to get the stats
	for organism in organism_list:
		#Number of nodes and edges
		node_nr = sum(1 for x in nodes.values() if x == organism)
		edge_nr = sum(1 for x in edges.values() if x == organism)
		#calculate denstiy
		density = 2 * edge_nr / (node_nr * (node_nr - 1))
		#Append to the report
		report_line = pd.DataFrame.from_dict({
			'organism' : [organism], 
			'nodes' : [node_nr], 
			'edges' : [edge_nr],
			'density': [density]})
		report = pd.concat([report, report_line], ignore_index = True)
	ensurePathExists(output_path)
	report.to_csv(output_path, index = False)