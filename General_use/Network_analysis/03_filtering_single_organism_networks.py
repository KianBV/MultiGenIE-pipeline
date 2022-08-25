#Filter out the organism networks from a multilayered network
import pandas as pd
import networkx as nx
import numpy as np
from utils import ensurePathExists, get_network_layer
import argparse

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

if __name__ == '__main__':

	#
	# Args
	parser = argparse.ArgumentParser()
	parser.add_argument("--input_file",type = str, help = "Path to and file name of the .gpickle file for the given inout network, eg. folder1/folder2/folder3/expression_multilayer_network.gpickle")
	parser.add_argument("--output_path", default = "C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/StringDB/expression_network/expression_multilayer_network",type = str, help = "Path to and name for the output file")
	parser.add_argument("--organisms", nargs='+', default = ["Homo_sapiens"],type = str, help = "For which organisms do you want to filter a network. Defaults to human and mouse.")
	parser.add_argument("--output_type", nargs = '+',default = ["gpickle", "graphml"], type = str, help = "Which file format should the output files be." )
	args = parser.parse_args()
    #
    #
    #
	all_organisms = args.organisms
	input_path = args.input_file
	output_path = args.output_path
	output_files = set(args.output_type)
    #
    #Import the network
    #
	print('Reading network')
	G = nx.read_gpickle(input_path)

	for organism in all_organisms:
		#Create a subgraph of that organism
		SG = get_network_layer(G, layer = organism)
			

		#Write in gpickle
		if "gpickle" in output_files:
			path_gpickle = output_path + organism + ".gpickle"
			ensurePathExists(path_gpickle)
			nx.write_gpickle(SG,path_gpickle)
			print("Exported the .gpickle")
		#Write in edgelist
		if "edgelist" in output_files:
			path_edgelist = output_path + organism + '.edgelist'
			ensurePathExists(path_gpickle)
			nx.write_edgelist(SG, path_edgelist)
			print("Exported the .edgelist")
		#Write in graphml
		if "graphml" in output_files:
			path_graphml = output_path + organism + '.graphml'
			ensurePathExists(path_graphml)
			nx.write_graphml(SG, path_graphml)
			print("Exported the .graphml")
		print("Done with {org}".format(org = organism))


	print('Done')