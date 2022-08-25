import pandas as pd
import networkx as nx
from utils import ensurePathExists, get_network_layer
import argparse


if __name__ == '__main__':
	#
	#Arguments
	#
	parser = argparse.ArgumentParser()
	celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
	parser.add_argument("--celltype", nargs='+', default='spermatocyte', type=str, choices=celltypes, help="All the cell types you want. Defaults to spermatocyte")
	parser.add_argument("--TPM", default=1, type=float, help="TPM cutoff. Defaults to 1")
	databases = ["GeneMANIA", "IntAct", "CPDB", "Mentha", "StringDB"]
	parser.add_argument("--databases",nargs='+', default="StringDB", type=str, choices = databases, help="Database. Defaults to StringDB")
	args = parser.parse_args()
	#
	#
	#
	all_cell_types = args.celltype
	TPM_cutoff = args.TPM
	all_databases = args.databases
	#
	report = pd.DataFrame(columns = ['database', 'organism', 'cell_type', 'nodes', 'edges', 'density'])

	for database in all_databases:
		for cell_type in all_cell_types:
			print("Starting {db} {c}".format(db = database, c = cell_type))
			#Load the graph
			G = nx.read_gpickle('C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/expression_network/expression_{c}_multilayer_network_TPM_cutoff_{tpm}.gpickle'.format(db = database, tpm = TPM_cutoff, c = cell_type))
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
				report_line = pd.DataFrame.from_dict({'database' : [database], 
					'organism' : [organism], 
					'cell_type' : [cell_type],
					'nodes' : [node_nr], 
					'edges' : [edge_nr],
					'density': [density]})
				report = pd.concat([report, report_line], ignore_index = True)
	report.to_csv('C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/Network_stats/Network_stats_output.txt', index = False)