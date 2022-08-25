#Date: 19.8.2022.
#Written by: Kian Bigović Villi
#Written in: Python 3.9.12
#Takes in multilayer graphs and outputs the weight ranks
import pandas as pd
import networkx as nx
from utils import ensurePathExists, get_network_layer
import argparse


if __name__ == '__main__':
	#
	# Args
	#
	parser = argparse.ArgumentParser()
	potential_celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
	parser.add_argument("--celltype", nargs = '+', default = ['enterocyte', 'neuron', 'muscle', 'spermatogonia', 'spermatocyte', 'spermatid'], type = str, choices = potential_celltypes, help = "Cell type. Defaults to all cell types.")
	potential_databases = ["GeneMANIA", "IntAct", "CPDB", "Mentha", "StringDB"]
	parser.add_argument("--database", nargs = '+', default = ['GeneMANIA', 'IntAct', 'CPDB', 'Mentha','StringDB'], type = str, choices = potential_databases, help = "Database. Defaults to all databases")
	parser.add_argument("--TPM", default = 1, type = float, help = "TPM cutoff. Defaults to 1")
	args = parser.parse_args()
    #
	all_cell_types = args.celltype  # spermatocyte or enterocyte
	TPM_cutoff = args.TPM
	all_databases = args.database
	#
	#Stats
	#
	for database in all_databases:
		for cell_type in all_cell_types:
			print("Starting {db} {c}".format(db = database, c = cell_type))
			#Load the graph gickle
			G = nx.read_gpickle('C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/expression_network/expression_{c}_multilayer_network_TPM_cutoff_{tpm}.gpickle'.format(db = database, tpm = TPM_cutoff, c = cell_type))
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

				
				#Add attributes to the file
				df_weights['databases'] = database
				df_weights['cell_type'] = cell_type

				#Export it
				path_weights = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/Network_stats/weights/TESTweights_{org}_{c}_network_{db}.txt'.format(org = organism, c = cell_type, db = database)
				ensurePathExists(path_weights)
				df_weights.to_csv(path_weights, index = False)

				#Calculate the rank
				df_weights['rank'] = df_weights['weight'].rank(method = 'min', ascending = False)
				df_weights = df_weights.loc[:, ['weight', 'databases', 'cell_type', 'rank']].drop_duplicates(subset = ['weight','rank'])
				
				#Export it
				path_weights = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/Network_stats/weights_ranks/TESTweights_ranks_{org}_{c}_network_{db}.txt'.format(org = organism, c = cell_type, db = database)
				ensurePathExists(path_weights)
				df_weights.to_csv(path_weights, index = False)
				print("Done for {org} {db} {c}".format(org = organism, c = cell_type, db = database))
