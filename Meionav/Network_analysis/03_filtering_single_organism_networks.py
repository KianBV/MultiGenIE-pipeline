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
	#
	parser = argparse.ArgumentParser()
	celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
	parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Defaults to spermatocyte")
	databases = ["GeneMANIA", "IntAct", "CPDB", "Mentha", "StringDB"]
	parser.add_argument("--database", default="StringDB", type=str, choices = databases, help="Database. Defaults to StringDB")
	parser.add_argument("--TPM", default=1, type=float, help="TPM cutoff. Defaults to 1")
	potential_organisms = ["Homo_sapiens","Mus_musculus","Drosophila_melanogaster","Rattus_norvegicus","Danio rerio", "Caenorhabditis_elegans", "Saccharomyces_cerevisiae"]
	parser.add_argument("--organisms", nargs='+', default = ["Homo_sapiens"],type = str, choices = potential_organisms, help = "For which organisms do you want to filter a network. Defaults to human and mouse.")

	args = parser.parse_args()
    #
    #
    #
	cell_type = args.celltype
	all_organisms = args.organisms
	TPM_cutoff = args.TPM
	database = args.database
    #
    #Import the network
    #
	print('Reading network')
	G = nx.read_gpickle('C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/expression_network/expression_{c}_multilayer_network_TPM_cutoff_{tpm}.gpickle'.format(db = database, tpm = TPM_cutoff, c = cell_type))

	for organism in all_organisms:
		#Create a subgraph of that organism
		SG = get_network_layer(G, layer = organism)
			
		path_gpickle = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/expression_network/Single_organism/{org}/expression_{c}_{org}_network_TPM_cutoff_{tpm}.gpickle'.format(db = database, tpm = TPM_cutoff, c = cell_type, org = organism)
		ensurePathExists(path_gpickle)
		nx.write_gpickle(SG,path_gpickle)
		print("Exported the {o} .edgelist".format(o = organism))
		
		path_graphml = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/expression_network/Single_organism/{org}/expression_{c}_{org}_network_TPM_cutoff_{tpm}.graphml'.format(db = database, tpm = TPM_cutoff, c = cell_type, org = organism)
		ensurePathExists(path_graphml)
		nx.write_graphml(SG, path_graphml)
		print("Exported the .graphml")

	print('Done with {db}'.format(db = database))