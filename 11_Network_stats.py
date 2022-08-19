import pandas as pd
import networkx as nx
from utils import ensurePathExists, get_network_layer

#Possible databases are GeneMANIA, IntAct, CPDB, Mentha, StringDB
all_databases = ['GeneMANIA', 'IntAct', 'CPDB', 'Mentha','StringDB']
#The cuttoff for the gene being expressed
TPM_cutoff = 1
#Potential cells are 'enterocyte', 'neuron', 'muscle', 'spermatogonia', 'spermatocyte', 'spermatid'
all_cell_types = ['enterocyte', 'neuron', 'muscle', 'spermatogonia', 'spermatocyte', 'spermatid']
all_organisms = ['Homo_sapiens', 'Mus_musculus', 'Drosophila_melanogaster']




report = pd.DataFrame(columns = ['database', 'organism', 'cell_type', 'nodes', 'edges', 'density'])

for database in all_databases:
	for cell_type in all_cell_types:
		print("Starting {db} {c}".format(db = database, c = cell_type))
		#Load the graph
		G = nx.read_gpickle('C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/DGE_network/DGE_{c}_multilayer_network_TPM_cutoff_{tpm}.gpickle'.format(db = database, tpm = TPM_cutoff, c = cell_type))
		#Gte the nodes and their layer
		nodes = nx.get_node_attributes(G, 'layer')
		#Get a df of nodes and layers
		df_nodes = pd.DataFrame.from_dict(nodes, orient='index', columns = [ 'organism'])
		#Get a list of organisms in the network
		organism_list = df_nodes['organism'].unique()
		#Loop over it to get the stats
		for organism in organism_list:
			SG = get_network_layer(G, layer = organism)
			nodes = SG.number_of_nodes()
			edges = SG.number_of_edges()
			density = nx.density(SG)
			#Append to the report
			report_line = pd.DataFrame.from_dict({'database' : [database], 
				'organism' : [organism], 
				'cell_type' : [cell_type],
				'nodes' : [nodes], 
				'edges' : [edges],
				'density': [density]})
			report = pd.concat([report, report_line], ignore_index = True)
report.to_csv('C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/Network_stats/Network_stats.txt', index = False)