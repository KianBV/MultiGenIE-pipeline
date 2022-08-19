#Filter out the organism networks from a multilayered network
import pandas as pd
import networkx as nx
import numpy as np
from utils import ensurePathExists

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

#Possible databases are GeneMANIA, IntAct, CPDB, Mentha, StringDB
database = 'StringDB'
#The cuttoff for the gene being expressed
TPM_cutoff = 1
#Potential cells are 'enterocyte', 'neuron', 'muscle', 'spermatogonia', 'spermatocyte', 'spermatid'
cell_type = 'neuron'

# path types are "DGE_{c}_multilayer_network_TPM_cutoff_{tpm}".format(c = cell_type, tpm = TPM_cutoff) for filtered netwoeks or Genome_multilayer_network for full ones
G = nx.read_gpickle('C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/DGE_network/DGE_{c}_multilayer_network_TPM_cutoff_{tpm}.gpickle'.format(db = database, tpm = TPM_cutoff, c = cell_type))
#G = read_gpickle('C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/DGE_network/Genome_multilayer_network.gpickle'.format(db = database))


#Get the nodes
nodes = nx.get_node_attributes(G, 'layer')
#Turn it into a dataframe
df_nodes = pd.DataFrame.from_dict(nodes, orient='index', columns = [ 'organism'])
organism_list = df_nodes['organism'].unique()

for organism in organism_list:
	#Filter the nodes by the appropriate organism
	org_nodes = df_nodes.loc[df_nodes['organism'] == organism].index
	#Create a subgraph of that organism
	SG = nx.subgraph(G, org_nodes)
		
	path_gpickle = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/DGE_network/Single_organism/{o}/DGE_{c}_organism_network_TPM_cutoff_{tpm}.gpickle'.format(db = database, tpm = TPM_cutoff, c = cell_type, o = organism)
	ensurePathExists(path_gpickle)
	nx.write_gpickle(SG,path_gpickle)
	print("Exported the {o} network".format(o = organism))
	
	path_graphml = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/DGE_network/Single_organism/{o}/DGE_{c}_organism_network_TPM_cutoff_{tpm}.graphml'.format(db = database, tpm = TPM_cutoff, c = cell_type, o = organism)
	ensurePathExists(path_graphml)
	nx.write_graphml(SG, path_graphml)
	print("Exported the .graphml")
print('Done with {db}'.format(db = database))