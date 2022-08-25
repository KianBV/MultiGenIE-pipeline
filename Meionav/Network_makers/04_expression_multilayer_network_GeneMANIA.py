#Date: 9.8.2022.
#Written by: Kian Bigović Villi
#The script inputs a geneMANIA edgelist of a set of organismns, changes all of the Gene IDs to ENSEMBL gene ID, and outputs a multilayered network of the GeneMANIA interactions baes on the input gene expression data.
#The full script takes around 4ish minutes on my computer (16GB RAM memory, 2,5 GHz CPU, 512 GB SSD) - This is much faster than making the full network and filtering it later.
import numpy as np
import pandas as pd
import networkx as nx
from utils import ensurePathExists, open_undefined_last_column_files
import time 
from itertools import chain, product
import argparse
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

if __name__ == '__main__':
	#
	#Arguments
	#
	parser = argparse.ArgumentParser()
	celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
	parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Defaults to spermatocyte")
	parser.add_argument("--TPM", default=1, type=float, help="TPM cutoff. Defaults to 1")
	potential_organisms = ["Homo_sapiens","Mus_musculus","Drosophila_melanogaster","Rattus_norvegicus","Danio rerio", "Caenorhabditis_elegans", "Saccharomyces_cerevisiae"]
	parser.add_argument("--organisms", nargs='+', default = ["Homo_sapiens", "Mus_musculus", "Drosophila_melanogaster"],type = str, choices = potential_organisms, help = "All the organisms to build the network. Defaults to human, mouse and fruit fly")
	potential_orthologs = ['metazoa','opisthokonta', 'Metazoa', 'Opisthokonta']
	parser.add_argument("--orthologs", default = "metazoa", type = str.lower, choices = potential_orthologs, help = "Defines the Eggnogg dataset for the gene orthology search. Defaults to metazoa, but if yeast or other fungi are included, use opisthokonta.")
	args = parser.parse_args()
	#
	#
	#
	cell_type = args.celltype
	TPM_cutoff = args.TPM
	all_organisms = args.organisms
	potential_tax_ids =['9606', '10090','7227', '10116', '7955', '6239']
	dict_tax_id = dict(zip(potential_organisms,potential_tax_ids))
	organism_tax_ids = list(map(dict_tax_id.get, all_organisms))
	orthologs_id = "33208" if args.orthologs == "metazoa" else "33154"
	#
	#Network creation
	#
	start = time.time()

	def id_getter(organism):
		#combines all the id files into a single file and indexes the organism
		df_id = pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/GeneMANIA/ID_files/raw/{org}_id.txt".format( org = organism))
		df_id['organism'] = organism
		return df_id

	#Create the ID tables and indices needed to change all IDS to ENSEMBL stable gene IDs
	print("Starting with ID files")
	df_id = pd.concat([id_getter(i) for i in all_organisms])
	df_ENS_id = df_id.loc[(df_id["Source"]=="Ensembl Gene ID")].drop_duplicates(subset=["Preferred_Name"])
	#Make the preferred name vs ENS id dictionary
	dict_ENS_id = dict(df_ENS_id.loc[:, ["Preferred_Name", "Name"]].values)
	#Make a df to link ENSEMBL Gene ID to all other IDs
	df_id["Preferred_Name"] = df_id["Preferred_Name"].map(dict_ENS_id)
	print("Done with ID files")



	def global_network_maker(organism):
		#Creates a network based from the organism inputted. It requires the GeneMANIA "identifier_mappings.txt" and "COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt" file for the appropriate organism
		#Inputs the organism name in a "Genus_description" format eg. "Homo_sapiens", adds nodes and edges to a global network variable "G"
		"""
		This function inputs a string containing the latin name of an organism, loads the GeneMANIA data of the associated organism and adds the appropriate nodes and edges
		to a global nx.Graph variable

		The organism is inputed in a the genus_descriptor format eg. Homo_sapiens

		It requires pathing to organism_whole_network.txt GeneMANIA data.
		It also requires an input file containing gene expression data, this is expected to be a single column of ENSEMBL stable Gene IDs.
		"""
		global dict_ENS_id
		#Creating the whole full genome network
		print("Creating the full network for {:}".format(organism))
		#Input the graph edge data and rename the weight variable (case sensitivity)
		df_edges=pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/GeneMANIA/whole_networks/{org}_whole_network.txt".format( org = organism)).rename({"Weight" : "weight"}, axis = 1)

		"""#Normalize the weights between 0 and 1
					max_value = df_edges['weight'].max()
					min_value = df_edges['weight'].min()
					df_edges['weight'] = (df_edges['weight']-min_value) / (max_value - min_value)
				"""
		#Changing the ID to the ENS ID
		df_edges["Gene_A"] = df_edges["Gene_A"].map(df_ENS_id.set_index('Preferred_Name')['Name'])
		df_edges["Gene_B"] = df_edges["Gene_B"].map(df_ENS_id.set_index('Preferred_Name')['Name'])


		#Input the DGE data
		df_dge = pd.read_csv("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/Test_data/{c}/{org}-FPKM-{c}.csv.gz".format(org = organism, c = cell_type))
		df_dge = df_dge.dropna(subset=["id_gene"]).loc[df_dge['TPM'] > TPM_cutoff,:]
		#filter the edges for the dge genes
		df_edges = df_edges.loc[df_edges['Gene_A'].isin(df_dge['id_gene']) & df_edges['Gene_B'].isin(df_dge['id_gene'])]
		#Gets the list of of the nodes in the network
		print("Adding nodes")
		all_nodes = pd.unique(df_edges[['Gene_A', 'Gene_B']].values.ravel( 'K'))
		#Adds the nodes to the graph
		G.add_nodes_from(all_nodes, layer = organism)
		print("Added nodes, starting with edges")
		
		#Add edge type 
		df_edges["type"]='intra'
		#Add organism tag
		df_edges["layer"] = organism
		#Index the edges
		df_edges = df_edges.set_index(['Gene_A','Gene_B'])
		edge_idx = df_edges.index.to_list()
		#Put the attributes inside a dictionary
		edge_attr = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_edges.to_dict('records')]
		#combine the edges 
		G.add_edges_from([(i, j, d) for (i, j), d in zip(edge_idx, edge_attr)])
		print(G)

	G = nx.Graph()

	for i in all_organisms:
		global_network_maker(i)
	print("All organisms added, starting with cross edges")
	#Create a df with all nodes and their respecitve organism from the graph data
	all_graph_nodes = pd.DataFrame(G.nodes(data = 'layer'), columns = ['node', 'organism'])



	#Make a df that links ENSEBL Protein and Gene id
	df_prot_id = df_id.loc[(df_id["Source"]=="Ensembl Protein ID")].rename( {"Name":"Protein_id"}, axis = 1)
	#Combine all of the ID pairs in a single dataframe and keep only those that are in our network
	df_prot_id = df_prot_id.loc[df_prot_id['Preferred_Name'].isin(all_graph_nodes['node'])]	
	#Map the EGGnog id
	df_prot_id['organism'] = df_prot_id['organism'].map(dict_tax_id)
	
	#Create the EGGnogg ID : ENSEMBL gene ID dictionary for mapping later
	dict_gene_prot = dict(zip(df_prot_id['organism'] + '.' + df_prot_id['Protein_id'], df_prot_id['Preferred_Name']))

	#Create a set of all of the EGGnogg IDs in our networks
	Egg_id_set = set(df_prot_id['organism'] + '.' + df_prot_id['Protein_id'])

	print("Done with IDs")
	print("Starting with EGGnogg members file")
	df_Egg = open_undefined_last_column_files(
	    "C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/EGGnogg_files/{id}_members.tsv.gz".format(id = orthologs_id),
	    n_fixed_cols=5,
	    names = ['family', 'id_eggnog', '_1', '_2', 'aliases', 'species'],
	    nrows = None
	    )
	#For keeping the columns and species we need
	df_Egg = df_Egg.set_index('id_eggnog')['aliases']
	wanted_species = frozenset(organism_tax_ids)


	def select_by_species(text, keeplist):
	# Only keep genes from species we are interested in (lower the search space)
	    return [i for i in text.split(',') if i.split('.', 1)[0] in keeplist]

	df_Egg = df_Egg.apply(select_by_species, args=(wanted_species,))

	def select_by_at_least_one_match(ilist, keeplist):
	        # Only keep genes that are found in any of our gene list (lower the search space)
	        wanted_genes = [i for i in ilist if i in keeplist]
	        return wanted_genes if len(wanted_genes) >= 1 else None

	df_Egg = df_Egg.apply(select_by_at_least_one_match, args=(Egg_id_set, ))
	df_Egg = df_Egg.dropna()


	def ortholog_edge_generator(row):
		#Takes in a series of EGGnogg orthogroups, filters them for the genes present in our network, and outputs their edges

		#Filter the ortholog groups by the genes within our network
		orthologs = [x for x in row if x in Egg_id_set]
		orthologs = list(map(dict_gene_prot.get, orthologs))
		#Make a list of combinations for the given orthogroup, note every gene has an edge with itself (will be filtered out later, it isn't as time intensive to keep them)
		edges = list(product(orthologs, repeat = 2))
		return edges

	print("Computed all the cross-edges")
	#Creates the edges
	all_ort_edges = df_Egg.apply(ortholog_edge_generator)
	#Create a list of edge tuples
	all_ort_edges = all_ort_edges.tolist()
	all_ort_edges = list(chain(*all_ort_edges))
	#Filter out the ones where the gene has an interacton with itself
	all_ort_edges = [i for i in all_ort_edges if i[0] != i[1]]

	#Since multiple orthologs from the same organism might be present in the same orthogroup, the list should be filtered for the interactions between genes of the same organism
	
	#Turn it into a dataframe
	df_all_ort_edges = pd.DataFrame(all_ort_edges, columns = ['Gene_A','Gene_B'])

	#Turn the gene into organism
	df_all_ort_edges['Gene_A_c']=df_all_ort_edges['Gene_A'].map(df_ENS_id.set_index('Preferred_Name')['organism'])
	df_all_ort_edges['Gene_B_c']=df_all_ort_edges['Gene_B'].map(df_ENS_id.set_index('Preferred_Name')['organism'])
	
	#Filter out the interactons where both genes belong to the same organism
	df_all_ort_edges = df_all_ort_edges.loc[df_all_ort_edges['Gene_A_c'] != df_all_ort_edges['Gene_B_c']]
	
	#Turn it back into a tuple list for easy add to the graph
	all_ort_edges = list(zip(df_all_ort_edges['Gene_A'],df_all_ort_edges['Gene_B']))

	#
	#Exporting
	#



	print("Adding the cross-edges")
	#Add the edges to our graph
	G.add_edges_from(all_ort_edges, type='cross', weight = 1)
	print("Added all the cross-edges")
	print("Full network complete")
	print(G)

	print('Exporting the full newtwork')
	#Write in gpickle
	path_gpickle = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/GeneMANIA/expression_network/expression_{c}_multilayer_network_TPM_cutoff_{tpm}.gpickle'.format( tpm = TPM_cutoff, c = cell_type)
	ensurePathExists(path_gpickle)
	nx.write_gpickle(G,path_gpickle)
	print("Exported the .gpickle")
	"""
	path_edgelist = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/GeneMANIA/expression_network/expression_{c}_multilayer_network_TPM_cutoff_{tpm}.edgelist'.format( tpm = TPM_cutoff, c = cell_type)
	ensurePathExists(path_gpickle)
	nx.write_edgelist(G, path_edgelist)
	print("Exported the .edgelist")

	path_graphml = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}GeneMANIA/expression_network/expression_{c}_multilayer_network_TPM_cutoff_{tpm}.graphml'.format( tpm = TPM_cutoff, c = cell_type)
	ensurePathExists(path_graphml)
	nx.write_graphml(G, path_graphml)
	print("Exported the .graphml")
	"""
	print("Done with exporting")
	end = time.time()
	print("Done, runtime : {:}".format(end-start))
