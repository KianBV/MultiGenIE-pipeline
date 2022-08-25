#Date: 15.8.2022.
#Written by: Kian Bigović Villi
#Written in: Python 3.9.12
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
	potential_organisms = ["Homo_sapiens","Mus_musculus","Drosophila_melanogaster","Rattus_norvegicus","Danio rerio", "Caenorhabditis_elegans", "Saccharomyces_cerevisiae"]
	parser.add_argument("--organisms", nargs='+', default = ["Homo_sapiens", "Mus_musculus", "Drosophila_melanogaster"],type = str, choices = potential_organisms, help = "All the organisms to build the network. Defaults to human, mouse and fruit fly")
	potential_orthologs = ['metazoa','opisthokonta', 'Metazoa', 'Opisthokonta']
	parser.add_argument("--orthologs", default = "metazoa", type = str.lower, choices = potential_orthologs, help = "Defines the Eggnogg dataset for the gene orthology search. Defaults to metazoa, but if yeast or other fungi are included, use opisthokonta.")
	parser.add_argument("--expression_path",type = str, help = "Path to the file that contains the gene expression data. Each file inside that folder should be renamed to x_expression_data.txt.gz, where x is the latin name of given organism.")
	parser.add_argument("--output_file", default = "C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/StringDB/expression_network/expression_multilayer_network",type = str, help = "Path to and name of the output files")
	potential_outputs = ["gpickle", "edgelist", "graphml"]
	parser.add_argument("--output_type", nargs = '+',default = ["gpickle", "graphml"], type = str, help = "Which file format should the output files be." )
	args = parser.parse_args()
	#
	#
	#
	all_organisms = args.organisms
	potential_tax_ids =['9606', '10090','7227', '10116', '7955', '6239']
	dict_tax_id = dict(zip(potential_organisms,potential_tax_ids))
	organism_tax_ids = list(map(dict_tax_id.get, all_organisms))
	orthologs_id = "33208" if args.orthologs == "metazoa" else "33154"
	expression_path = args.expression_path
	output_path = args.output_file
	output_files = set(args.output_type)
	#
	#Network creation
	#
	start = time.time()
	def global_network_maker(organism):
		"""
		This function inputs a string containing the latin name of an organism, loads the StringDB data of the associated organism, adds the appropriate nodes and edges
		to a global nx.Graph variable, and outputs the ids used in mapping the file.

		The organism is inputed in a the genus_descriptor format eg. Homo_sapiens

		It requires pathing to organism.protein.aliases.vX.txt.gz and organism.protein.links.vX.txt.gz files, with X being the StringDB version
		It also requires an input file containing gene expression data, this is expected to be a single column of ENSEMBL stable Gene IDs.
		"""
		print('Starting with {org}'.format(org = organism))
		#Import STRING IDs
		df_id = pd.read_table('C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/StringDB/{org}.protein.aliases.v11.5.txt.gz'.format(org = organism)).dropna()
		#Keep only the ENSEMBL gene identifiers
		df_id = df_id.loc[df_id['source'] == 'Ensembl_gene']
		#For later identification
		df_id['organism'] = organism
		print("Processing StringDB data")
		#Import string data
		df_edges = pd.read_csv("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/StringDB/{org}.protein.links.v11.5.txt.gz".format(org = organism), sep=' ')
		#Calculate the score
		df_edges['combined_score'] = df_edges['combined_score']/1000
		#Import the Expression data and filter it
		df_dge = pd.read_csv(expression_path + "/{org}_expression_data.txt.gz".format(org = organism))
		#Remap the IDs
		df_edges['protein1'] = df_edges['protein1'].map(df_id.set_index('#string_protein_id')['alias'])
		df_edges['protein2'] = df_edges['protein2'].map(df_id.set_index('#string_protein_id')['alias'])

		#Rename the columns
		df_edges = df_edges.rename(columns = {'protein1' : 'Gene_A', 'protein2' : 'Gene_B', 'combined_score' : 'weight'})
		#Remove those not expressed
		df_edges = df_edges.loc[df_edges['Gene_A'].isin(df_dge['id_gene']) & df_edges['Gene_B'].isin(df_dge['id_gene'])]
		df_edges['organism'] = organism
		df_all_nodes = pd.concat([df_edges.drop(['Gene_A', 'weight'], axis = 1).rename(columns = {'Gene_B' : 'node'}), df_edges.drop(['Gene_B', 'weight'], axis = 1).rename(columns = {'Gene_A' : 'node'})])
		df_all_nodes = df_all_nodes.drop_duplicates().rename(columns = {'organism' : 'layer'})
		G.add_nodes_from(df_all_nodes['node'])
		#Add the attributes
		attrs = df_all_nodes.set_index('node').T.to_dict()
		nx.set_node_attributes(G, attrs)
		print("Added nodes")
		print("Adding edges")

		#Add type
		df_edges["type"]='intra'
		df_edges["layer"] = organism
		#Index the edges
		df_edges = df_edges.set_index(['Gene_A','Gene_B'])
		edge_idx = df_edges.index.to_list()
		#Put the attributes inside a dictionary
		edge_attr = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_edges.to_dict('records')]
		#combine the edges 
		G.add_edges_from([(i, j, d) for (i, j), d in zip(edge_idx, edge_attr)])
		print("Added edges for {org}".format(org = organism))
		print(G)
		#FIlter the ID to inlcude only the expressed genes
		df_id = df_id.loc[df_id['alias'].isin(df_all_nodes['node'])]
		#Return the ID
		return df_id




	G = nx.Graph()
	#Get all the IDs in one frame
	df_all_id = pd.concat([global_network_maker(i) for i in all_organisms])
	#Make a dictionary for linking eggnogg/String ids with Gene IDs
	dict_gene_prot = dict(zip(df_all_id['#string_protein_id'], df_all_id['alias']))
	#Create a set of all interesting proteins
	Egg_id_set = set(df_all_id['#string_protein_id'])
	
	#
	# Making the multilayer network
	#

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
	df_all_ort_edges['organism_A']=df_all_ort_edges['Gene_A'].map(df_all_id.set_index('alias')['organism'])
	df_all_ort_edges['organism_B']=df_all_ort_edges['Gene_B'].map(df_all_id.set_index('alias')['organism'])
	
	#Filter out the interactons where both genes belong to the same organism
	df_all_ort_edges = df_all_ort_edges.loc[df_all_ort_edges['organism_A'] != df_all_ort_edges['organism_B']]
	
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
	if "gpickle" in output_files:
		path_gpickle = output_path + ".gpickle"
		ensurePathExists(path_gpickle)
		nx.write_gpickle(G,path_gpickle)
		print("Exported the .gpickle")
	#Write in edgelist
	if "edgelist" in output_files:
		path_edgelist = output_path + '.edgelist'
		ensurePathExists(path_gpickle)
		nx.write_edgelist(G, path_edgelist)
		print("Exported the .edgelist")
	#Write in graphml
	if "graphml" in output_files:
		path_graphml = output_path + '.graphml'
		ensurePathExists(path_graphml)
		nx.write_graphml(G, path_graphml)
		print("Exported the .graphml")

	print("Done with exporting")
	end = time.time()
	print("Done, runtime : {:}".format(end-start))