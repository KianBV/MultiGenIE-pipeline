#Date: 13.8.2022.
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
	potential_tax_ids =[9606, 10090,7227, 10116, 7955, 6239]
	dict_tax_id = dict(zip(potential_organisms,potential_tax_ids))
	organism_tax_ids = list(map(dict_tax_id.get, all_organisms))
	orthologs_id = "33208" if args.orthologs == "metazoa" else "33154"
	#
	#Network creation
	#
	start = time.time()

	print("Starting")
	G = nx.Graph()
	print("Inputing full Mentha edge data")
	#Input egde data
	df_data = pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/Mentha/Full_mentha", sep = ";")
	#Filter out the intra edges of our organims of interest
	df_edges = df_data.loc[df_data['Taxon A'].isin(organism_tax_ids) & df_data['Taxon B'].isin(organism_tax_ids) & (df_data['Taxon A'] == df_data['Taxon B']) ,:]
	#Take the edge data into it
	df_edges = df_edges.loc[:, ['Protein A', 'Protein B', 'Taxon A', 'Score']].rename(columns = {'Protein A' : 'Gene_A', 'Protein B' : 'Gene_B', 'Taxon A' : 'organism', 'Score' : 'weight'})
	#organism id to name
	dict_organism = dict(zip(organism_tax_ids, all_organisms))
	df_edges['organism'] = df_edges['organism'].map(dict_organism)
	print("Filtered edge data by organisms")
	print("Inputing IDs")
	def id_getter(organism):
		#Fetch the ENSEMBL ID vs UniProt ID files
		if organism in ["Drosophila_melanogaster","Caenorhabditis_elegans"]:
			#Drosophila and Caenorhabditis done have protein ID versions
			df_id = pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/IDs/{x}_ENSEMBL_ID.txt".format(x = organism))
			df_id = df_id.rename(columns = {'Gene stable ID' : 'Preferred_Name', 'Protein stable ID' : 'Protein_id', 'UniProtKB Gene Name ID': 'Uni_Protein_id', 'UniProtKB/Swiss-Prot ID' : "Uni_Protein_id"})
			df_id['organism'] = organism
		else:
			#The rest have them, and they are used in Eggnogg
			df_id = pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/IDs/{x}_ENSEMBL_ID.txt".format(x = organism))
			df_id = df_id.rename(columns = {'Gene stable ID' : 'Preferred_Name', 'Protein stable ID' : 'Protein_id', 'UniProtKB Gene Name ID': 'Uni_Protein_id'})
			df_id['organism'] = organism
		return df_id

	#Combine all ID's into a single file
	df_all_id = pd.concat([id_getter(i) for i in all_organisms]).dropna(subset = ['Uni_Protein_id'])

	#Mapping df
	df_map_id = df_all_id.drop_duplicates(subset = "Uni_Protein_id")

	#Map it to the edge file
	df_edges['Gene_A'] = df_edges['Gene_A'].map(df_map_id.set_index('Uni_Protein_id')['Preferred_Name'])
	df_edges['Gene_B'] = df_edges['Gene_B'].map(df_map_id.set_index('Uni_Protein_id')['Preferred_Name'])
	df_edges = df_edges.dropna()
	print("Done with IDs")
	print("Print inputing DEGs")
	#Loads the DEGs
	def deg_getter(organism):
		#Input the DGE data
		df_dge = pd.read_csv("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/Test_data/{c}/{org}-FPKM-{c}.csv.gz".format(org = organism, c = cell_type))
		#Take only the DGE genes
		df_dge = df_dge.dropna(subset=["id_gene"]).loc[df_dge['TPM'] > TPM_cutoff,:]
		return df_dge

	df_all_dge = pd.concat([deg_getter(i) for i in all_organisms])
	#Filter out the DEGs
	df_edges = df_edges.loc[df_edges['Gene_A'].isin(df_all_dge['id_gene']) & df_edges['Gene_B'].isin(df_all_dge['id_gene'])]
	print("Filtered by Degs")
	print("Adding nodes")
	#Get a list of all nodes
	df_all_nodes = pd.concat([df_edges.drop(['Gene_A', 'weight'], axis = 1).rename(columns = {'Gene_B' : 'node'}), df_edges.drop(['Gene_B', 'weight'], axis = 1).rename(columns = {'Gene_A' : 'node'})])
	df_all_nodes = df_all_nodes.drop_duplicates().rename(columns = {'organism' : 'layer'})




	#Adds nodes
	G.add_nodes_from(df_all_nodes['node'])
	#Add the attributes
	attrs = df_all_nodes.set_index('node').T.to_dict()
	nx.set_node_attributes(G, attrs)
	print("Added nodes")
	print("Adding edges")
	#Drop the organism column
	df_edges = df_edges.drop(['organism'], axis = 1)
	df_edges["type"]='intra'
	#Index the edges
	df_edges = df_edges.set_index(['Gene_A','Gene_B'])
	edge_idx = df_edges.index.to_list()
	#Put the attributes inside a dictionary
	edge_attr = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_edges.to_dict('records')]
	#combine the edges 
	G.add_edges_from([(i, j, d) for (i, j), d in zip(edge_idx, edge_attr)])
	print("Added edges")
	print(G)

	print("Setting up EGGnogg IDs")
	all_graph_nodes = df_all_nodes.rename(columns = {'layout' : 'organism'})

	df_prot_id = df_all_id.loc[df_all_id['Preferred_Name'].isin(all_graph_nodes['node'])].copy()
	#Create the organism, EGGnog id dictionary
	organism_tax_ids = [str(i) for i in organism_tax_ids]
	dict_org = dict(zip(all_organisms, organism_tax_ids))
	#Map the EGGnog id
	df_prot_id['organism'] = df_prot_id['organism'].map(dict_org)
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
	
	#remove duplicate mapping
	df_all_id = df_all_id.drop_duplicates(subset = "Preferred_Name")

	#Map the gene to the organism
	df_all_ort_edges['Gene_A_c']=df_all_ort_edges['Gene_A'].map(df_all_id.set_index('Preferred_Name')['organism'])
	df_all_ort_edges['Gene_B_c']=df_all_ort_edges['Gene_B'].map(df_all_id.set_index('Preferred_Name')['organism'])
	
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
	path_gpickle = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/Mentha/expression_network/expression_{c}_multilayer_network_TPM_cutoff_{tpm}.gpickle'.format(tpm = TPM_cutoff, c = cell_type)
	ensurePathExists(path_gpickle)
	nx.write_gpickle(G,path_gpickle)
	print("Exported the .gpickle")
	"""
	path_edgelist = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/Mentha/expression_network/expression_{c}_multilayer_network_TPM_cutoff_{tpm}.edgelist'.format(tpm = TPM_cutoff, c = cell_type)
	ensurePathExists(path_gpickle)
	nx.write_edgelist(G, path_edgelist)
	print("Exported the .edgelist")

	path_graphml = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/Mentha/expression_network/expression_{c}_multilayer_network_TPM_cutoff_{tpm}.graphml'.format(tpm = TPM_cutoff, c = cell_type)
	ensurePathExists(path_graphml)
	nx.write_graphml(G, path_graphml)
	print("Exported the .graphml")
	"""
	print("Done with exporting")
	end = time.time()
	print("Done, runtime : {:}".format(end-start))