#Date: 19.8.2022.
#Written by: Kian Bigović Villi
#Written in: Python 3.9.12
import numpy as np
import pandas as pd
import networkx as nx
from utils import ensurePathExists, open_undefined_last_column_files
import time 
from itertools import chain, product
pd.options.mode.chained_assignment = None  # default='warn'

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

if __name__ == '__main__':
	#
	#Arguments
	#
	parser = argparse.ArgumentParser()
	potential_organisms = ["Homo_sapiens","Mus_musculus","Drosophila_melanogaster","Rattus_norvegicus", "Caenorhabditis_elegans", "Saccharomyces_cerevisiae"]
	parser.add_argument("--organisms", nargs='+', default = ["Homo_sapiens", "Mus_musculus", "Drosophila_melanogaster"],type = str, choices = potential_organisms, help = "All the organisms to build the network. Defaults to human, mouse and fruit fly")
	potential_orthologs = ['metazoa','opisthokonta', 'Metazoa', 'Opisthokonta']
	parser.add_argument("--orthologs", default = "metazoa", type = str.lower, choices = potential_orthologs, help = "Defines the Eggnogg dataset for the gene orthology search. Defaults to metazoa, but if yeast or other fungi are included, use opisthokonta.")	
	args = parser.parse_args()
	#
	#
	#
	column_set = ['#ID(s) interactor A', 'ID(s) interactor B', 'Alt. ID(s) interactor A', 'Alt. ID(s) interactor B', 'Taxid interactor A', 'Taxid interactor B', 'Confidence value(s)', 'Interaction type(s)']
	all_organisms = args.organisms
	#Keep only the organism tax IDs of loaded organisms
	potential_tax_ids =['9606', '10090','7227','10116', '6239']
	dict_tax_id = dict(zip(potential_organisms,potential_tax_ids))
	organism_tax_ids = list(map(dict_tax_id.get, all_organisms))
	#Keep only the IntAct identifiers of loaded organisms
	potential_organism_set = ['taxid:9606(human)|taxid:9606(Homo sapiens)',
				'taxid:10090(mouse)|taxid:10090(Mus musculus)',
				'taxid:7227(drome)|taxid:7227("Drosophila melanogaster (Fruit fly)")',
				'taxid:10116(rat)|taxid:10116("Rattus norvegicus (Rat)")',
				'taxid:6239(caeel)|taxid:6239(Caenorhabditis elegans)']
	dict_org_set = dict(zip(potential_organisms, potential_organism_set))
	organism_set = list(map(dict_org_set.get, all_organisms))
	orthologs_id = "33208" if args.orthologs == "metazoa" else "33154"
	#
	#Network creation
	#
	start = time.time()
	G = nx.Graph()
	print("Starting.")
	print("Loading edge data")
	#Load edge data
	df_intact = pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/Intact/intact.txt")
	print("Proccesing edge data")
	#Reduce search space
	df_intact = df_intact.loc[df_intact['Taxid interactor A'].isin(organism_set) & df_intact['Taxid interactor B'].isin(organism_set), column_set]

	#Keep only intraspecies interactions
	df_intact = df_intact.loc[df_intact['Taxid interactor A'] == df_intact['Taxid interactor B']]

	#Get scores as numbers
	df_intact['weight'] = df_intact['Confidence value(s)'].str.extract('miscore:(.*)', expand = True)
	df_intact['weight'] = df_intact['weight'].astype('float')

	#Take the average of all scores for a certain pair of interactions
	df_edges = df_intact.groupby(['#ID(s) interactor A','ID(s) interactor B','Taxid interactor B'])['weight'].mean().reset_index()
	df_edges = df_edges.rename(columns = {'#ID(s) interactor A':'Gene_A','ID(s) interactor B':'Gene_B', 'Taxid interactor B':"organism"})

	#Get a list of all genes
	df_id = pd.concat([df_intact.loc[:, ['#ID(s) interactor A','Alt. ID(s) interactor A', 'Taxid interactor A']].rename(columns = {'#ID(s) interactor A':'gene_id','Alt. ID(s) interactor A':'alt_id','Taxid interactor A':'organism'}),
						df_intact.loc[:, ['ID(s) interactor B','Alt. ID(s) interactor B', 'Taxid interactor B']].rename(columns = {'ID(s) interactor B':'gene_id','Alt. ID(s) interactor B':'alt_id','Taxid interactor B':'organism'}) ])
	#Drop duplicates
	df_id = df_id.drop_duplicates()

	#Get IDs that have ensembl ID as a main ID, and those that have ensembl ID as their alternative ID
	df_id_1 = df_id.loc[df_id['gene_id'].str.contains('ensembl')]
	df_id_2 = df_id.loc[df_id['alt_id'].str.contains('ensembl')]

	#Extract the ENSEMBL ID from the string
	df_id_1['ENS_id'] = df_id_1['gene_id'].str.extract('ensembl:(.*)', expand = True)
	df_id_2['ENS_id'] = df_id_2['alt_id'].str.extract('ensembl:(.*)', expand = True)

	#Remove the IDs that come after the ENSEMBL ID
	df_id_2['ENS_id'] = df_id_2['ENS_id'].str.extract('^([^|]*).*', expand = True)
	
	#Remove the ID versions
	df_id_1['ENS_id'] = df_id_1['ENS_id'].str.extract('(^[^\\.]*)', expand = True)
	df_id_2['ENS_id'] = df_id_2['ENS_id'].str.extract('^([^.]*).*', expand = True)

	#Combine all of the IDs
	df_ENS_id = pd.concat([df_id_1.dropna(),df_id_2.dropna()])

	#Add the organism to each interaction
	dict_org_id = dict(zip(organism_set, all_organisms))

	df_ENS_id['organism'] = df_ENS_id['organism'].map(dict_org_id)
	df_ENS_id = df_ENS_id.drop_duplicates(subset = 'gene_id')

	#Get all edges into their ENSEMBL IDs
	df_edges['Gene_A'] = df_edges['Gene_A'].map(df_ENS_id.set_index('gene_id')['ENS_id'])
	df_edges['Gene_B'] = df_edges['Gene_B'].map(df_ENS_id.set_index('gene_id')['ENS_id'])
	
	#Remove the NaNs
	df_edges = df_edges.dropna()
	#Now the file has both ENSEMBL gene, protein and transcript IDs which all have to be mapped to Gene IDs

	def id_getter(organism):
		#Fetch the ENSEMBL ID vs UniProt ID files
		if organism in ["Drosophila_melanogaster","Caenorhabditis_elegans"]:
			#Drosophila and Caenorhabditis done hasve protein ID versions
			df_id = pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/IDs/Intact/{org}_ENSEMBL_ID.txt".format(org = organism))
			df_id = df_id.rename(columns = {'Gene stable ID' : 'Preferred_Name', 'Protein stable ID' : 'Protein_id', 'Transcript stable ID': 'Transcript_id'})
			df_id['organism'] = organism
		else:
			#The rest have them, and they are used in Eggnogg
			df_id = pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/IDs/Intact/{org}_ENSEMBL_ID.txt".format(org = organism)).drop('Protein stable ID version', axis = 1)
			df_id = df_id.rename(columns = {'Gene stable ID' : 'Preferred_Name', 'Protein stable ID' : 'Protein_id', 'Transcript stable ID': 'Transcript_id'})
			df_id['organism'] = organism
		return df_id

	#Combine all ID's into a single file
	df_all_id = pd.concat([id_getter(i) for i in all_organisms])

	#Separate out the transcript and protein IDs
	df_trans_id = df_all_id.loc[:, ["Preferred_Name", "Transcript_id"]].dropna().rename(columns = {"Transcript_id":"alt_id"})
	df_prot_id = df_all_id.loc[:,["Preferred_Name", "Protein_id"]].dropna().rename(columns = {"Protein_id":"alt_id"})

	#Combine all alternative ids (transcript and protein) into a single file
	df_dict_id = pd.concat([df_trans_id,df_prot_id]).drop_duplicates(subset = 'alt_id')

	#Map the edge data into ENSEMBL IDs 
	df_edges['Gene_A'] = df_edges['Gene_A'].map(df_dict_id.set_index('alt_id')['Preferred_Name'])
	df_edges['Gene_B'] = df_edges['Gene_B'].map(df_dict_id.set_index('alt_id')['Preferred_Name'])

	#Map the organism name
	df_edges['organism'] = df_edges['organism'].map(dict_org_id)




	#Get all nodes into one column
	df_all_nodes = pd.concat([df_edges.drop(['Gene_A', 'weight'], axis = 1).rename(columns = {'Gene_B' : 'node'}), df_edges.drop(['Gene_B', 'weight'], axis = 1).rename(columns = {'Gene_A' : 'node'})])
	#Drop duplicates and add the organism id
	df_all_nodes = df_all_nodes.drop_duplicates().rename(columns = {'organism' : 'layer'})
	#Add the nodes to the graph
	G.add_nodes_from(df_all_nodes['node'])
	#Add the node attributes
	attrs = df_all_nodes.set_index('node').T.to_dict()
	nx.set_node_attributes(G, attrs)
	print("Added nodes")
	print("Adding edges")


	#Drop organism column
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

	all_graph_nodes = df_all_nodes.rename(columns = {'layout' : 'organism'})
	

	df_prot_id = df_all_id.loc[:, ['Preferred_Name','Protein_id','organism']].dropna()
	#Keep the IDs that are in out nodes
	df_prot_id = df_prot_id.loc[df_prot_id['Preferred_Name'].isin(all_graph_nodes['node'])]

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
	
	#Map the gene to the organism
	df_all_ort_edges['Gene_A_c']=df_all_ort_edges['Gene_A_c'].map(df_all_id.set_index('alias')['organism'])
	df_all_ort_edges['Gene_B_c']=df_all_ort_edges['Gene_B_c'].map(df_all_id.set_index('alias')['organism'])
	
	#Filter out the interactons where both genes belong to the same organism
	df_all_ort_edges = df_all_ort_edges.loc[df_all_ort_edges['Gene_A_c'] != df_all_ort_edges['Gene_B_c']]
	
	#Turn it back into a tuple list for easy add to the graph
	all_ort_edges = list(zip(df_all_ort_edges['Gene_A'],df_all_ort_edges['Gene_B']))
	
	#
	#Exporting
	#
	print("Adding the cross-edges")
	#Add the edges to our graph
	print(G)
	G.add_edges_from(all_ort_edges, type='cross')
	print("Added all the cross-edges")
	print("Full network complete")
	print(G)


	print('Exporting the full newtwork')
	#Write in gpickle
	path_gpickle = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/IntAct/genome_network/Genome_multilayer_network.gpickle'
	ensurePathExists(path_gpickle)
	nx.write_gpickle(G,path_gpickle)
	print("Exported the .gpickle")

	path_edgelist = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/IntAct/genome_network/Genome_multilayer_network.edgelist'
	ensurePathExists(path_gpickle)
	nx.write_edgelist(G, path_edgelist)
	print("Exported the .edgelist")

	path_graphml = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/IntAct/genome_network/Genome_multilayer_network.graphml'
	ensurePathExists(path_graphml)
	nx.write_graphml(G, path_graphml)
	print("Exported the .graphml")

	print("Done with exporting")
	end = time.time()
	print("Done, runtime : {:}".format(end-start))
