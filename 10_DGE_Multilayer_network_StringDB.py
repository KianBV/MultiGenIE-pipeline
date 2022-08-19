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
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


#The cuttoff for the gene being expressed
TPM_cutoff = 1
#Potential organsims are "Homo_sapiens","Mus_musculus","Drosophila_melanogaster","Rattus_norvegicus","Danio rerio", "Caenorhabditis_elegans"
all_organisms = ["Homo_sapiens","Mus_musculus","Drosophila_melanogaster"]
#Potential IDs are '9606', '10090','7227', '10116', '7955', '6239'
organism_Egg_ids = ['9606', '10090', '7227']
#Potential cells are 'enterocyte', 'neuron', 'muscle', 'spermatogonia', 'spermatocyte', 'spermatid'
cell_type = 'spermatid'
database = 'StringDB'

start = time.time()
def global_network_maker(organism):

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
	df_dge = pd.read_csv("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/Test_data/{c}/{org}-FPKM-{c}.csv.gz".format(org = organism, c = cell_type))
	df_dge = df_dge.loc[df_dge['TPM']>TPM_cutoff]

	#Make the StringDB to ENS identifier
	dict_id = dict(zip(df_id['#string_protein_id'], df_id['alias']))

	#Remap the IDs
	df_edges['protein1'] = df_edges['protein1'].map(dict_id)
	df_edges['protein2'] = df_edges['protein2'].map(dict_id)
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
	df_id = df_id.loc[df_id['alias'].isin(df_dge['id_gene'])]
	#Return the ID
	return df_id




G = nx.Graph()
df_all_id = pd.concat([global_network_maker(i) for i in all_organisms])
dict_gene_prot = dict(zip(df_all_id['#string_protein_id'], df_all_id['alias']))
Egg_id_set = set(df_all_id['#string_protein_id'])


print("Done with IDs")
print("Starting with EGGnogg members file")
df_Egg = open_undefined_last_column_files(
    "C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/EGGnogg_files/33208_members.tsv.gz",
    n_fixed_cols=5,
    names = ['family', 'id_eggnog', '_1', '_2', 'aliases', 'species'],
    nrows = None
    )
#For keeping the columns and species we need
df_Egg = df_Egg.set_index('id_eggnog')['aliases']
wanted_species = frozenset(organism_Egg_ids)


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
#Make an aditional column to map the gene to their organism
df_all_ort_edges['Gene_A_c']=df_all_ort_edges['Gene_A']
df_all_ort_edges['Gene_B_c']=df_all_ort_edges['Gene_B']
#Map the gene to the organism
dict_filter = dict(zip(df_all_id['alias'],df_all_id['organism']))
df_all_ort_edges['Gene_A_c']=df_all_ort_edges['Gene_A_c'].map(dict_filter)
df_all_ort_edges['Gene_B_c']=df_all_ort_edges['Gene_B_c'].map(dict_filter)
#Filter out the interactons where both genes belong to the same organism
df_all_ort_edges = df_all_ort_edges.loc[df_all_ort_edges['Gene_A_c'] != df_all_ort_edges['Gene_B_c']]
#Turn it back into a tuple list for easy add to the graph
all_ort_edges = list(zip(df_all_ort_edges['Gene_A'],df_all_ort_edges['Gene_B']))




print("Adding the cross-edges")
#Add the edges to our graph
G.add_edges_from(all_ort_edges, type='cross', weight = 1)
print("Added all the cross-edges")
print("Full network complete")
print(G)

print('Exporting the full newtwork')
#Write in gpickle
path_gpickle = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/DGE_network/DGE_{c}_multilayer_network_TPM_cutoff_{tpm}.gpickle'.format(db = database, tpm = TPM_cutoff, c = cell_type)
ensurePathExists(path_gpickle)
nx.write_gpickle(G,path_gpickle)
print("Exported the .gpickle")
"""
path_edgelist = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/DGE_network/DGE_{c}_multilayer_network_TPM_cutoff_{tpm}.edgelist'.format(db = database, tpm = TPM_cutoff, c = cell_type)
ensurePathExists(path_gpickle)
nx.write_edgelist(G, path_edgelist)
print("Exported the .edgelist")

path_graphml = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/Gephi/graphml/DGE_{c}_multilayer_network_TPM_cutoff_{tpm}.graphml'.format(db = database, tpm = TPM_cutoff, c = cell_type)
ensurePathExists(path_graphml)
nx.write_graphml(G, path_graphml)
print("Exported the .graphml")
"""
print("Done with exporting")
end = time.time()
print("Done, runtime : {:}".format(end-start))