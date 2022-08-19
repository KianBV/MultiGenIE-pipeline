#Date: 8.8.2022.
#Written by: Kian Bigović Villi
#The script inputs a geneMANIA edgelist of a set of organismns, changes all of the Gene IDs to ENSEMBL gene ID, and outputs a multilayered network of the GeneMANIA interactions. It can also subgraph based upon input gene expression data.
#The full script takes around 12ish minutes on my computer (16GB RAM memory, 2,5 GHz CPU, 512 GB SSD)
import numpy as np
import pandas as pd
import networkx as nx
from utils import ensurePathExists, open_undefined_last_column_files
import time 
from itertools import chain, product

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
database = "GeneMANIA"
#Potential organsims are "Homo_sapiens","Mus_musculus","Drosophila_melanogaster","Rattus_norvegicus","Danio_rerio","Caenorhabditis_elegans"
all_organisms = ["Homo_sapiens","Mus_musculus","Drosophila_melanogaster","Rattus_norvegicus","Danio_rerio","Caenorhabditis_elegans"]
organism_Egg_ids = ['9606', '10090','7227', '10116', '7955', '6239']


def global_network_maker(organism):
	#Creates a network based from the organism inputted. It requires the GeneMANIA "identifier_mappings.txt" and "COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt" file for the appropriate organism
	#Inputs the organism name in a "Genus_description" format eg. "Homo_sapiens", adds nodes and edges to a global network variable "G"

	global G, dict_ENS_id
	#Creating the whole full genome network
	print("Creating the full network for {:}".format(organism))
	#Input the graph edge data and rename the weight variable (case sensitivity)
	df_edges=pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/{db}/whole_networks/{x}_whole_network.txt".format(db = database, x = organism)).rename({"Weight" : "weight"}, axis = 1)
	#Chaging the ID to the ENS ID
	df_edges["Gene_A"] = df_edges["Gene_A"].map(dict_ENS_id)
	df_edges["Gene_B"] = df_edges["Gene_B"].map(dict_ENS_id)
	#Gets the list of of the nodes in the network
	print("Adding nodes")
	all_nodes = pd.unique(df_edges[['Gene_A', 'Gene_B']].values.ravel( 'K'))
	#Adds the nodes to the graph
	G.add_nodes_from(all_nodes, layer = organism)
	print("Added nodes, starting with edges")
	#Add edge type 
	df_edges["type"]='intra'
	#Index the edges
	df_edges = df_edges.set_index(['Gene_A','Gene_B'])
	edge_idx = df_edges.index.to_list()
	#Put the attributes inside a dictionary
	edge_attr = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_edges.to_dict('records')]
	#combine the edges 
	G.add_edges_from([(i, j, d) for (i, j), d in zip(edge_idx, edge_attr)])
	print(G)
start = time.time()
G = nx.Graph()



def id_getter(organism):
	#combines all the id files into a single file and indexes the organism
	df_id = pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/{db}/ID_files/raw/{x}_id.txt".format(db = database, x = organism))
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


for i in all_organisms:
	global_network_maker(i)
print("Add organisms added, starting with cross edges")
#Create a df with all nodes and their respecitve organism from the graph data
all_graph_nodes = pd.DataFrame(G.nodes(data = 'layer'), columns = ['node', 'organism'])


#Make a df that links ENSEBL Protein and Gene id
df_prot_id = df_id.loc[(df_id["Source"]=="Ensembl Protein ID")].rename( {"Name":"Protein_id"}, axis = 1)
#Combine all of the ID pairs in a single dataframe and keep only those that are in our network
df_prot_id = df_prot_id.loc[df_prot_id['Preferred_Name'].isin(all_graph_nodes['node'])]	
#Create the organism, EGGnog id dictionary
dict_organism = dict(zip(all_organisms, organism_Egg_ids))
#Map the EGGnog id
df_prot_id['organism'] = df_prot_id['organism'].map(dict_organism)
#Create the EGGnogg ID : ENSEMBL gene ID dictionary for mapping later
dict_gene_prot = dict(zip(df_prot_id['organism'] + '.' + df_prot_id['Protein_id'], df_prot_id['Preferred_Name']))
#Create a set of all of the EGGnogg IDs in our networks
Egg_id_set = set(df_prot_id['organism'] + '.' + df_prot_id['Protein_id'])

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
dict_filter = dict(zip(df_id['Preferred_Name'],df_id['organism']))
df_all_ort_edges['Gene_A_c']=df_all_ort_edges['Gene_A_c'].map(dict_filter)
df_all_ort_edges['Gene_B_c']=df_all_ort_edges['Gene_B_c'].map(dict_filter)
#Filter out the interactons where both genes belong to the same organism
df_all_ort_edges = df_all_ort_edges.loc[df_all_ort_edges['Gene_A_c'] != df_all_ort_edges['Gene_B_c']]
#Turn it back into a tuple list for easy add to the graph
all_ort_edges = list(zip(df_all_ort_edges['Gene_A'],df_all_ort_edges['Gene_B']))



print("Adding the cross-edges")
#Add the edges to our graph
G.add_edges_from(all_ort_edges, type='cross')
print("Added all the cross-edges")
print("Full network complete")
print(G)

print('Exporting the full newtwork')
#Write in gpickle
path_gpickle = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/Total_network/Genome_multilayer_network.gpickle'.format(db = database)
ensurePathExists(path_gpickle)
nx.write_gpickle(G,path_gpickle)
path_edgelist = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/Total_network/Genome_multilayer_network.edgelist'.format(db = database)
ensurePathExists(path_gpickle)
nx.write_edgelist(G, path_edgelist)

end = time.time()
print("Done, runtime : {:}".format(end-start))


#The DGE code, input the DGEs and makes a subgraph based on them.
#Note that it will most definately be faster to make a smaller network from the start than to make a big one and filter it - See 01_DGE_Multilayer_network.py



'''#Create a dictionary to map all possible gene IDs to ENSEMBL Gene ID
dict_dge = dict(zip(df_id['Name'],df_id['Preferred_Name']))


def dge_getter(organism):
	#Get the DGE data for all organisms, convert the ID and combine it into a single dataframe

	df_dge = pd.read_csv("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/Test_data/spermatocyte_FPKM/{:}-FPKM-spermatocyte.csv.gz".format(organism))
	df_dge["id_gene"] = df_dge["id_gene"].map(dict_dge)
	df_dge = df_dge.dropna(subset=["id_gene"])
	return df_dge

all_dge = pd.concat([dge_getter(i) for i in dge_organisms])
#Filter out the DGEs, this can be done beforehand (then the step would be ignored)
all_dge = all_dge.loc[all_dge['TPM']>2,:]
#Make a subgraph out of the selected DGEs
SG = G.subgraph(all_dge['id_gene'])
print(SG) 

#Exporting the subgraph
print('Exporting the dge newtwork')
#Write in gpickle
path_gpickle = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/DGE_network/DGE_global_network.gpickle'.format(db = database)
ensurePathExists(path_gpickle)
nx.write_gpickle(SG,path_gpickle)
path_edgelist = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/DGE_network/DGE_global_network.edgelist'.format(db = database)
ensurePathExists(path_gpickle)
nx.write_edgelist(SG, path_edgelist)

end = time.time()
print("Done, runtime : {:}".format(end-start))'''
