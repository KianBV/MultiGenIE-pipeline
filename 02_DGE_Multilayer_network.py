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
database = "GeneMANIA"
TPM_cutoff = 60
#Potential organsims are "Homo_sapiens","Mus_musculus","Drosophila_melanogaster","Rattus_norvegicus","Danio_rerio","Caenorhabditis_elegans"

all_organisms = ["Homo_sapiens","Mus_musculus","Drosophila_melanogaster"]

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


	#Input the DGE data
	df_dge = pd.read_csv("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/Test_data/spermatocyte_FPKM/{:}-FPKM-spermatocyte.csv.gz".format(organism))
	#Change the DGE IDs
	df_dge["id_gene"] = df_dge["id_gene"].map(dict_dge)
	#Remove the NAs and filter the dge data
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
#Make a df that links ENSEBL Protein and Gene id
df_prot_id = df_id.loc[(df_id["Source"]=="Ensembl Protein ID")].rename( {"Name":"Protein_id"}, axis = 1)
print("Done with ID files")
#Create a dict that will match all possible gene IDs to the GeneMANIA preffered ID
dict_dge = dict(zip(df_id['Name'],df_id['Preferred_Name']))


for i in all_organisms:
	global_network_maker(i)
print("Add organisms added, starting with cross edges")
#Create a df with all nodes and their respecitve organism from the graph data
all_graph_nodes = pd.DataFrame(G.nodes(data = 'layer'), columns = ['node', 'organism'])



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

print("Adding the cross-edges")
#Add the edges to our graph
G.add_edges_from(all_ort_edges, type='cross')
print("Added all the cross-edges")
print("Full network complete")
print(G)

print('Exporting the full newtwork')
#Write in gpickle
path_gpickle = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/DGE_network/02_DGE_global_network_TPM_cutoff_{tpm}.gpickle'.format(db = database, tpm = TPM_cutoff)
ensurePathExists(path_gpickle)
nx.write_gpickle(G,path_gpickle)
print("Exported the .gpickle")

path_edgelist = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/DGE_network/02_DGE_global_network_TPM_cutoff_{tpm}.edgelist'.format(db = database, tpm = TPM_cutoff)
ensurePathExists(path_gpickle)
nx.write_edgelist(G, path_edgelist)
print("Exported the .edgelist")

path_graphml = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/{db}/Gephi/02_DGE_global_network_TPM_cutoff_{tpm}.graphml'.format(db = database, tpm = TPM_cutoff)
ensurePathExists(path_graphml)
nx.write_graphml(G, path_graphml)
print("Exported the .graphml")

print("Done with exporting")
end = time.time()
print("Done, runtime : {:}".format(end-start))
