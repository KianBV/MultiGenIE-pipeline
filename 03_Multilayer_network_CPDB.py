import numpy as np
import pandas as pd
import networkx as nx
from utils import ensurePathExists, open_undefined_last_column_files
import time 
from itertools import chain, product
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


start = time.time()

TPM_cutoff = 60
#Potential organsims are "Homo_sapiens","Mus_musculus","Drosophila_melanogaster","Rattus_norvegicus","Caenorhabditis_elegans"

all_organisms = ["Homo_sapiens","Mus_musculus"]
organism_Egg_ids = ['9606', '10090']


df_id_hs = pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/IDs/Homo_sapiens_ENSEMBL_ID.txt")
print('Starting with the network')
G = nx.Graph()
print("Creating the full network for Homo sapiens")
#Loads the database
df_data_hs = pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/CPDB/ConsensusPathDB_Homo_sapiens_PPI.gz", skiprows = 1)
#Removes the interactions lacking a confidence score
df_data_hs = df_data_hs.dropna()
#Save the ENSEMBL ID
df_edges_hs = df_data_hs['interaction_participants__ensembl_gene'].str.split(',', 1, expand = True).rename(columns = {0:'Gene_A',1:'Gene_B'})
#Define the edge weight
df_edges_hs['weight'] = df_data_hs['interaction_confidence']
#Remove pairs that have missing values
df_edges_hs = df_edges_hs.replace(r'^\s*$', np.nan, regex=True).dropna()
#Keep only the ones that code for a protein
df_edges_hs = df_edges_hs.loc[df_edges_hs['Gene_A'].isin(df_id_hs['Gene stable ID']) & df_edges_hs['Gene_B'].isin(df_id_hs['Gene stable ID'])]
#Obtains the list of all nodes
all_nodes_hs = pd.Series(pd.unique(df_edges_hs[['Gene_A', 'Gene_B']].values.ravel('K')))
print("Adding nodes")
#Adds nodes
G.add_nodes_from(all_nodes_hs, layer = 'Homo_sapiens')
print("Added nodes, starting with edges")

#Add edge type 
df_edges_hs["type"]='intra'
#Index the edges
df_edges_hs = df_edges_hs.set_index(['Gene_A','Gene_B'])
edge_idx_hs = df_edges_hs.index.to_list()
#Put the attributes inside a dictionary
edge_attr_hs = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_edges_hs.to_dict('records')]
#combine the edges 
G.add_edges_from([(i, j, d) for (i, j), d in zip(edge_idx_hs, edge_attr_hs)])
print(G)
print("Homo sapines added, starting with cross edges")

#Mus Musculus data - a bit more complicated because it first has to be mapped to ENSEMBL ID 

#Import the ENSEMBL : UniProt ID file
df_id_mm = pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/IDs/Mus_musculus_ENSEMBL_ID.txt")
#Input the Uniprot Entry Name : Gene ID file
df_uni_id_mm = pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/IDs/mm_uniprot_id.gz")
#Input the CPDB data
df_data_mm = pd.read_table("C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/Raw_data/CPDB/ConsensusPathDB_Mus_musculus_PPI.gz", skiprows = 1)
#Remove missing values
df_data_mm = df_data_mm.dropna()


#Splits the interactions column into two
df_edges_mm = df_data_mm['interaction_participants'].str.split(',', 1, expand = True).rename(columns = {0:'Gene_A',1:'Gene_B'})
#Add the weights
df_edges_mm['weight'] = df_data_mm['interaction_confidence']


#Make the dictionary to transform UniProt EntryName to GeneID
dict_in_UNI = dict(zip(df_uni_id_mm['Entry Name'], df_uni_id_mm['Entry']))
#Map it to the edges
df_edges_mm["Gene_A"] = df_edges_mm["Gene_A"].map(dict_in_UNI)
df_edges_mm["Gene_B"] = df_edges_mm["Gene_B"].map(dict_in_UNI)
#Make ENSEMBL : UniProt ID dictionary to map it all back to ENSEMBL ID
dict_ENS_UNI = dict(set(zip(df_id_mm['UniProtKB Gene Name ID'],df_id_mm['Gene stable ID'])))
df_edges_mm["Gene_A"] = df_edges_mm["Gene_A"].map(dict_ENS_UNI)
df_edges_mm["Gene_B"] = df_edges_mm["Gene_B"].map(dict_ENS_UNI)
#Remove NaNs
df_edges_mm = df_edges_mm.dropna()

#Add Nodes
#Gets all mm nodes
all_nodes_mm = pd.unique(df_edges_mm[['Gene_A', 'Gene_B']].values.ravel('K'))
print("Adding nodes")
#Adds nodes
G.add_nodes_from(all_nodes_mm, layer = 'Mus_musculus')
print("Added nodes, starting with edges")

#Add edge type 
df_edges_mm["type"]='intra'
#Index the edges
df_edges_mm = df_edges_mm.set_index(['Gene_A','Gene_B'])
edge_idx_mm = df_edges_mm.index.to_list()
#Put the attributes inside a dictionary
edge_attr_mm = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_edges_mm.to_dict('records')]
#combine the edges 
G.add_edges_from([(i, j, d) for (i, j), d in zip(edge_idx_mm, edge_attr_mm)])
print(G)
print("Mus_musculus added, starting with cross edges")

print("All organisms added, starting with cross edges")
#Create a df with all nodes and their respecitve organism from the graph data
all_graph_nodes = pd.DataFrame(G.nodes(data = 'layer'), columns = ['node', 'organism'])


#Adds organism tags
df_id_hs['organism'] = 'Homo_sapiens'
df_id_mm['organism'] = 'Mus_musculus'
#Combine HS and MM IDs and keep only those with protein IDs
df_prot_id = pd.concat([df_id_hs, df_id_mm]).dropna(subset = ['Protein stable ID'])
#Add the protin IDs and their versions to the list of proteins
df_prot_id = pd.concat([df_prot_id.drop('Protein stable ID', axis = 1).rename(columns = {'Protein stable ID version' : 'Protein stable ID'}), df_prot_id.drop('Protein stable ID version', axis = 1)])
df_prot_id = df_prot_id.rename(columns = {'Gene stable ID':'Preferred_Name','Protein stable ID':'Protein_id'})
#Use only the ones in our graph
df_prot_id = df_prot_id.loc[df_prot_id['Preferred_Name'].isin(all_graph_nodes['node'])]	

#Create the organism, EGGnog id dictionary
dict_organism = dict(zip(all_organisms, organism_Egg_ids))
#Map the EGGnog id
df_prot_id['organism'] = df_prot_id['organism'].map(dict_organism)
#Create the EGGnogg ID : ENSEMBL gene ID dictionary for mapping later
dict_gene_prot = dict(zip(df_prot_id['organism'] + '.' + df_prot_id['Protein_id'], df_prot_id['Preferred_Name']))
#Create a set of all of the EGGnogg IDs in our networks
Egg_id_set = set(df_prot_id['organism'] + '.' + df_prot_id['Protein_id'])
print(Egg_id_set)



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
dict_filter = dict(zip(df_prot_id['Preferred_Name'],df_prot_id['organism']))
df_all_ort_edges['Gene_A_c']=df_all_ort_edges['Gene_A_c'].map(dict_filter)
df_all_ort_edges['Gene_B_c']=df_all_ort_edges['Gene_B_c'].map(dict_filter)
#Filter out the interactons where both genes belong to the same organism
df_all_ort_edges = df_all_ort_edges.loc[df_all_ort_edges['Gene_A_c'] != df_all_ort_edges['Gene_B_c']]
#Turn it back into a tuple list for easy add to the graph
all_ort_edges = list(zip(df_all_ort_edges['Gene_A'],df_all_ort_edges['Gene_B']))




print("Adding the cross-edges")
#Add the edges to our graph
print(G)
G.add_edges_from(all_ort_edges, type='cross')
print("Added all the cross-edges")
print("Full network complete")
print(G)
print(G)

print('Exporting the full newtwork')
#Write in gpickle
path_gpickle = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/CPDB/DGE_network/03_global_network.gpickle'
ensurePathExists(path_gpickle)
nx.write_gpickle(G,path_gpickle)
print("Exported the .gpickle")

path_edgelist = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/CPDB/DGE_network/03_global_network.edgelist'
ensurePathExists(path_gpickle)
nx.write_edgelist(G, path_edgelist)
print("Exported the .edgelist")

path_graphml = 'C:/Users/Kian/Desktop/Kian_Praksa/IGC/databases/results/CPDB/Gephi/graphml/03_global_network_.graphml'
ensurePathExists(path_graphml)
nx.write_graphml(G, path_graphml)
print("Exported the .graphml")

print("Done with exporting")
end = time.time()
print("Done, runtime : {:}".format(end-start))