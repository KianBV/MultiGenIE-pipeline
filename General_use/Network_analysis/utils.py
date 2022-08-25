import os
import gzip
from io import StringIO
import pandas as pd
import networkx as nx
from collections import defaultdict, Counter


def transpose_variable_across_layers(G, variable, combination='sum'):
    """ Method to transpose results in one layer to another. Note duplicate results are either summed or set to majority."""
    dict_i_values = {i: d[variable] for i, d in G.nodes(data=True) if d.get(variable, None) is not None}
    dict_j_values = defaultdict(list)
    for i, v in dict_i_values.items():
        cross_edges = [j for _i, j, d in G.edges(i, data=True) if d.get('type', None) == 'cross']
        for j in cross_edges:
            dict_j_values[j].append(v)
    # Combine multiple values
    if combination == 'sum':
        dict_j_values = {k: sum(l) for k, l in dict_j_values.items()}
    elif combination == 'majority':
        dict_j_values = {k: Counter(l).most_common()[0][0] for k, l in dict_j_values.items()}
    else:
        TypeError("Combination must be either 'sum', or 'majority'.")
    # Set attributes to network
    nx.set_node_attributes(G, values=dict_j_values, name=variable)
    return G


def get_network_layer(G, layer=''):
    return G.subgraph([n for n, d in G.nodes(data=True) if (d.get('layer') == layer)]).copy()


def get_network_by_attribute(G, attribute='', value=''):
    return G.subgraph([n for n, d in G.nodes(data=True) if (d.get(attribute) == value)]).copy()


def get_network_largest_connected_component(G):
    largest_cc = max(nx.connected_components(G), key=len)
    return G.subgraph(largest_cc).copy()


def open_undefined_last_column_files(filepath, skiprows=0, n_fixed_cols=None, sep='\t', *args, **kwargs):
    """ Some StringDB files need manual parsing to be loaded as a pandas DataFrame."""
    with gzip.open(filepath, 'rt') as f:
        ios = u''
        # Skip header
        for i in range(skiprows):
            _ = f.readline()
        # Loop file content
        for i, line in enumerate(f, start=0):
            sline = line.split(sep)
            ios += u'\t'.join(sline[:n_fixed_cols]) + u'\t' + sline[-1]

        return pd.read_csv(StringIO(ios), sep='\t', encoding='utf-8', *args, **kwargs)


def ensurePathExists(path):
    dirname = os.path.dirname(path)
    if not os.path.exists(dirname):
        os.makedirs(dirname)