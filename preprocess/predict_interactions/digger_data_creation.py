# Files I need to update:
# DDI.pkl
# DomainG.pkl
# gid2name.pkl

import os
import pickle
import mygene
import pandas as pd
import networkx as nx


# load predicted interactions
def read_interactions(file: str, third_col=False):
    interactions = []
    with open(file, 'r') as f:
        f.readline()
        for line in f.readlines():
            line = line.strip().split("\t")
            if third_col:
                assoc = (line[0], line[1], line[2]) if line[0] < line[1] else (line[1], line[0], line[2])
            else:
                assoc = (line[0], line[1]) if line[0] < line[1] else (line[1], line[0])
            interactions.append(assoc)
    return interactions


# load ddi/ppi graph and add new domains to DomainG.pkl
def load_graph(file_path: str, ppi_ddi=False):
    with open(file_path, 'rb') as g:
        ddi_G = pickle.load(g)
        print(ddi_G)
    if ppi_ddi:
        return set((e[0], e[1]) if e[0].split("/")[1] < e[1].split("/")[1] else (e[1], e[0]) for e in ddi_G.edges)
    return set((e[0], e[1]) if e[0] < e[1] else (e[1], e[0]) for e in ddi_G.edges)


def combine_domaing_graph(graph: nx.Graph, interactions: set):
    confidence_map = {'Gold': 'high', 'Silver': 'mid', 'Bronze': 'low'}
    # add "original" attribute to every edge in the graph
    for edge in graph.edges:
        graph[edge[0]][edge[1]]['confidence'] = 'original'
    # add new interactions to the graph with respective confidence
    for edge in interactions:
        if not graph.has_edge(edge[0], edge[1]):
            graph.add_edge(edge[0], edge[1], confidence=confidence_map[edge[2]])
    return graph


def extend_domaing(path_graph, path_predicted):
    original_domain_g = load_graph(path_graph, ppi_ddi=True)
    original_graph = nx.Graph(original_domain_g)

    ddi_predicted = set(read_interactions(path_predicted, third_col=True))
    predicted_graph = combine_domaing_graph(original_graph, ddi_predicted)
    print("Predicted PPI/DDI interactions:", len(set(ddi_predicted)))
    print("Combined graph:", predicted_graph)
    # check how many edges are gold, silver, bronze
    gold = 0
    silver = 0
    bronze = 0
    for edge in predicted_graph.edges:
        if predicted_graph[edge[0]][edge[1]]['confidence'] == 'high':
            gold += 1
        elif predicted_graph[edge[0]][edge[1]]['confidence'] == 'mid':
            silver += 1
        elif predicted_graph[edge[0]][edge[1]]['confidence'] == 'low':
            bronze += 1
    print(f"Gold: {gold}\nSilver: {silver}\nBronze: {bronze}")
    # save as DomainG_ext.pkl
    return predicted_graph


def unknown_gene_ids(known_ids: set, graph: nx.Graph):
    unknown = set()
    for node in graph.nodes:
        gene_id = node.split("/")[0]
        if gene_id not in known_ids:
            unknown.add(gene_id)
    return unknown


def remove_unknown_nodes(graph: nx.Graph, known: set):
    nodes_to_remove = []
    for node in graph.nodes:
        gene_id = node.split("/")[0]
        if gene_id in known:
            continue
        nodes_to_remove.append(node)
    graph.remove_nodes_from(nodes_to_remove)
    return graph


def update_gid2name(gid2name_path: str, mart_table: str, graph: nx.Graph):
    mg = mygene.MyGeneInfo()
    mart_table = pd.read_csv(mart_table, sep='\t')
    # convert entrez id to int
    mart_table = mart_table.astype({'NCBI gene (formerly Entrezgene) ID': str})
    mart_table['NCBI gene (formerly Entrezgene) ID'] = mart_table['NCBI gene (formerly Entrezgene) ID'].str.split('.').str[0]
    gid_map = mart_table.set_index('NCBI gene (formerly Entrezgene) ID')['Gene name'].to_dict()
    gid2name = pickle.load(open(gid2name_path, 'rb'))
    new_genes = unknown_gene_ids(set(gid_map.keys()), graph)
    print(f"New genes: {len(new_genes)}")
    online_query = []
    for gene in new_genes:
        if gene in gid_map and gene not in gid2name:
            gid2name[gene] = gid_map[gene]
        else:
            online_query.append(gene)
    print(f"Querying online: {len(online_query)}")
    # query online
    out = mg.querymany(online_query, scopes='entrezgene', fields='symbol', species='human')
    for entry in out:
        if 'symbol' in entry:
            gid2name[entry['query']] = entry['symbol']
    print(f"gid2name length: {len(gid2name)}")
    return gid2name


def extend_data(input_path: str, output_path: str, mart_table: str, predicted_interactions: str):
    graph = extend_domaing(f'{input_path}/DomainG.pkl', predicted_interactions)
    gid2name = update_gid2name(f'{input_path}/gid2name.pkl', mart_table, graph)

    domain_g = remove_unknown_nodes(graph, set(gid2name.keys()))
    print(f"Graph after removing unknown nodes: {domain_g}")
    pickle.dump(domain_g, open(f"{output_path}/DomainG_ext.pkl", 'wb'))
    pickle.dump(gid2name, open(f'{output_path}/gid2name_ext.pkl', 'wb'))


if __name__ == '__main__':
    extend_data('../data_digger', '../data_digger',
                '../sourcedata/mart_export.txt', '../predicted_ddi_ppi.tsv.tsv')
    print("Done, exiting.")
