# domain_G, DDI, gid2name
import os.path
import pickle
import networkx as nx


def load_obj(organism, name):
    with open(f'../container/domain/data/{organism}/{name}.pkl', 'rb') as f:
        return pickle.load(f)


def save_obj(organism, name, obj):
    path = f'../container/domain/data/{organism}/{name}.pkl'
    # don't overwrite existing backup files, just to be sure
    if "bak" in name:
        while os.path.isfile(path):
            try:
                num = int(path[-5])
            except ValueError:
                num = 1
            path = path.replace(".pkl", f"_{num}.pkl")
    with open(path, 'wb') as f:
        pickle.dump(obj, f)


def read_ddi_tsv(filename):
    interactions = []
    interactions_ddi = []
    with open(filename, 'r') as f:
        for line in f.readlines():
            interactions.append(line.strip().split("\t"))
            # for the DDI graph
            try:
                line = line.strip().split("\t")
                partner1 = line[0].split("/")[1]
                partner2 = line[1].split("/")[1]
                interactions_ddi.append((partner1, partner2, line[2]))
            except IndexError:
                continue
    return interactions[1:], interactions_ddi


def annotate_graph(orig_graph, pred_graph):
    confidence_map = {'Gold': 'high', 'Silver': 'mid', 'Bronze': 'low'}
    print(f"Original graph has {len(orig_graph.edges)} edges")
    print(f"Extended graph has {len(pred_graph.edges)} edges")
    new_graph = nx.Graph()
    # get DDI_original edges
    for edge in orig_graph.edges:
        new_graph.add_edge(edge[0], edge[1], confidence='original')

    added_edges = 0
    skipped = 0
    for edge in pred_graph.edges(data=True):
        # check if edge is already in DDI_annotated
        if new_graph.has_edge(edge[0], edge[1]) or new_graph.has_edge(edge[0], edge[1]):
            skipped += 1
            continue
        new_graph.add_edge(edge[0], edge[1], confidence=confidence_map[edge[2]['confidence']])
        added_edges += 1

    print("Added", added_edges, "edges")
    print("Skipped", skipped, "edges")
    print(f"New graph has {len(new_graph.edges)} edges\n")
    return new_graph


def new_graph(edge_list):
    G = nx.Graph()
    for edge in edge_list:
        G.add_edge(edge[0], edge[1], confidence=edge[2])
    print(G)
    return G


def dummy_attribute(graph):
    # add 'original' attribute to all edges
    for edge in graph.edges(data=True):
        graph[edge[0]][edge[1]]['confidence'] = 'original'
    return graph


def add_predicted_nodes(graph, predicted_ddis_path, ppi_graph):
    mapping = {'Gold': 'high', 'Silver': 'mid', 'Bronze': 'low'}
    with open(predicted_ddis_path, 'r') as f:
        # skip header
        f.readline()
        for line in f.readlines():
            line = line.strip().split("\t")
            prot_1 = line[0].split("/")[0]
            prot_2 = line[1].split("/")[0]
            # if there is no edge between the proteins, a DDI can't exist
            if (not graph.has_node(line[0]) or not graph.has_node(line[1])) and (not ppi_graph.has_edge(prot_1, prot_2)):
                continue
            if not graph.has_edge(line[0], line[1]):
                graph.add_edge(line[0], line[1], confidence=mapping[line[2]])
    return graph


def main(organism, backup=True):
    # load PPI graph
    ppi_graph = load_obj(organism, 'PPI')
    print(f"PPI graph has {len(ppi_graph.edges)} edges and {len(ppi_graph.nodes)} nodes")

    graphs = ['DomainG', 'DDI']
    for graph in graphs:
        # load original graph
        try:
            original_graph = load_obj(organism, graph)
        except FileNotFoundError:
            continue
        print(f"Original graph has {len(original_graph.edges)} edges and {len(original_graph.nodes)} nodes")
        annotated_graph = dummy_attribute(original_graph)
        if backup:
            save_obj(organism, f"{graph}.bak", annotated_graph)
        # add predicted nodes
        extended_graph = add_predicted_nodes(annotated_graph,
                                             f'resultdata/predicted_ddi_ppi_alt.tsv', ppi_graph)
        print(f"Extended {graph} graph has {len(extended_graph.edges)} edges and {len(extended_graph.nodes)} nodes")
        save_obj(organism, graph, extended_graph)


if __name__ == '__main__':
    main('Mus musculus[mouse]')

    # domain_g_original: nx.Graph = load_obj('DomainG_original')
    # predicted_edges, predicted_edges_ddi = read_ddi_tsv('../domain/data/Homo sapiens[human]/predicted_ddi_ppi_alt.tsv')
    # domain_g = new_graph(predicted_edges)
    # domain_g = annotate_graph(domain_g_original, domain_g)
    # save_obj("Domain_G_ext", domain_g)

