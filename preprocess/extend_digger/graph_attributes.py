# domain_G, DDI, gid2name
import pickle
import networkx as nx


def load_obj(organism, name):
    with open(f'../../container/domain/data/{organism}/{name}.pkl', 'rb') as f:
        return pickle.load(f)


def save_obj(organism, name, obj):
    with open(f'../../container/domain/data/{organism}/{name}_ext.pkl', 'wb') as f:
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


def add_predicted_nodes(graph, predicted_ddis_path, graph_type):
    with open(predicted_ddis_path, 'r') as f:
        for line in f.readlines():
            line = line.strip().split("\t")
            if not graph.has_edge(line[0], line[1]):
                graph.add_edge(line[0], line[1], confidence=line[2])
    return graph


def main(organism):
    graphs = ['DomainG', 'DDI']
    for graph in graphs:
        # load original graph
        try:
            original_graph = load_obj(organism, graph)
        except FileNotFoundError:
            continue
        annotated_graph = dummy_attribute(original_graph)
        # add predicted nodes
        extended_graph = add_predicted_nodes(annotated_graph,
                                             f'../resultdata/predicted_ddi_ppi.tsv', graph)

        save_obj(organism, graph, extended_graph)


if __name__ == '__main__':
    main('Homo Sapiens[human]')

    # domain_g_original: nx.Graph = load_obj('DomainG_original')
    # predicted_edges, predicted_edges_ddi = read_ddi_tsv('../domain/data/Homo sapiens[human]/predicted_ddi_ppi.tsv')
    # domain_g = new_graph(predicted_edges)
    # domain_g = annotate_graph(domain_g_original, domain_g)
    # save_obj("Domain_G_ext", domain_g)

