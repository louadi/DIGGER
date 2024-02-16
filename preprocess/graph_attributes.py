# domain_G, DDI, gid2name
import pickle
import networkx as nx

def load_obj(name):
    with open('../domain/data/Mus musculus[mouse]/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

def save_obj(name, obj):
    with open('../domain/data/Mus musculus[mouse]/' + name + '.pkl', 'wb') as f:
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


def read_interactions():
    all_interactions = []
    for i in ['gold', 'silver', 'bronze']:
        with open(f'/mnt/d/programming/bachelor_projects/ppidm-check/resultdata/interactions_{i}', 'r') as f:
            print(f"Interactions for {i} confidence")
            for line in f.readlines():
                # replace the first letter in i with the uppercase version
                all_interactions.append(line.strip().split("\t") + [i[0].upper() + i[1:]])
    return all_interactions


def dummy_attribute(graph):
    # add 'original' attribute to all edges
    for edge in graph.edges(data=True):
        graph[edge[0]][edge[1]]['confidence'] = 'original'
    return graph


if __name__ == '__main__':
    graphs = ['DomainG', 'DDI']
    for graph in graphs:
        # load original graph
        original_graph = load_obj(graph)
        annotated_graph = dummy_attribute(original_graph)
        save_obj(graph + '_annotated', original_graph)

    # domain_g_original: nx.Graph = load_obj('DomainG_original')
    # predicted_edges, predicted_edges_ddi = read_ddi_tsv('../domain/data/Homo sapiens[human]/predicted_ddi_ppi.tsv')
    # domain_g = new_graph(predicted_edges)
    # domain_g = annotate_graph(domain_g_original, domain_g)
    # save_obj("Domain_G_ext", domain_g)

