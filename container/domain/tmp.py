import os.path
import pickle
import re
import networkx as nx
import mygene

def load_obj(name):
    with open('data/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


def vis_node_(node, DomainG):
    DomainG = DomainG
    G = nx.Graph()
    if DomainG.has_node(node):
        # copy over the edges with their attributes
        for edge in DomainG.edges(node, data=True):
            G.add_edge(edge[0], edge[1], **edge[2])

    # get the amount of edges that are predicted and how many are original
    predicted = 0
    original = 0
    for edge in G.edges(data=True):
        if edge[2]['origin'] == 'predicted':
            predicted += 1
        else:
            original += 1
    print(f"Node {node} has {predicted} predicted edges and {original} original edges")

    g = nx.Graph()
    g.add_node(node)
    for n in G.nodes():
        domain = n.split("/")[1]
        if n != node:
            # add predicted attribute if it is predicted in G
            if (G.has_edge(n, node) or G.has_edge(node, n)) \
                    and (G[node][n]['origin'] == 'predicted' or G[n][node]['origin'] == 'predicted'):
                if not g.has_edge(n, domain):
                    g.add_edge(n, domain, origin='predicted')
                if not g.has_edge(node, domain):
                    g.add_edge(node, domain, origin='predicted')
            else:
                g.add_edge(n, domain, origin='original')
                g.add_edge(node, domain, origin='original')

    # get the amount of edges that are predicted and how many are original
    predicted = 0
    original = 0
    for edge in g.edges(data=True):
        if edge[2]['origin'] == 'predicted':
            predicted += 1
        else:
            original += 1

    print(f"Node {node} has {predicted} predicted edges and {original} original edges")


def graph_test(DomainG, node):
    G = nx.Graph()
    edge_attr = set()
    for edge in DomainG.edges(data=True):
        edge_attr.add(edge[2]['confidence'])
    print(edge_attr)
    if DomainG.has_node(node):
        G.add_edges_from(DomainG.edges(node, data=True))

        # copy over the edges with their attributes
        # for edge in DomainG.edges(node, data=True):
        #     G.add_edge(edge[0], edge[1], **edge[2])
    # check if edge data is there
    print(G.edges(data=True))
    return G

    g = nx.Graph()
    g.add_node(node)
    for n in G.nodes():
        domain = n.split("/")[1]
        if n == node:
            continue
        # add predicted attribute if it is predicted in G
        if G.has_edge(n, node) and (G[node][n]['confidence'] != 'original' or G[n][node]['confidence'] != 'original'):
            print(G[node][n]['confidence'], G[n][node]['confidence'])
            if not g.has_edge(n, domain):
                g.add_edge(n, domain, confidence=G[n][node]['confidence'])
            if not g.has_edge(node, domain):
                g.add_edge(node, domain, confidence=G[n][node]['confidence'])

        else:
            g.add_edge(n, domain, confidence='original')
            g.add_edge(node, domain, confidence='original')


# returns true if a graph1 is a subgraph of graph2
def is_subgraph(graph1, graph2):
    for i in graph1.nodes():
        if i not in graph2.nodes():
            return False
    for i in graph1.edges():
        if i not in graph2.edges():
            return False
    return True


# get node and edge overlap between two graphs
def graph_overlap(graph1, graph2):
    node_overlap = 0
    for i in graph1.nodes():
        if i in graph2.nodes():
            node_overlap += 1
    edge_overlap = 0
    for i in graph1.edges():
        if i in graph2.edges():
            edge_overlap += 1
    return node_overlap, edge_overlap


def remove_nan_nodes(graph):
    for node in list(graph.nodes):  # Create a copy of the nodes
        if "nan" in node:
            graph.remove_node(node)
    return graph


def extract_confidence(domainG):
    confidences = set()
    for edge in domainG.edges(data=True):
        confidences.add(edge[2].get('confidence', 'original'))
    print(confidences)


def Ensemb_to_entrez(genes, organsim='human'):
    mg = mygene.MyGeneInfo()
    out = mg.querymany(genes, scopes="ensembl.gene", fields='entrezgene', species=organsim, verbose=False)
    translated = {}
    for result in out:
        if 'entrezgene' in result:
            translated[result['query']] = result['entrezgene']
    return translated


if __name__ == '__main__':
    pass
    # ppi_graph: nx.Graph = pickle.load(open('../domain/data/Homo sapiens[human]/DomainG.pkl', 'rb'))
    # print(len(ppi_graph.nodes))
    # ppi_graph = remove_nan_nodes(ppi_graph)
    # print(len(ppi_graph.nodes))
    # pickle.dump(ppi_graph, open('../domain/data/Homo sapiens[human]/DomainG_up.pkl', 'wb'))


    # edges_domainV = {'test': [1,2,3,4], 'original': [1,2,3,4]}
    # print(len(edges_domainV.get('original')) > 70 if edges_domainV.get('original') else True)
    # # load human DDI graph
    # DomainG_human = pickle.load(open('data/Homo sapiens[human]/DomainG.pkl', 'rb'))
    # extract_confidence(DomainG_human)
    # DomainG_mouse = pickle.load(open('data/Mus musculus[mouse]/DomainG.pkl', 'rb'))
    # print(len(DomainG_mouse.edges), len(DomainG_mouse.nodes))
    #
    # # check if graphs aver overlapping
    # print(is_subgraph(DomainG_mouse, DomainG_human))
    # node, edge = graph_overlap(DomainG_mouse, DomainG_human)
    # print(f"Edges in both: {edge}, Nodes in both: {node}")
    # #domaing = pickle.load(open('data/Homo sapiens[human]/DomainG.pkl', 'rb'))
    # #graph_test(domaing, '10114/PF00069')
