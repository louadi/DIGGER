from domain.Process.load_data import DomainG_all, PPI_all, g2d_all
import domain.Process.process_data as pr
import domain.Process.exonstodomain as exd


def create_subgraph(id_names: list, organism: str):
    """
    Create a subgraph containing the protein nodes as well as the domain nodes and each
    edge between them.
    """
    PPI = PPI_all[organism]
    g2d = g2d_all[organism]
    DomainG = DomainG_all[organism]

    g = exd.nx.Graph()
    # list contains protein confirmed by both PPI and DDI (for visualization)
    confirmed_proteins = []
    id_set = set([name[1] for name in id_names])
    missing_domains = set()
    for name in id_set:
        entrez_id = name
        if PPI.has_node(entrez_id):
            # filter PPI edges to only include ones where the destination node is in the list
            filtered_edges = [(source, target, data) for source, target, data in PPI.edges(entrez_id, data=True) if target in id_set]
            g.add_edges_from(filtered_edges)

        # search for the protein domains interactions:
        # before coming here need to check if the protein have known domains
        protein_domains = g2d[entrez_id]

        for domain in protein_domains:
            node = entrez_id + '/' + domain
            g.add_node(node)

            if domain not in [x[2] for x in id_names if x[1] == name][0]:
                missing_domains.add(node)
                continue

            # add domain interactions to the graph:
            if DomainG.has_node(node):
                filtered_edges = [(source, target, data) for source, target, data in DomainG.edges(node, data=True) if
                                  target.split("/")[0] in id_set]
                g.add_edges_from(filtered_edges)

        edges = []

        # link domains with their proteins
        for node in g.nodes():
            # a node is domain node:
            if len(node.split('/')) == 2:
                edges.append((node.split('/')[0], node))
                confirmed_proteins.append((node.split('/')[0]))

        g.add_edges_from(edges)
        confirmed_proteins = list(set(confirmed_proteins))
    return g, confirmed_proteins, missing_domains


def vis_nodes_many(graph, id_names, confirmed_proteins, missing_domains):
    nodes = {'original': [], 'high': [], 'mid': [], 'low': []}
    edges = {'original': [], 'high': [], 'mid': [], 'low': []}
    transcript_names = [id_name[0] for id_name in id_names]
    entrez_ids = [id_name[1] for id_name in id_names]

    for node in graph.nodes:
        if len(node.split('/')) == 2:
            nodes['original'].append(domain_node(node, missing_domains))
        else:
            try:
                nodes['original'].append(protein_node(node, entrez_ids, transcript_names))
            except:
                continue

    for edge in graph.edges(data=True):
        source = edge[0]
        target = edge[1]
        confidence = edge[2].get('confidence', 'original')
        if source.split('/')[0] not in entrez_ids or target.split('/')[0] not in entrez_ids:
            continue
        color = edge_color(edge)
        edges[confidence].append(f'{{from: "{source}", to: "{target}", value: "1", dashes: false, physics: true, {color}}},')

    return nodes, edges


def protein_node(node, entrez_ids, transcript_names):
    node_entrez_id = node
    label = transcript_names[entrez_ids.index(node_entrez_id)]
    return f'{{id: "{node}", label: "{label}", group: "protein", physics: true, source: "PPI", value: "4"}},'


def domain_node(node, missing_domains):
    label = node.split('/')[1]
    color = ''
    if node in missing_domains:  color = 'color: missing,'
    return f'{{id: "{node}", label: "{label}", {color} group: "domain", physics: true, source: "DDI", value: "4"}},'


def edge_color(edge):
    e1 = edge[0].split('/')
    e2 = edge[1].split('/')
    # protein to protein
    if len(e1) == 1 and len(e2) == 1:
        edge_color = ''

        option = 'length: PR_LENGTH,' + edge_color + ' width: WIDTH_SCALE * 4'

    # domain to domain
    elif len(e1) == 2 and len(e2) == 2:
        try:
            if edge[2]['confidence'] != 'original':
                color = 'LIGHTGREEN'
            else:
                color = 'GREEN'
        except KeyError:
            color = 'GREEN'
        option = f'length: PR_DM, color: {color}, width: WIDTH_SCALE'

    # domain to protein
    else:
        option = 'length: LENGTH_domain, color: YELLOW, width: WIDTH_SCALE * 2'
    return option