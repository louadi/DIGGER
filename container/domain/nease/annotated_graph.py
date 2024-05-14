# This file is responsible for filtering the PPI graph based on DDI, ELM and PDB
import pickle
import timeit

import networkx as nx
import pickle5
import pandas as pd
# load files


def load_file(file):
    with open(file, 'rb') as f:
        try:
            return pickle.load(f)
        except:
            return pickle5.load(f)


# files I need: graph, ELM_interactions, pdb, Human
def load_needed_files(organism):
    # domaing = load_file(f'../../../container/domain/data/{organism}/DomainG.pkl')
    domaing = load_file(f'data/{organism}/graph.pkl')
    ppi_graph = load_file(f'data/{organism}/PPI.pkl')
    ELM_interactions = load_file(f'data/{organism}/ELM_interactions')
    pdb = load_file(f'data/{organism}/pdb')
    trivial_name = organism.split('[')[1][:-1].capitalize()
    Organsim = load_file(f'data/{organism}/{trivial_name}')
    return domaing, ppi_graph, ELM_interactions, pdb, Organsim


def ppi_interactions(ppi_graph):
    interactions = set()
    for edge in ppi_graph.edges:
        sorted_edge = tuple(sorted(edge))
        interactions.add(sorted_edge)
    return interactions


def filter_by_ddi(ddi_graph, ppis):
    print(len(ddi_graph.edges))
    ddi_supported_ppis = set()
    for edge in ddi_graph.edges(data=True):
        entrez_only = [edge[0].split("/")[0], edge[1].split("/")[0]]
        entrez_only_ordered = tuple(sorted(entrez_only))
        ddi_supported_ppis.add(entrez_only_ordered)
    return ppis & ddi_supported_ppis


def filter_by_elm(elm, ppis):
    # get sorted tuples
    elm_sorted = [tuple(sorted([str(row[0]), str(row[1])]))
                  for row in zip(elm['Interator gene 1'], elm['Interator gene 2'])]
    elm = set(elm_sorted)
    return elm & ppis


def filter_by_pdb(pdb, organism, ppis):
    # convert organism['gene name'] and organims['entrez'] to dictionary
    organism['NCBI gene ID'] = organism['NCBI gene ID'].astype(str)
    # gene_to_entrez = organism.set_index('Gene name')['NCBI gene ID'].to_dict()
    gene_to_entrez = pdb.set_index('symbol')['entrezgene'].to_dict()

    # convert pdb['symbol'] to entrez and put it in an ordered tuple with pdb['entrezgene']
    # pdb_entrez = [tuple(sorted([str(gene_to_entrez.get(row['symbol'], None)), str(row['entrezgene'])]))
    #               for index, row in pdb.iterrows()]
    pdb_entrez = [tuple(sorted([str(gene_to_entrez.get(row['Gene name'], None)), str(row['NCBI gene ID'])]))
                  for index, row in organism.iterrows()]

    return ppis & set(pdb_entrez)


def filter_ppi_graph(ppis, ddi_graph, elm, pdb, Organism):
    # filter ppi_graph
    print("Filtering PPI graph, current number of edges: ", len(ppis), end=" ")
    ddi_filtered = filter_by_ddi(ddi_graph, ppis)
    # print("Filtered by ddi: ", len(ddi_filtered))
    # filter ELM_interactions
    elm_filtered = filter_by_elm(elm, ppis)
    # print("Filtered by elm: ", len(elm_filtered))
    # filter pdb
    pdb_filtered = filter_by_pdb(pdb, Organism, ppis)
    # print("Filtered by pdb: ", len(pdb_filtered))
    ppis = ddi_filtered | elm_filtered | pdb_filtered
    print(f"Number of edges after filtering: {len(ppis)}, (ddi: {len(ddi_filtered)}, elm: {len(elm_filtered)}, "
          f"pdb: {len(pdb_filtered)})")
    return ppis


def pathway_node_degree(annotated_ppi_graph, pathway_entrez):
    # get the degree of each node in the pathway
    degree = sum([val for _, val in annotated_ppi_graph.degree(pathway_entrez)])
    return degree


if __name__ == '__main__':
    organism = "Mus musculus[mouse]"
    domaing, ppi_graph, ELM_interactions, pdb, Organsim = load_needed_files(organism)
    ppis = ppi_interactions(ppi_graph)
    ppis = filter_ppi_graph(ppis, domaing, ELM_interactions, pdb, Organsim)
    # load some original files
    organism_true = load_file(f'/home/elias/NEASE/nease/data/database/Mouse')
    pdb_true = load_file(f'/home/elias/NEASE/nease/data/database/pdb_mouse')

    print("Number of edges after filtering: ", len(ppis))
