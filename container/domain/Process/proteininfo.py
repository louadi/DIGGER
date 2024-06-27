import re

from domain.Process import process_data as pr
from domain.Process import exonstodomain as exd
from domain.Process import network_analysis as nt

import os
from django.conf import settings
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('agg')

# Global image path
images_path = os.path.join(settings.MEDIA_ROOT, 'images/')
if not os.path.exists(images_path):
    os.makedirs(images_path)


def get_protein_info(ID, organism):
    #Get data about the input transcript or the protein
    info, trID = ID_mapper(ID, organism)
    if info == 0: return 0
    (
        domains, unique, exons, text1, domainshtml, Text_nodes, text_edges, tran_name, gene_name, Ensemble_geneID,
        entrezID,
        gene_description, exons, droped1, droped2, exons_in_interface, co_partners) = info
    #save Image of protein Structure

    Protein_structure(trID, exons, domains, images_path, trID, exons_in_interface)

    return (domains, unique, exons, text1, domainshtml, Text_nodes, text_edges, tran_name, gene_name, Ensemble_geneID,
            entrezID, gene_description, exons, droped1, droped2, trID, images_path, co_partners)


def get_protein_info2(ID):
    #Get data about the input transcript or the protein
    info, trID = ID_mapper(ID)
    if info == 0: return 0
    (domains, unique, exons, text1, domainshtml, Text_nodes, text_edges, tran_name, gene_name, Ensemble_geneID,
     entrezID, gene_description, exons, droped1, droped2, exons_in_interface, co_partners) = info

    return (domains, unique, exons, text1, domainshtml, Text_nodes, text_edges, tran_name, gene_name, Ensemble_geneID,
            entrezID, gene_description, exons, droped1, droped2, trID, co_partners)


def ID_mapper(ID, organism):
    # The use input ID of a transcript or protein:
    #if the ID is a transcript ID:
    if re.match(r"ENS\w*T\d+$", ID):
        return exd.input_transcript(ID, organism), ID
    trID = nt.pr_to_tr(ID, organism)
    if re.match(r"ENS\w*P\d+$", ID):
        return exd.input_transcript(trID, organism), trID


def Visualize_transciript(exon_table, domain_table, exons_in_interface):
    features1 = []
    features2 = []
    for st, e, rank, idd in zip(exon_table["CDS start"], exon_table["CDS end"], exon_table["Exon rank in transcript"],
                                exon_table["Exon stable ID"]):
        if not np.isnan(st):
            if idd not in exons_in_interface:
                features1.append(
                    GraphicFeature(ax=1, start=st / 3, end=e / 3, color="#ffd700", label=str(rank)))
                fend = e / 3
            else:
                features1.append(
                    GraphicFeature(ax=1, start=st / 3, end=e / 3, color="#FF9200", label=str(rank)))
                fend = e / 3
    domain_table = domain_table[["Pfam ID", "Pfam start", "Pfam end", "Interactions mediated by the domain"]]
    domain_table = domain_table.drop_duplicates()
    for st, e, i in zip(domain_table["Pfam start"], domain_table["Pfam end"], domain_table["Pfam ID"]):
        if not np.isnan(st):
            features2.append(
                GraphicFeature(ax=2, start=st, end=e, color="#ffcccc", label=i))
    return features1, features2, fend


def Protein_structure(ID, exons, domains, path, trID, exons_in_interface):
    #save Image of protein Structure
    features1, features2, fend = Visualize_transciript(exons, domains, exons_in_interface)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(16, 3.5))

    record = GraphicRecord(sequence_length=fend, features=features1, )
    record.plot(ax=ax1, figure_width=23, with_ruler=False)

    record = GraphicRecord(sequence_length=fend, features=features2, )
    record.plot(ax=ax2, figure_width=23, with_ruler=True, annotate_inline=True)

    ax1.title.set_text('Coding Exons')
    ax1.title.set_position([.5, -0.4])
    ax2.title.set_text('Pfam Domains')
    ax2.title.set_position([.5, -0.5])

    fig.savefig(path + trID, bbox_inches='tight')
    plt.clf()
    plt.close(fig)
    return
