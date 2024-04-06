import os
import re
import timeit

import pandas as pd

mart_table = pd.read_csv('../sourcedata/mart_export.txt', sep='\t', dtype={'Gene stable ID': str,
                                                                        'UniProtKB/Swiss-Prot ID': str,
                                                                        'NCBI gene (formerly Entrezgene) ID': str})

mart_table.rename(columns={'Gene stable ID': 'GeneID', 'UniProtKB/Swiss-Prot ID': 'UniProtID',
                           'NCBI gene (formerly Entrezgene) ID': 'EntrezID'}, inplace=True)


uniprot_to_entrez_dict = mart_table.set_index('UniProtID')['EntrezID'].to_dict()

pattern = re.compile(r"(PF.*\d)\t(.*)_(.*)\t(PF.*)")


def get_dom_seq(source_address: str):
    # first get seq_dom and then invert it
    seqDom = dict()
    for src in ['pfam-seq-sp', 'pfam-seq-tr']:
        file1 = open(source_address + src, 'r')
        for line in file1:
            line_sp = line.rstrip().split('\t')
            dom = line_sp[0]
            seq = line_sp[1]
            if seq in seqDom:
                seqDom[seq].add(dom)
            else:
                seqDom[seq] = set()
                seqDom[seq].add(dom)
        print(f"Read {src}")
    # reverse the dictionary
    dom_seq = dict()
    for seq in seqDom:
        for dom in seqDom[seq]:
            if dom in dom_seq:
                dom_seq[dom].add(seq)
            else:
                dom_seq[dom] = {seq}
    return dom_seq


def read_interactions(file_path: str, third_col=None):
    interactions = set()
    with open(file_path, 'r') as f:
        # skip header
        f.readline()
        for line in f.readlines():
            line = line.strip().split("\t")
            if third_col is not None:
                if line[2] not in third_col:
                    continue
            assoc = (line[0], line[1]) if line[0] < line[1] else (line[1], line[0])
            # assoc = (line[0], line[1])
            interactions.add(assoc)
    return interactions


def add_classification(file_path: str, classifications_info: str):
    print("Reading the current PPI/DDI interactions")
    interactions = []
    with open(file_path, 'r') as f:
        for line in f.readlines():
            domain1_ppi = line.split("\t")[0]
            domain2_ppi = line.strip().split("\t")[1]
            domain1 = domain1_ppi.split("/")[1]
            domain2 = domain2_ppi.split("/")[1]
            acc = (domain1_ppi, domain2_ppi) if domain1 < domain2 else (domain2_ppi, domain1_ppi)
            interactions.append(acc)

    print("Reading all the DDI classes")
    output_ddis = []
    tmp = set()
    classes = dict()
    with open(classifications_info, 'r') as f:
        for line in f.readlines():
            if not line.startswith("PF"):
                continue
            line_split = line.strip().split("\t")
            domain1 = line_split[0]
            domain2 = line_split[1]
            ddi_class = line_split[-2]
            classes[(domain1, domain2)] = ddi_class

    for i, j in interactions:
        output_ddis.append((i, j, classes[(i.split("/")[1], j.split("/")[1])]))
        tmp.add((i, j))

    if len(output_ddis) != len(interactions):
        print(f"Something has gone wrong, not all DDIs in output anymore ({len(output_ddis)}/{len(interactions)}), "
              f"aborting")
        print(list(set(interactions) - set(tmp))[:5])
        return
    print("Writing the classes to file")
    with open(file_path, 'w') as f:
        f.write('domain_1\tdomain_2\tclass\n')
        for interact in output_ddis:
            f.write(f"{interact[0]}\t{interact[1]}\t{interact[2]}\n")


def process_interaction(dd_interaction: tuple, dom_seq: dict, all_ppis: set):
    try:
        proteins_a = dom_seq[dd_interaction[0]]
        proteins_b = dom_seq[dd_interaction[1]]
    except KeyError:
        return set()

    # replace uniprot ids with entrez ids
    proteins_a = {(protein, uniprot_to_entrez_dict.get(protein)) for protein in proteins_a}
    proteins_b = {(protein, uniprot_to_entrez_dict.get(protein)) for protein in proteins_b}

    dd_pp_interaction = set()
    for prot_a in proteins_a:
        for prot_b in proteins_b:
            # no mapping = can't be used
            if prot_a[1] is None or prot_b[1] is None:
                continue
            # not a ppi = can't have a ddi
            if (prot_a[0], prot_b[0]) not in all_ppis and (prot_b[0], prot_a[0]) not in all_ppis:
                continue
            dd_pp_interaction.add(tuple([f"{prot_a[1]}/{dd_interaction[0]}", f"{prot_b[1]}/{dd_interaction[1]}"]))
    return dd_pp_interaction


def main():
    print("Starting the DDI & PPI combination")
    if not os.path.isfile("../resultdata/source_combined"):
        print("Combining source files for available PPIs")
        os.system("cat ../sourcedata/source* | sort | uniq > ../resultdata/source_combined")
    file_path = "../resultdata/source_combined"
    inter_predicted = read_interactions("../resultdata/result-all")
    dom_seq = get_dom_seq("../sourcedata/")
    all_ppis = read_interactions(file_path)
    print("Read all the necessary data")

    start = timeit.default_timer()

    combined_interactions = set()
    tmp_interactions = []
    c = 0

    for dd_interaction in inter_predicted:
        if c % 1000 == 0:
            print(f"Processed {c}/{len(inter_predicted)} interactions in {round((timeit.default_timer() - start) / 60, 3)} minutes")
        dd_pp_interaction = process_interaction(dd_interaction, dom_seq, all_ppis)
        tmp_interactions.append(dd_pp_interaction)
        c += 1

    for interaction in tmp_interactions:
        combined_interactions.update(interaction)

    print(f"Read {len(combined_interactions)} interactions in {round((timeit.default_timer() - start) / 60, 3)} minutes")

    with open("../resultdata/predicted_ddi_ppi.tsv", 'w') as f:
        for ddi in combined_interactions:
            f.write(f'{ddi[0]}\t{ddi[1]}\n')

    add_classification("../resultdata/predicted_ddi_ppi.tsv", "../resultdata/result-all")


if __name__ == '__main__':
    main()
