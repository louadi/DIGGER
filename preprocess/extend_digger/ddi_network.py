import multiprocessing
import os
from multiprocessing import current_process
import re
import timeit

import pandas as pd

mart_table = pd.read_csv('sourcedata/mart_export.txt', sep='\t', dtype={'Gene stable ID': str,
                                                                        'UniProtKB/Swiss-Prot ID': str,
                                                                        'NCBI gene (formerly Entrezgene) ID': str})

mart_table.rename(columns={'Gene stable ID': 'GeneID', 'UniProtKB/Swiss-Prot ID': 'UniProtID',
                           'NCBI gene (formerly Entrezgene) ID': 'EntrezID'}, inplace=True)


uniprot_to_entrez_dict = mart_table.set_index('UniProtID')['EntrezID'].to_dict()

pattern = re.compile(r"(PF.*\d)\t(.*)_(.*)\t(PF.*)")

inter_predicted = set()

seq_dom = dict()


def get_seq_dom(source_address: str):
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
    return seqDom


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


def process_line(line, predicted_int, uniprot_entrez_map):
    line = line.strip()
    match = pattern.match(line)
    if match is None:
        return None

    domain1, prot1, prot2, domain2 = match.groups()

    try:
        if domain1 not in seq_dom[prot1] or domain2 not in seq_dom[prot2]:
            return None
    except KeyError:
        return None

    if '_' in domain1 or '_' in domain2:
        return None

    interaction = (domain1, domain2) if domain1 < domain2 else (domain2, domain1)
    if interaction not in predicted_int:
        return None

    prot1 = uniprot_entrez_map.get(prot1)
    prot2 = uniprot_entrez_map.get(prot2)
    if prot1 is None or prot2 is None:
        return None
    return f"{prot1}/{domain1}", f"{prot2}/{domain2}"


def worker(chunk):
    global chunks
    results = set()
    for line in chunk:
        result = process_line(line, inter_predicted, uniprot_to_entrez_dict)
        if result is not None:
            results.add(result)
    print(f"chunk complete by worker {current_process().name}")
    return results


def main():
    print("Starting the DDI & PPI combination")
    if not os.path.isfile("resultdata/source_pfam_combined"):
        os.system("cat resultdata/source*pfam | sort | uniq > resultdata/source_pfam_combined")
    file_path = "resultdata/source_pfam_combined"
    global inter_predicted
    global seq_dom
    inter_predicted = read_interactions("resultdata/result-all")
    seq_dom = get_seq_dom("sourcedata/")

    start = timeit.default_timer()

    num_processes = multiprocessing.cpu_count()
    print(f"initializing {num_processes} cores")
    with multiprocessing.Pool(processes=num_processes) as pool, open(file_path, 'r') as f:
        lines = f.readlines()
        chunk_size = len(lines) // num_processes  # Split into roughly equal chunks
        chunks = [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]
        print(f"chunk size: {chunk_size:,} lines per chunk")

        results_sets = pool.map(worker, chunks)

    # Combine results from all chunks
    combined_results = [result for sublist in results_sets for result in sublist]

    # get the results
    interactions = set()
    for i in combined_results:
        interactions.add(i)

    print(f"Read {len(interactions)} interactions in {round((timeit.default_timer() - start) / 60, 3)} minutes from "
          f"{file_path}")

    with open("resultdata/predicted_ddi_ppi.tsv", 'w') as f:
        for ddi in interactions:
            f.write(f'{ddi[0]}\t{ddi[1]}\n')

    add_classification("resultdata/predicted_ddi_ppi.tsv", "resultdata/result-all")


if __name__ == '__main__':
    main()
