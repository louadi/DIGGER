import datetime
import os
import sys
import timeit
from itertools import combinations, chain
from math import sqrt


def read_chain_dom(source_address):
    seqpdbchain = dict()
    pdbchainDom = dict()
    file1 = open(source_address + 'pdb_chain_pfam.tsv', 'r')
    file1.readline()
    file1.readline()
    start1 = timeit.default_timer()
    for line in file1:
        line_sp = line.rstrip().split('\t')
        dom = line_sp[3]
        seq = line_sp[2]
        pdbchain = line_sp[0] + '_' + line_sp[1]
        if pdbchain in pdbchainDom:
            pdbchainDom[pdbchain].append(dom)
        else:
            pdbchainDom[pdbchain] = list()
            pdbchainDom[pdbchain].append(dom)

        if seq in seqpdbchain:
            seqpdbchain[seq].add(pdbchain)
        else:
            seqpdbchain[seq] = set()
            seqpdbchain[seq].add(pdbchain)
    print("Read pdb_chain_pfam.tsv")

    seqDom = dict()
    file1 = open(source_address + 'pfam-seq-sp', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        dom = line_sp[0]
        seq = line_sp[1]
        if seq in seqDom:
            seqDom[seq].append(dom)
        else:
            seqDom[seq] = list()
            seqDom[seq].append(dom)
    print("Read pfam-seq-sp")

    file1 = open(source_address + 'pfam-seq-tr', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        dom = line_sp[0]
        seq = line_sp[1]
        if seq in seqDom:
            seqDom[seq].append(dom)
        else:
            seqDom[seq] = list()
            seqDom[seq].append(dom)
    print("Read pfam-seq-tr")

    file1.close()
    print("Iterating through files took", round(timeit.default_timer() - start1, 1), "seconds")
    print("Length of seqDom, seqpdbchain, pdbchainDom:", len(seqDom), len(seqpdbchain), len(pdbchainDom))
    return seqDom, seqpdbchain, pdbchainDom


def similarity_calculator_interaction(source, domain, seqDom, seqpdbchain, pdbchainDom, source_address, result_address,
                                      redo=False):
    if not redo and os.path.exists(result_address + domain + '-' + domain + '-interaction-' + source[8:]):
        print("Already done", source)
        return
    start = datetime.datetime.now()
    print(f"SIMILARITY Calculator for {source} Interactions")

    result = open(result_address + source + 'pfam', 'w')
    file1 = open(source_address + source, 'r')
    each_interaction_seq_seq = dict()
    interaction_dict = dict()

    for line in file1:
        if '_' in line:
            continue
        line_sp = line.rstrip().split('\t')
        item1 = line_sp[0]
        item2 = line_sp[1]

        # if (item2, item1) in used:
        #     continue
        domains1 = []
        domains2 = []

        if item1 in seqDom and item2 in seqDom:
            domains1: list = seqDom[item1]
            domains2: list = seqDom[item2]
            # comment out the lines below for the ordered data
            if item1 in seqpdbchain and item2 in seqpdbchain:
                pdbchains1 = seqpdbchain[item1]
                pdbchains2 = seqpdbchain[item2]
                for pc in pdbchains1:
                    domains1 = domains1 + pdbchainDom[pc]
                for pc in pdbchains2:
                    domains2 = domains2 + pdbchainDom[pc]


        else:
            # comment out the lines below for the ordered data
            if item1 in seqpdbchain and item2 in seqpdbchain:
                pdbchains1 = seqpdbchain[item1]
                pdbchains2 = seqpdbchain[item2]
                for pc in pdbchains1:
                    domains1 = domains1 + pdbchainDom[pc]
                for pc in pdbchains2:
                    domains2 = domains2 + pdbchainDom[pc]

            else:
                continue

        domains1 = list(set(domains1))
        domains2 = list(set(domains2))

        subsets1 = chain(*map(lambda x: combinations(domains1, x), range(1, 3)))
        subsets2 = chain(*map(lambda x: combinations(domains2, x), range(1, 3)))
        subsets1 = list(subsets1)
        subsets2 = list(subsets2)

        # this is for ordered data, used later when extending DIGGER graphs
        seq_seq = item1 + '_' + item2 if item1 < item2 else item2 + '_' + item1
        # seq_seq = item1 + '_' + item2

        for set1 in subsets1:
            for set2 in subsets2:
                interaction1 = ''
                interaction2 = ''
                set1 = set(set1)
                set2 = set(set2)
                for datum in set1:
                    interaction1 += datum + '_'
                for datum in set2:
                    interaction2 += datum + '_'

                # change this to have it in domain_protein ordering instead of <
                # interaction = (interaction1[:-1], interaction2[:-1])
                if interaction1 < interaction2:
                    interaction = (interaction1[:-1], interaction2[:-1])
                else:
                    interaction = (interaction2[:-1], interaction1[:-1])

                result.write(str(interaction[0]) + '\t' + seq_seq + '\t' + str(interaction[1]) + '\n')

                if interaction[0] in each_interaction_seq_seq:
                    each_interaction_seq_seq[interaction[0]].add(seq_seq)
                else:
                    each_interaction_seq_seq[interaction[0]] = set()
                    each_interaction_seq_seq[interaction[0]].add(seq_seq)

                if interaction[1] in each_interaction_seq_seq:
                    each_interaction_seq_seq[interaction[1]].add(seq_seq)
                else:
                    each_interaction_seq_seq[interaction[1]] = set()
                    each_interaction_seq_seq[interaction[1]].add(seq_seq)

                if interaction in interaction_dict:
                    interaction_dict[interaction].add(seq_seq)
                else:
                    interaction_dict[interaction] = set()
                    interaction_dict[interaction].add(seq_seq)

    file1.close()
    result.close()

    result_file = open(result_address + domain + '-' + domain + '-interaction-' + source[8:], 'w')

    for interaction in interaction_dict:
        interctor1 = interaction[0]
        interctor2 = interaction[1]
        intersection = interaction_dict[interaction]
        nominator = len(intersection)
        denom1 = len(each_interaction_seq_seq[interctor1])
        denom2 = len(each_interaction_seq_seq[interctor2])

        similarity = float(nominator) / (sqrt(denom1) * sqrt(denom2))
        result_file.write(f"{interctor1}\t{interctor2}\t{denom1}\t{denom2}\t{nominator}\t{similarity}\n")
    result_file.close()
    end = datetime.datetime.now()
    print("finished running similarity calculations for", source)
    print("Running Time: " + str(end - start) + "\n")
