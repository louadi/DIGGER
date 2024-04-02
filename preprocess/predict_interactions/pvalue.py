import datetime
import math
from math import sqrt
import sys


def pvalue_calculation(source, seqDom, pdbchainDom, source_address, result_address):
    start = datetime.datetime.now()
    print("P-value calculation for %s (around 2 secs)" % source)

    allPPI = set()
    file1 = open(source_address + source, 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        item1 = line_sp[0]
        item2 = line_sp[1]
        if item1 < item2:
            seq_seq = item1 + '_' + item2
        else:
            seq_seq = item2 + '_' + item1

        allPPI.add(seq_seq)

    N = len(allPPI)

    print(N)

    domain_seq_seq = dict()
    for seq_seq in allPPI:
        seq_split = seq_seq.split('_')
        if len(seq_split) > 2:
            if len(seq_split) == 3:
                print("Make sure there are no _ in the sequence names")
                sys.exit(1)
            seq1 = seq_split[0] + '_' + seq_split[1]
            seq2 = seq_split[2] + '_' + seq_split[3]
        else:
            seq1 = seq_split[0]
            seq2 = seq_split[1]
        if seq1 in seqDom:
            domains = seqDom[seq1]
            for domain in domains:
                if domain in domain_seq_seq:
                    domain_seq_seq[domain].add(seq_seq)
                else:
                    domain_seq_seq[domain] = set()
                    domain_seq_seq[domain].add(seq_seq)

        if seq2 in seqDom:
            domains = seqDom[seq2]
            for domain in domains:
                if domain in domain_seq_seq:
                    domain_seq_seq[domain].add(seq_seq)
                else:
                    domain_seq_seq[domain] = set()
                    domain_seq_seq[domain].add(seq_seq)

            if seq2 in pdbchainDom:
                domains = pdbchainDom[seq2]
                for domain in domains:
                    if domain in domain_seq_seq:
                        domain_seq_seq[domain].add(seq_seq)
                    else:
                        domain_seq_seq[domain] = set()
                        domain_seq_seq[domain].add(seq_seq)

    source_name = source[8:]
    final_result = open(result_address + 'newpvalue-' + source_name, 'w')

    file1 = open(result_address + 'pfam-pfam-interaction-calculated', 'r')
    head = file1.readline().rstrip().split("\t")
    # get index where the source name is in the header
    source_index = head.index(source_name)
    c = 0
    for line in file1:
        lineSplit = line.rstrip().split("\t")
        dom1 = lineSplit[0]
        dom2 = lineSplit[1]
        AssocScore = lineSplit[-1]
        source_score = 0
        if source_index < len(lineSplit):
            source_score = lineSplit[source_index]

        Ne = 0
        Md = 0
        Kde = 0
        if dom1 in domain_seq_seq:
            Ne = len(domain_seq_seq[dom1])
        if dom2 in domain_seq_seq:
            Md = len(domain_seq_seq[dom2])
        if dom1 in domain_seq_seq and dom2 in domain_seq_seq:
            temp = domain_seq_seq[dom1] & domain_seq_seq[dom2]
            Kde = len(temp)

        if Kde == 0 or source_score == '0':
            final_result.write(str(dom1) + "\t" + str(dom2) + "\t" + str(AssocScore) + "\t" + "NA" + "\n")
            c += 1
            # print("-")
            continue
        Nlist = []
        Nelist = []
        Mdlist = []
        #         Kdelist = []

        NMinusMdlist = []
        NMinusNelist = []

        minMdNe = min(Ne, Md)

        coFactorNe = sqrt(2 * (math.pi) * Ne)
        coFactorMd = sqrt(2 * (math.pi) * Md)
        coFactorN = sqrt(2 * (math.pi) * N)
        NminusNe = N - Ne
        NminusMd = N - Md
        coFactorNminusNe = sqrt(2 * (math.pi) * NminusNe)
        coFactorNminusNMd = sqrt(2 * (math.pi) * NminusMd)
        headCoFactors = coFactorNe * coFactorMd * coFactorNminusNe * coFactorNminusNMd
        headCoFactors = math.log10(headCoFactors)

        logNe = math.log10(Ne / (math.e))
        logNe = logNe * Ne
        logMd = math.log10(Md / (math.e))
        logMd = logMd * Md
        logNMinusNe = math.log10(NminusNe / (math.e))
        logNMinusNe = logNMinusNe * NminusNe
        logNMinusMd = math.log10(NminusMd / (math.e))
        logNMinusMd = logNMinusMd * NminusMd
        headLog = logNe + logMd + logNMinusNe + logNMinusMd

        p_value = 0.0
        for i in range(Kde, minMdNe + 1):
            coFactorI = sqrt(2 * (math.pi) * i)
            NeMinusI = Ne - i
            MdMinusI = Md - i
            NMinusMdNePlusI = N - Md - Ne + i
            if NeMinusI == 0:
                coFactorNeMinusI = 1
            else:
                coFactorNeMinusI = sqrt(2 * (math.pi) * NeMinusI)

            if MdMinusI == 0:
                coFactorMdMinusI = 1
            else:
                coFactorMdMinusI = sqrt(2 * (math.pi) * MdMinusI)
            coFactorNMinusMdNePlusI = sqrt(2 * (math.pi) * NMinusMdNePlusI)

            tailCoFactor = coFactorNeMinusI * coFactorI * coFactorMdMinusI * coFactorNMinusMdNePlusI * coFactorN
            tailCoFactor = math.log10(tailCoFactor)

            logI = math.log10(i / (math.e))
            logI = logI * i
            if NeMinusI == 0:
                logNeMinusI = 1
            else:
                logNeMinusI = math.log10(NeMinusI / (math.e))
                logNeMinusI = logNeMinusI * NeMinusI
            if MdMinusI == 0:
                logMdMinusI = 1
            else:
                logMdMinusI = math.log10(MdMinusI / (math.e))
                logMdMinusI = logMdMinusI * MdMinusI
            logN = math.log10(N / (math.e))
            logN = logN * N
            logNMinusMdNePlusI = math.log10(NMinusMdNePlusI / (math.e))
            logNMinusMdNePlusI = logNMinusMdNePlusI * NMinusMdNePlusI

            tailLog = logI + logNeMinusI + logMdMinusI + logN + logNMinusMdNePlusI
            result = headLog + headCoFactors - tailLog - tailCoFactor
            result = 10 ** result
            p_value += result

        # print(p_value)
        if p_value > 1:
            p_value = "1*"
        final_result.write(str(dom1) + "\t" + str(dom2) + "\t" + str(AssocScore) + "\t" + str(p_value) + "\n")
        c += 1
        continue

    final_result.close()
    end = datetime.datetime.now()
    print("Running Time: " + str(end - start) + "\n")


def accumulate_pvalues(sources, result_address):
    source_names = [x[8:] for x in sources]
    pvalue_sources = {x: {} for x in source_names}
    caps_score = dict()

    file1 = open(result_address + 'pfam-pfam-interaction-calculated', 'r')
    for line in file1:
        line_sp = line.rstrip().split("\t")
        caps_score[(line_sp[0], line_sp[1])] = line_sp[-1]

    for source in source_names:
        with open(result_address + 'newpvalue-' + source, 'r') as file1:
            for line in file1:
                line_sp = line.rstrip().split("\t")
                pvalue_sources[source][(line_sp[0], line_sp[1])] = line_sp[3]

    header = 'd1\td2\tcaps_score\t' + '\t'.join(source_names) + '\n'
    with open(result_address + 'newpvalue-all', 'w') as result:
        result.write(header)
        for item in pvalue_sources[source_names[0]]:
            p_value_string = '\t'.join([pvalue_sources[x][item] for x in source_names])
            result.write(item[0] + '\t' + item[1] + '\t' + caps_score[item] + '\t' + p_value_string + '\n')
    result.close()


def gold_silver_bronze(result_address):
    gs = set()
    gs_file = open(result_address + 'pfam-pfam-interaction-goldstandard', 'r')
    for line in gs_file:
        line_sp = line.rstrip().split("\t")
        gs.add((line_sp[0], line_sp[1]))
        gs.add((line_sp[1], line_sp[0]))

    calculated_dict = {}
    length = 0
    calculated = open(result_address + 'pfam-pfam-interaction-calculated', 'r')
    head = calculated.readline().rstrip().split("\t")
    for line in calculated:
        length += 1
        line_sp = line.rstrip().split("\t")
        calculated_dict[(line_sp[0], line_sp[1])] = dict()
        calculated_dict[(line_sp[0], line_sp[1])]['caps'] = line_sp[-1]
        for i in range(2, len(head) - 1):
            calculated_dict[(line_sp[0], line_sp[1])][head[i]] = line_sp[i]

    # write the results
    result = open(result_address + 'result-all', 'w')
    pv = open(result_address + 'newpvalue-all', 'r')
    pvalue_all_head = pv.readline().rstrip().split("\t")
    available_sources = pvalue_all_head[3:]
    header = 'D1\tD2\tSCORE\t' + '\t'.join([f"{x.upper()}_SCORE\t{x.upper()}_PV" for x in available_sources]) + '\tCLASS\tINTERPRO\n'
    result.write(header)
    d_gold_distinct = set()
    d_silver_distinct = set()
    d_bronze_distinct = set()
    for line in pv:
        line_sp = line.rstrip().split("\t")
        d1 = line_sp[0]
        d2 = line_sp[1]
        gs_flag = 'No'
        if (d1, d2) in gs:
            gs_flag = 'Yes'
        caps_score = line_sp[2]
        if caps_score != calculated_dict[(d1, d2)]['caps']:
            print("Something wrong!")
            print(d1, d2, caps_score, calculated_dict[(d1, d2)]['caps'])
            exit()
        # replace this with something that works for all sources
        pvalues_item = {}
        for i, source in enumerate(available_sources):
            pvalues_item[source] = line_sp[i + 3]
        critical_val = 0.05 / length

        how_many_significant = 0
        how_many = 0

        for source in available_sources:
            if pvalues_item[source] != 'NA':
                how_many += 1
                if pvalues_item[source] != '1*' and float(pvalues_item[source]) <= critical_val:
                    how_many_significant += 1

        if how_many >= 4 and how_many == how_many_significant:
            quality = 'Gold'
            d_gold_distinct.add((d1, d2))
        elif how_many < 4 and how_many == how_many_significant:
            quality = 'Silver'
            d_silver_distinct.add((d1, d2))
        elif how_many_significant > 0:
            quality = 'Bronze'
            d_bronze_distinct.add((d1, d2))
        else:
            continue

        line_to_write = f"{d1}\t{d2}\t{caps_score}"
        for source in available_sources:
            line_to_write += f"\t{calculated_dict[(d1,d2)][source]}\t{pvalues_item[source]}"
        line_to_write += f"\t{quality}\t{gs_flag}\n"
        result.write(line_to_write)

    result.close()

    print('total predicted DDIs: gold:', len(d_gold_distinct),
          'silver:', len(d_silver_distinct), 'bronze:', len(d_bronze_distinct))


def one_to_one(result_address):
    result = open(result_address + 'pfam-pfam-interaction-calculated_tuple', 'r')
    doms = []
    interactions_dict = dict()
    result.readline()
    for line in result:
        line_sp = line.rstrip().split("\t")
        doms.append(line_sp[0])
        doms.append(line_sp[1])
        if line_sp[0] in interactions_dict:
            interactions_dict[line_sp[0]].add(line_sp[1])
        else:
            interactions_dict[line_sp[0]] = set()
            interactions_dict[line_sp[0]].add(line_sp[1])

        if line_sp[1] in interactions_dict:
            interactions_dict[line_sp[1]].add(line_sp[0])
        else:
            interactions_dict[line_sp[1]] = set()
            interactions_dict[line_sp[1]].add(line_sp[0])

    interactions = set()
    result = open(result_address + 'pfam-pfam-interaction-calculated_tuple', 'r')
    one_to_one = open(result_address + 'one_to_one_tuple', 'w')
    # line = result.readline()
    # one_to_one.write(line)
    for line in result:
        line_sp = line.rstrip().split("\t")
        if len(interactions_dict[line_sp[0]]) == 1 and len(interactions_dict[line_sp[1]]) == 1:
            interactions.add((line_sp[0], line_sp[1]))
            one_to_one.write(line)

    # print((ok_doms))
    print(len(interactions))
