import datetime


def clean_3did_kbdock_domine_downloaded_files(source_address, result_address):
    start = datetime.datetime.now()
    print("After Downloading Newest version of 3did and KBDOCK, and DOMINE data are to be cleaned.")
    file1 = open(source_address + '3did_flat', 'r')
    result = open(result_address + '3did', 'w')
    associations = set()
    for line in file1:
        if 'PF' not in line:
            continue
        line_sp = line.split('\t')
        PF1 = line_sp[len(line_sp) - 2].lstrip()
        PF1 = PF1[1:8]
        PF2 = line_sp[len(line_sp) - 1]
        PF2 = PF2[:7]
        for line in file1:
            if line.startswith('//'):
                break
            if not line.startswith('#=3D'):
                continue
            line_sp = line.split('\t')
            temp = line_sp[2].split(':')
            chain1 = temp[0]
            temp = line_sp[3].split(':')
            chain2 = temp[0]
            if chain1 != chain2:
                if PF1 < PF2:
                    associations.add((PF1, PF2))
                else:
                    associations.add((PF2, PF1))

    for (pf1, pf2) in associations:
        result.write(pf1 + '\t' + pf2 + '\n')

    result.close()

    file1 = open(source_address + 'INTERACTION.txt', 'r')
    result = open(source_address + 'domine', 'w')
    associations = set()
    for line in file1:
        if 'PF' not in line:
            continue
        line_sp = line.split('|')
        PF1 = line_sp[0]
        PF2 = line_sp[1]
        if PF1 < PF2:
            associations.add((PF1, PF2))
        else:
            associations.add((PF2, PF1))

    for (pf1, pf2) in associations:
        result.write(pf1 + '\t' + pf2 + '\n')
    result.close()

    end = datetime.datetime.now()
    print("Running Time: " + str(end - start) + "\n")


def kbdock_union_3did(source_address, result_address):
    interactions_3did = set()
    pfam_3did = set()
    file1 = open(result_address + '3did', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        item1 = line_sp[0]
        item2 = line_sp[1]
        interactions_3did.add((item1, item2))
        pfam_3did.add(item1)
        pfam_3did.add(item2)

    union = interactions_3did
    print(len(union))
    result = open(result_address + 'unionkbdock3did', 'w')
    for item in union:
        flag = ''
        if item in interactions_3did:
            flag = '3did'
        result.write(item[0] + '\t' + item[1] + '\t' + flag + '\n')

    result.close()


def checkOneOnedomain(result_address):
    did = set()
    file1 = open(result_address + '3did', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        item1 = line_sp[0]
        item2 = line_sp[1]
        if item1 < item2:
            did.add((item1, item2))
        else:
            did.add((item2, item1))

    kbdock = set()
    file1 = open(result_address + 'kbdock', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        item1 = line_sp[0]
        item2 = line_sp[1]
        if item1 < item2:
            kbdock.add((item1, item2))
        else:
            kbdock.add((item2, item1))

    intersection = kbdock & did
    print(len(did))
    print(len(kbdock))
    print(len(intersection))

    caps_score = dict()
    caps = set()
    file1 = open(result_address + 'pfam-pfam-interaction-calculated', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            caps.add((line_sp[0], line_sp[1]))
            caps_score[(line_sp[0], line_sp[1])] = line_sp[9]
        else:
            caps.add((line_sp[1], line_sp[0]))
            caps_score[(line_sp[1], line_sp[0])] = line_sp[9]

    print(len(caps))

    all_interactions = set()
    for source in ['source1_intact', 'source2_mint', 'source3_dip', 'source4_biogrid', 'source5_string-exp',
                   'source5_string-rest', 'source6_sifts', 'source7_hprd']:
        # for source in ['source5_string-exp', 'source5_string-rest']:
        file1 = open(result_address + source, 'r')
        for line in file1:
            if '_' in line:
                continue
            line_sp = line.rstrip().split('\t')
            seq1 = line_sp[0]
            seq2 = line_sp[1]
            all_interactions.add((seq1, seq2))
    print(len(all_interactions))

    seqDom = dict()
    file1 = open(result_address + 'pfam-seq-sp', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        dom = line_sp[0]
        seq = line_sp[1]
        if seq in seqDom:
            seqDom[seq].append(dom)
        else:
            seqDom[seq] = list()
            seqDom[seq].append(dom)

    file1 = open(result_address + 'pfam-seq-tr', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        dom = line_sp[0]
        seq = line_sp[1]
        if seq in seqDom:
            seqDom[seq].append(dom)
        else:
            seqDom[seq] = list()
            seqDom[seq].append(dom)

    counter = set()
    counter2 = set()
    pf_lst = set()
    seq_lst = set()
    for (seq1, seq2) in all_interactions:
        if seq1 in seqDom and seq2 in seqDom and len(seqDom[seq1]) == 1 and len(seqDom[seq2]) == 1:
            if seq1 < seq2:
                counter.add((seq1, seq2))
            else:
                counter.add((seq2, seq1))

            seq_lst.add(seq1)
            seq_lst.add(seq2)

            pf1 = seqDom[seq1][0]
            pf2 = seqDom[seq2][0]
            if (pf1, pf2) in caps or (pf2, pf1) in caps:
                if pf1 < pf2:
                    counter2.add((pf1, pf2))
                else:
                    counter2.add((pf2, pf1))
                pf_lst.add(pf1)
                pf_lst.add(pf2)
                # print(seq1, seq2, pf1, pf2)
                # raw_input()

    print(len(counter))
    print(len(seq_lst))
    print(len(counter2))
    print(len(pf_lst))
    print(len(counter2 & intersection))

    res_file = open(result_address + 'result', 'w')
    file1 = open(result_address + 'result-all', 'r')
    line = file1.readline()
    res_file.write(line.rstrip() + "\tSINGLE_DOMAIN_SEQUENCES\n")
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if (line_sp[0], line_sp[1]) in counter2 or (line_sp[1], line_sp[0]) in counter2:
            res_file.write(line.rstrip() + "\tSingle\n")
        else:
            res_file.write(line.rstrip() + "\tMultiple\n")

    res_file.close()

    print(len(caps))
