import operator

import datetime
import pickle
import sys
import itertools
import timeit

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import random
from math import sqrt

from typing import Tuple

# reproducibility yay
random.seed(42)


def read_interactions(file_path: str):
    interactions_3did = set()
    pfam_3did = set()
    file1 = open(file_path, 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        item1 = line_sp[0]
        item2 = line_sp[1]
        if item1 < item2:
            interactions_3did.add((item1, item2))
        else:
            interactions_3did.add((item2, item1))
        pfam_3did.add(item1)
        pfam_3did.add(item2)
    return interactions_3did, pfam_3did


def interactions(file_path: str, info: dict, info_tuple_multiple: dict, source: str) -> Tuple[set, set, dict, dict]:
    interactions_intact = set()
    pfam_intact = set()
    file1 = open(file_path, 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            interaction = (line_sp[0], line_sp[1])
        else:
            interaction = (line_sp[1], line_sp[0])

        if '_' not in line:
            interactions_intact.add(interaction)
            if interaction[0] in info:
                if interaction[1] in info[interaction[0]]:
                    info[interaction[0]][interaction[1]][source] = float(line_sp[5])
                else:
                    info[interaction[0]][interaction[1]] = dict()
                    info[interaction[0]][interaction[1]][source] = float(line_sp[5])
            else:
                info[interaction[0]] = dict()
                info[interaction[0]][interaction[1]] = dict()
                info[interaction[0]][interaction[1]][source] = float(line_sp[5])

            pfam_intact.add(line_sp[0])
            pfam_intact.add(line_sp[1])
        else:
            if interaction[0] in info_tuple_multiple:
                if interaction[1] in info_tuple_multiple[interaction[0]]:
                    info_tuple_multiple[interaction[0]][interaction[1]][source] = float(line_sp[5])
                else:
                    info_tuple_multiple[interaction[0]][interaction[1]] = dict()
                    info_tuple_multiple[interaction[0]][interaction[1]][source] = float(line_sp[5])
            else:
                info_tuple_multiple[interaction[0]] = dict()
                info_tuple_multiple[interaction[0]][interaction[1]] = dict()
                info_tuple_multiple[interaction[0]][interaction[1]][source] = float(line_sp[5])
    return interactions_intact, pfam_intact, info, info_tuple_multiple


def coef_score(coefficients: list, interaction: list, sources: list):
    result = []
    for coef, src in zip(coefficients, sources):
        result.append(coef * interaction[src])
    return result


def extract_info(relevant_pfams: set, score_info: dict):
    result = {}
    for pfam_id in relevant_pfams:
        try:
            result[pfam_id] = score_info[pfam_id]
        except KeyError:
            result[pfam_id] = 0
    pickle.dump(result, open('../pickles/info_scores.pickle', 'wb'))
    print("Wrote scores to pickle")


def plot_loss_curve(losses: list):
    plt.plot(losses)
    plt.xlabel('Iterations')
    plt.ylabel('Loss')
    plt.show()


def plot_score_density(positive_scores: list, negative_scores: list, v_lines=False):
    plt.style.use('ggplot')
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    sns.kdeplot(negative_scores, color='orange', ax=ax1)
    if v_lines:
        ax1.axvline(x=0.01586, linestyle='dotted', color='#0097ff', label='threshold (0.01586)')
        ax1.axvline(x=0.03, linestyle='dotted', color='#fc1b80', label='threshold (0.03)')
    ax1.set_ylabel('Negative train - Density')

    sns.kdeplot(positive_scores, color='blue', ax=ax2)
    if v_lines:
        ax2.axvline(x=0.01586, linestyle='dotted', color='#0097ff', label='threshold')
        ax2.axvline(x=0.03, linestyle='dotted', color='#fc1b80')
    ax2.set_ylabel('Positive train - Density')

    plt.xlim(0, 1)

    # Adding labels to the plot
    plt.xlabel('Score')
    plt.suptitle('Score density in negatives vs. positives')

    # Displaying the plot
    plt.tight_layout()
    ax1.legend()
    plt.show()


def plot_density_by_source(source_scores: dict):
    plt.style.use('ggplot')
    fig, ax = plt.subplots((len(source_scores)), 1, sharex=True, figsize=(10, 10))

    for i, source in enumerate(source_scores):
        # sns.kdeplot(source_scores[source], ax=ax[i], label=source)
        sns.histplot(source_scores[source], ax=ax[i], color='blue')
        ax[i].set_ylabel(source)
        ax[i].set_yscale('log')
    plt.xlabel('Score')
    plt.show()


def plot_auc_boxes(all_auc: dict):
    plt.style.use('ggplot')
    fig, ax = plt.subplots()
    ax.boxplot(all_auc.values())
    ax.set_xticklabels(all_auc.keys())
    plt.xlabel('Number of sources')
    plt.ylabel('AUC')
    plt.show()


def calculate_coefficient_auc(coefficients: tuple, info: dict, gold_standard: set, background_data: list, source_names: list):
    coefficients = list(coefficients)
    # print(coefficients)
    # make sure coefficients are integers between 1-100
    # assert all(1 <= coef <= 100 for coef in coefficients)
    # background_data = random.sample(background_data, len(gold_standard))
    background_data = random.sample(background_data, len(gold_standard))
    len_background = len(background_data)
    coef_sum = sum(coefficients)
    all_data_positive_negative = {}
    for datum in gold_standard:
        try:
            dom_score = coef_score(coefficients, info[datum[0]][datum[1]], source_names)
            all_data_positive_negative[datum] = sum(dom_score) / coef_sum
        except KeyError:
            all_data_positive_negative[datum] = 0

    for datum in background_data:
        try:
            dom_score = coef_score(coefficients, info[datum[0]][datum[1]], source_names)
            all_data_positive_negative[datum] = sum(dom_score) / coef_sum
        except KeyError:
            pass

    sorted_all_data_positive_negative = sorted(all_data_positive_negative.items(),
                                               key=operator.itemgetter(1), reverse=True)
    yindex = 0.0
    area_under_curve = 0.0

    for datum in sorted_all_data_positive_negative:
        if datum[0] in gold_standard:
            yindex += 1
        else:
            area_under_curve += (yindex / len(gold_standard)) * (1.0 / len_background)
    return area_under_curve


def calculate_coefficient_sum(coefficients: tuple, info: dict, gold_standard: set, background_data: list, source_names: list):
    # calculate the sum of the scores for the given coefficients
    coefficients = list(coefficients)
    background_data = random.sample(background_data, len(gold_standard))
    len_background = len(background_data)
    coef_sum = sum(coefficients)
    all_positive = 0
    all_background = 0
    for datum in gold_standard:
        try:
            dom_score = coef_score(coefficients, info[datum[0]][datum[1]], source_names)
            all_positive += sum(dom_score) / coef_sum
        except KeyError:
            pass

    for datum in background_data:
        try:
            dom_score = coef_score(coefficients, info[datum[0]][datum[1]], source_names)
            all_background += sum(dom_score) / coef_sum
        except KeyError:
            pass

    all_positive = all_positive / len(gold_standard)
    all_background = all_background / len_background
    return all_positive - all_background


# get random coefficients and calculate AUC
def best_coefficients_rand(info: dict, gold_standard: set, background_data: list, source_names: list):
    num_sources = len(source_names)
    best_coeffs = tuple([1] * num_sources)
    tested = set()
    all_auc = []
    max_auc = [0.0]
    max_coefs = [best_coeffs]
    iterations = 10_000
    i = 0

    while i < iterations:
        coefficients = np.random.randint(1, 101, num_sources)
        coefficients = tuple(coefficients)

        # make sure we don't try the same coefficients again
        if coefficients in tested:
            continue
        # Calculate AUC
        auc = calculate_coefficient_auc(coefficients, info, gold_standard, background_data, source_names)
        if auc > min(max_auc):
            auc_index = max_auc.index(min(max_auc))
            if len(max_auc) > 500:
                max_auc.remove(min(max_auc))
                max_coefs.pop(auc_index)
            max_auc.append(auc)
            max_coefs.append(coefficients)
            best_coeffs = coefficients
        all_auc.append(auc)
        i += 1

    # for each of the best auc get the coefficient_sum and choose the highest one
    final_auc = 0.0
    final_coefs = None
    max_sum = 0.0
    for coef in max_coefs:
        coef_sum = calculate_coefficient_sum(coef, info, gold_standard, background_data, source_names)
        if coef_sum > max_sum:
            final_coefs = coef
            max_sum = coef_sum
            final_auc = calculate_coefficient_auc(coef, info, gold_standard, background_data, source_names)

    print(f"Using {num_sources} sources")
    print(f"Max AUC: {max(all_auc)}\nMin AUC: {min(all_auc)}\nFinal AUC: {final_auc}")
    print(f"Variance: {np.var(all_auc)}")
    print(f"Standard deviation: {np.std(all_auc)}\n")
    print(f"Best coefficients: {final_coefs}")
    return {f'coef{i}': coef for i, coef in enumerate(final_coefs, 1)}, all_auc


def create_wrong_assocations(sources, source_address, result_address):
    start = datetime.datetime.now()
    print("Create NEGATIVE assocaitions from all inputs (around ? mins)")

    # removed kbdock from here

    file1 = open(result_address + '3did', 'r')
    gs = set()
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            interaction = (line_sp[0], line_sp[1])
        else:
            interaction = (line_sp[1], line_sp[0])
        gs.add(interaction)

    common_factors = set()
    dom_common_factors = dict()
    all_interactions_DDI = set()
    # for source in ['source1_intact']:
    for source in sources:
        file1 = open(result_address + source + 'pfam', 'r')
        # interaction_score = dict()
        for line in file1:
            line_sp = line.rstrip().split('\t')
            cf = hash(line_sp[1])
            dom1 = line_sp[0]
            dom2 = line_sp[2]
            if '_' in dom1 or '_' in dom2:
                continue

            all_interactions_DDI.add((dom1, dom2))

            if 'string' in source:
                continue

            common_factors.add(cf)

            if dom1 not in dom_common_factors:
                dom_common_factors[dom1] = set()
            dom_common_factors[dom1].add(cf)

            if dom2 not in dom_common_factors:
                dom_common_factors[dom2] = set()
            dom_common_factors[dom2].add(cf)

    wrongcf_foreach_dom = dict()

    counter = 0
    used = set()
    result_file = open(result_address + 'negative_set', 'w')

    # for every domain that has interactions
    for dom in dom_common_factors:
        counter += 1
        # print(counter)
        # get the interactions of that domain
        number_cf_for_dom = dom_common_factors[dom]
        # remove them from all interactions in general
        all_possible_cf_for_dom = common_factors - number_cf_for_dom
        # randomly get the same amount of interactions from the rest (= node degree)
        wrong_cfs = random.sample(all_possible_cf_for_dom, len(number_cf_for_dom))
        # save it
        wrongcf_foreach_dom[dom] = wrong_cfs

    counter = 0
    for dom1 in dom_common_factors:
        for dom2 in dom_common_factors:
            counter += 1
            if dom1 < dom2:
                interaction = (dom1, dom2)
            elif dom1 > dom2:
                interaction = (dom2, dom1)
            else:
                continue

            if interaction in used:
                continue
            used.add(interaction)

            # making sure the interactions are not in the positive set already
            if interaction in gs or interaction in all_interactions_DDI:
                continue

            # set of all wrong common factors (hashes)
            wrong_cfs1: set = wrongcf_foreach_dom[dom1]

            wrong_cfs2: set = wrongcf_foreach_dom[dom2]

            intersection = set(wrong_cfs1) & set(wrong_cfs2)
            nominator = len(intersection)
            denom1 = len(set(wrong_cfs1))
            denom2 = len(set(wrong_cfs2))

            similarity = float(nominator) / (sqrt(denom1) * sqrt(denom2))
            if similarity != 0:
                result_file.write(f"{dom1}\t{dom2}\t{denom1}\t{denom2}\t{nominator}\t{similarity}\n")
    result_file.close()

    end = datetime.datetime.now()
    print("Running Time: " + str(end - start) + "\n")


def assign_interaction(sources, result_address):
    start = datetime.datetime.now()
    print("Filtering associations for Interactions(around ? mins)")
    interactions_3did, pfam_3did = read_interactions(result_address + '3did')
    # interactions_kbdock, pfam_kbdock = read_interactions(result_address + 'kbdock')
    # interactions_domine, pfam_domine = read_interactions(result_address + 'domine')

    info = dict()
    info_tuple_multiple = dict()

    interactions_sources = {}
    pfam_sources = {}
    source_names = []

    for source in sources:
        source_name = source[8:]
        interaction, pfam, info, info_tuple_multiple = interactions(
            result_address + 'pfam-pfam-interaction-' + source_name, info, info_tuple_multiple, source_name)
        print(f"Got {len(interaction)} interactions for source: {source_name}")
        source_names.append(source_name)
        interactions_sources[source_name] = interaction
        pfam_sources[source_name] = pfam

    for item1 in info:
        for item2 in info[item1]:
            for source in source_names:
                info[item1][item2].setdefault(source, 0)

    print("Getting gold_standard interactions: " + str(datetime.datetime.now() - start) + "\n")
    # find the best interactions that are well-supported by any of the sources
    gold_standard = set()
    for source in interactions_sources.keys():
        gold_standard = gold_standard.union(interactions_sources[source])
    background_data = list(gold_standard - (gold_standard.intersection(interactions_3did)))
    gold_standard = gold_standard.intersection(interactions_3did)

    scores = {}
    for datum in gold_standard:
        try:
            for key in info[datum[0]][datum[1]]:
                if key in scores:
                    scores[key].append(info[datum[0]][datum[1]][key])
                else:
                    scores[key] = [info[datum[0]][datum[1]][key]]
        except KeyError:
            continue
    # plot_density_by_source(scores)

    train_set = random.sample(gold_standard, int(len(gold_standard) / 2))

    test_set = gold_standard - set(train_set)
    print("Length of gold standard:", len(gold_standard))
    print("Length of train set:", len(train_set))
    print("Length of test set:", len(test_set))
    print("Length of background data:", len(background_data))

    print("Calculating coefficients, this will take a while...")
    best_coefs, all_auc = best_coefficients_rand(info, gold_standard, background_data, source_names)

    all_data_scores = dict()
    best_coefs = [x for x in best_coefs.values()]
    coef_summation = sum(best_coefs)
    for item1 in info:
        for item2 in info[item1]:
            dom_score = coef_score(best_coefs, info[item1][item2], source_names)
            if item1 in all_data_scores:
                all_data_scores[item1][item2] = sum(dom_score) / coef_summation
            else:
                all_data_scores[item1] = dict()
                all_data_scores[item1][item2] = sum(dom_score) / coef_summation

    # finding negative set with low-scoring
    print("finding low-scoring associations where score cap be anything")
    neg_model = 1

    gold_standard_negative_set = set()

    if neg_model == 1:
        negatives = set()
        negatives_score = dict()
        negative_file = open(result_address + 'negative_set', 'r')
        for line in negative_file:
            line_sp = line.rstrip().split('\t')
            # if float(line_sp[5]) > 0.04:
            #     continue
            negatives.add((line_sp[0], line_sp[1]))
            negatives_score[(line_sp[0], line_sp[1])] = float(line_sp[5])
            # print(line_sp[5])

        print("Negatives", len(negatives))
        gold_standard_negative_set = random.sample(negatives, len(gold_standard))

    elif neg_model == 2:
        max_pair = ''
        max_score = 0.0
        for datum1 in info:
            for datum2 in info[datum1]:
                if (datum1, datum2) in gold_standard:
                    continue
                flag_ok_for_negative = 0
                for source in source_names:
                    if info[datum1][datum2][source] > 0:
                        flag_ok_for_negative += 1

                if flag_ok_for_negative >= 6:
                    if len(gold_standard_negative_set) < len(gold_standard):
                        gold_standard_negative_set.add((datum1, datum2))
                        max_score = 0.0
                        for item in gold_standard_negative_set:
                            if all_data_scores[item[0]][item[1]] > max_score:
                                max_score = all_data_scores[item[0]][item[1]]
                                max_pair = item
                    else:
                        if all_data_scores[datum1][datum2] < max_score:
                            max_score = all_data_scores[datum1][datum2]
                            # print(max_pair)
                            gold_standard_negative_set.remove(max_pair)
                            gold_standard_negative_set.add((datum1, datum2))
                            max_score = 0.0
                            for item in gold_standard_negative_set:
                                if all_data_scores[item[0]][item[1]] > max_score:
                                    max_score = all_data_scores[item[0]][item[1]]
                                    max_pair = item
        print('negative set max score' + str(max_score))

    else:
        gold_standard_negative_set = random.sample(background_data, len(gold_standard))

    train_negative_set = random.sample(gold_standard_negative_set, int(len(gold_standard) / 2))
    test_negative_set = set(gold_standard_negative_set) - set(train_negative_set)

    print("Length of gold standard negative:", len(gold_standard_negative_set))
    print("Length of negative train:", len(train_negative_set))
    print("Length of negative test:", len(test_negative_set))
    # print(train_negative_set)
    negative_scores = []
    for datum in train_negative_set:
        negative_scores.append(negatives_score[datum])
    positive_scores = []
    for datum in train_set:
        positive_scores.append(all_data_scores[datum[0]][datum[1]])

    plot_score_density(positive_scores, negative_scores)

    best_fmeasure = 0
    best_threshold = 1000
    best_fmeasure_test = 0

    # plot distribution of scores for positives and negatives
    # calculating best Threshold and best F-measure
    count = 1
    for threshold in range(30, 1, -1):
        threshold = float(threshold) / 1000
        count_for_train = 0
        count_for_test = 0
        count_for_train_negative = 0
        count_for_test_negative = 0
        count_all_found = 0
        # Check the training set of interpro and see whether interproTrain association is found or not AND they are more
        # than THRESHOLD SCORE
        for datum in train_set:
            interaction = datum
            score = 0
            flag = False
            if datum[0] in all_data_scores:
                if datum[1] in all_data_scores[datum[0]]:

                    score = all_data_scores[datum[0]][datum[1]]
                    if score >= threshold:
                        # flag = True
                        count_for_train += 1

        for datum in test_set:
            interaction = datum
            score = 0
            # flag = False
            if datum[0] in all_data_scores:
                if datum[1] in all_data_scores[datum[0]]:

                    score = all_data_scores[datum[0]][datum[1]]
                    if score >= threshold:
                        # flag = True
                        count_for_test += 1

        for datum in train_negative_set:
            score = 0
            # flag = False
            if neg_model == 1:
                score = negatives_score[datum]
            else:
                if datum[0] in info:
                    if datum[1] in info[datum[0]]:
                        coef_summation = float(sum(best_coefs))
                        score = sum(coef_score(best_coefs, info[datum[0]][datum[1]], source_names)) / coef_summation
            if score >= threshold:
                # flag = True
                count_for_train_negative += 1


        for datum in test_negative_set:
            score = 0
            # flag = False
            if neg_model == 1:
                score = negatives_score[datum]
            else:
                if datum[0] in info:
                    if datum[1] in info[datum[0]]:
                        coef_summation = float(sum(best_coefs))
                        score = sum(coef_score(best_coefs, info[datum[0]][datum[1]], source_names)) / coef_summation
            if score >= threshold:
                # flag = True
                count_for_test_negative += 1

        # Check the found associations and see whether they are more than THRESHOLD SCORE
        for datum1 in info:
            for datum2 in info[datum1]:
                flag = False
                coef_summation = float(sum(best_coefs))
                score = sum(coef_score(best_coefs, info[datum1][datum2], source_names)) / coef_summation
                if score >= threshold:
                    flag = True
                    count_all_found += 1

        tp_train = count_for_train
        fn_train = len(train_set) - tp_train
        fp_train = count_for_train_negative
        f_score_train = (float(tp_train) * 2) / ((2 * tp_train) + fn_train + fp_train)

        tp_test = count_for_test
        fn_test = len(test_set) - tp_test
        fp_test = count_for_test_negative
        f_score_test = (float(tp_test) * 2) / ((2 * tp_test) + fn_test + fp_test)

        if best_fmeasure < f_score_train:
            best_fmeasure = f_score_train
            best_threshold = threshold
            best_fmeasure_test = f_score_test

        bold = '\033[1m'
        end = '\033[0m'
        print(f"############### Iteration {count} #################")
        print('Training Set'.ljust(50), "| Testing Set")
        print(f"T\t\tPredicted".ljust(45), "| T\t\tPredicted")
        print("T\t\tPos\t   Neg".ljust(45), "| T\t\tPos\t   Neg")
        print(f"T Pos\t{tp_train:4} | {fn_train}".ljust(48), f"| T Pos\t{tp_test} | {fn_test}")
        print(f"T Neg\t{fp_train:4} | / ".ljust(48), f"| T Neg\t{fp_test:4} | / ")
        print(f"T F_score: {round(f_score_train, 5)} | Best F_measure: {bold + str(round(best_fmeasure, 5)) + end}".ljust(58),
              f"| T F_score: {round(f_score_test, 5)} | Best F_measure: {bold + str(round(best_fmeasure_test, 5)) + end}")
        print("-" * 80)
        print(
            f"Current threshold: {threshold} | Best threshold: {bold + str(best_threshold) + end} | Count all found: {count_all_found}\n")
        count += 1

    if best_fmeasure_test < 0.8:
        print(f"Best threshold ({best_fmeasure_test}) is not very good, this may result in poor quality interactions.")
        continue_flag = input("Do you want to continue anyway? (y/n): ")
        if continue_flag == 'n':
            sys.exit()

    result_calculated = open(result_address + 'pfam-pfam-interaction-calculated', 'w')
    result_merged = open(result_address + 'pfam-pfam-interaction-merged', 'w')
    header = 'domain1\tdomain2\t' + '\t'.join(source_names) + '\tscore\n'
    result_calculated.write(header)
    result_merged.write(header)
    for datum1 in info:
        for datum2 in info[datum1]:
            score = all_data_scores[datum1][datum2]
            source_infos = '\t'.join([str(info[datum1][datum2][source]) for source in source_names])
            result_string = datum1 + '\t' + datum2 + '\t' + source_infos + '\t' + str(score) + '\n'
            result_merged.write(result_string)
            if score >= best_threshold:
                flag = True
                result_calculated.write(result_string)

    result_gold_standard = open(result_address + 'pfam-pfam-interaction-goldstandard', 'w')
    header = 'domain1\tdomain2\t' + '\t'.join(source_names) + '\tscore\tflag\n'
    result_gold_standard.write(header)
    for datum in gold_standard:
        if datum in train_set:
            flag = 'train-'
        else:
            flag = 'test-'

        try:
            score = all_data_scores[datum[0]][datum[1]]
        except KeyError:
            score = 0

        if score >= best_threshold:
            flag = flag + 'yes'
        else:
            flag = flag + 'no'

        sources = '\t'.join([str(info[datum[0]][datum[1]][source]) for source in source_names])
        result_gold = datum[0] + '\t' + datum[1] + '\t' + sources + '\t' + str(score) + '\t' + flag + '\n'
        result_gold_standard.write(result_gold)

    result_negative = open(result_address + 'pfam-pfam-interaction-negative', 'w')
    header = 'domain1\tdomain2\tscore\tflag\n'
    result_negative.write(header)
    for datum in gold_standard_negative_set:
        if datum in train_negative_set:
            flag = 'train-'
        else:
            flag = 'test-'

        if neg_model == 1:
            score = negatives_score[datum]
        else:
            coef_summation = sum(best_coefs)
            score = sum(coef_score(best_coefs, info[datum[0]][datum[1]], source_names)) / coef_summation

        if score >= best_threshold:
            flag = flag + 'yes'
        else:
            flag = flag + 'no'

        result_negative.write(
            datum[0] + '\t' + datum[1] + '\t' + str(score) + '\t' + flag + '\n')

    for item1 in info_tuple_multiple:
        for item2 in info_tuple_multiple[item1]:
            for source in source_names:
                info_tuple_multiple[item1][item2].setdefault(source, 0)

    result_calculated_tuple = open(result_address + 'pfam-pfam-interaction-calculated_tuple', 'w')
    result_merged_tuple = open(result_address + 'pfam-pfam-interaction-merged_tuple', 'w')
    for datum1 in info_tuple_multiple:
        for datum2 in info_tuple_multiple[datum1]:
            # flag = False
            coef_summation = sum(best_coefs)
            score = sum(coef_score(best_coefs, info_tuple_multiple[datum1][datum2], source_names)) / coef_summation

            source_infos = '\t'.join([str(info_tuple_multiple[datum1][datum2][source]) for source in source_names])
            result_string = datum1 + '\t' + datum2 + '\t' + source_infos + '\t' + str(score) + '\n'

            result_merged_tuple.write(result_string)

            if score >= best_threshold:
                flag = True
                result_calculated_tuple.write(result_string)

    result_calculated.close()
    result_merged.close()
    result_calculated_tuple.close()
    result_merged_tuple.close()
    result_gold_standard.close()
    result_negative.close()
    end = datetime.datetime.now()
    print("Running Time: " + str(end - start) + "\n")
