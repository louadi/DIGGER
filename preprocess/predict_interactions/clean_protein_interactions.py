# This script cleans protein interaction data from various sources that can then be used in the ppidm-check project
# found at https://github.com/Eeeeelias/ppidm-check. This enables adding predicted domain-domain interactions
import os
import mygene
import pandas as pd


def mygene_query(ensembl_prot_ids, field='uniprot', scope='protein'):
    mg = mygene.MyGeneInfo()
    out = mg.querymany(ensembl_prot_ids, scopes=f'ensembl.{scope}', fields=field, species='mouse')
    mygene_mapping = {}
    for item in out:
        try:
            if 'notfound' in item:
                continue
            if 'uniprot' in item:
                mygene_mapping[item['_id']] = item['uniprot']['Swiss-Prot']
        except KeyError:
            continue
    return mygene_mapping


def clean_mitab(source, target, interactors: list):
    print("Reading:", source)
    interactions = set()
    if 'dip' in target:
        df = pd.read_csv(source, sep='\t', index_col=False)
    else:
        df = pd.read_csv(source, sep='\t')

    try:
        df = df.rename(columns={interactors[0]: 'Alt IDs Interactor A',
                                interactors[1]: 'Alt IDs Interactor B'})
    except KeyError:
        print(f"Columns {interactors} not found, make sure they're spelled correctly")
        exit(1)

    for index, row in df.iterrows():
        interactor_a = row['Alt IDs Interactor A'].split('|')
        interactor_b = row['Alt IDs Interactor B'].split('|')
        curr_interact = []
        for i in interactor_a:
            if 'uniprot' in i:
                curr_interact.append(i.split(':')[1])
        for i in interactor_b:
            if 'uniprot' in i:
                curr_interact.append(i.split(':')[1])
        if len(curr_interact) == 2:
            interactions.add((curr_interact[0], curr_interact[1]))
    print("# Interactions:", len(interactions))
    with open(target, 'w') as f:
        for i in interactions:
            f.write(i[0].replace("_", "-") + '\t' + i[1].replace("_", "-") + '\n')
    return interactions


def clean_mint(source, target):
    print("Reading:", source)
    interactions = set()
    df = pd.read_csv(source, sep='\t', header=None)

    for index, row in df.iterrows():
        try:
            interactor_a = row[0].split(':')[1]
            interactor_b = row[1].split(':')[1]
        except:
            continue
        interactions.add((interactor_a, interactor_b))

    print("# Interactions", len(interactions))
    with open(target, 'w') as f:
        for i in interactions:
            f.write(i[0].replace("_", "-") + '\t' + i[1].replace("_", "-") + '\n')
    return interactions


def clean_string(source, target, mapping):
    print("Reading:", source)
    interactions = set()
    df = pd.read_csv(source, sep=' ')
    df['protein1'] = df['protein1'].str[6:]
    df['protein2'] = df['protein2'].str[6:]

    mapping_df = pd.read_csv(mapping, sep='\t')
    mapping_df = mapping_df[['Protein stable ID', 'UniProtKB/Swiss-Prot ID']]
    # remove rows where there is no mapping
    mapping_df = mapping_df.dropna()
    mapping_dict = mapping_df.set_index('Protein stable ID').squeeze().to_dict()

    # get ids that are not in the mapping
    not_mapped = set(df['protein1']) - set(mapping_dict.keys())
    print("Not mapped by biomart: ", len(not_mapped))
    print("Trying mygene mapping")
    mygene_mapping = mygene_query(list(not_mapped))
    print("Mapped by mygene: ", len(mygene_mapping))
    mapping_dict.update(mygene_mapping)

    for index, row in df.iterrows():
        try:
            interactor_a = mapping_dict[row['protein1']]
            interactor_b = mapping_dict[row['protein2']]
        except:
            continue
        interactions.add((interactor_a, interactor_b))
    print("# Interactions", len(interactions))

    with open(target, 'w') as f:
        for i in interactions:
            f.write(i[0].replace("_", "-") + '\t' + i[1].replace("_", "-") + '\n')
    return interactions


def clean_mippie(source, target, mapping):
    print("Reading:", source)
    interactions = set()
    df = pd.read_csv(source, sep='\t')
    mapping_df = pd.read_csv(mapping, sep='\t')
    mapping_df = mapping_df[['NCBI gene (formerly Entrezgene) ID', 'UniProtKB/Swiss-Prot ID']]
    # remove rows where there is no mapping
    mapping_df = mapping_df.dropna()

    mapping_df['NCBI gene (formerly Entrezgene) ID'] = mapping_df['NCBI gene (formerly Entrezgene) ID'].astype(int)
    mapping_dict = mapping_df.set_index('NCBI gene (formerly Entrezgene) ID').squeeze().to_dict()

    for index, row in df.iterrows():
        try:
            interactor_a = mapping_dict[row['entrezA']]
            interactor_b = mapping_dict[row['entrezB']]
        except:
            continue
        interactions.add((interactor_a, interactor_b))
    print("# Interactions", len(interactions))

    with open(target, 'w') as f:
        for i in interactions:
            f.write(i[0].replace("_", "-") + '\t' + i[1].replace("_", "-") + '\n')
    return interactions


def clean_homology(source, target, mapping):
    print("Reading:", source)
    interactions = set()

    mapping_df = pd.read_csv(mapping, sep='\t')
    mapping_df = mapping_df[['Transcript stable ID', 'UniProtKB/Swiss-Prot ID']]
    # remove rows where there is no mapping
    mapping_df = mapping_df.dropna()
    mapping_dict = mapping_df.set_index('Transcript stable ID').squeeze().to_dict()

    df = pd.read_csv(source, sep='\t')
    unmatched = (set(df['Input Protein ID A']) | set(df['Input Protein ID B'])) - set(mapping_dict.keys())
    matches = (set(df['Input Protein ID A']) | set(df['Input Protein ID B'])) & set(mapping_dict.keys())
    print("Unmatched: ", len(unmatched), "Matches: ", len(matches))
    online_mapped = mygene_query(list(unmatched), scope='transcript')
    print("Mapped by mygene: ", len(online_mapped))
    mapping_dict.update(online_mapped)

    for index, row in df.iterrows():
        interactor_a = row['Input Protein ID A']
        interactor_b = row['Input Protein ID B']
        try:
            interactor_a = mapping_dict[interactor_a]
            interactor_b = mapping_dict[interactor_b]
        except KeyError:
            continue
        interactions.add((interactor_a, interactor_b))
    print("# Interactions", len(interactions))
    with open(target, 'w') as f:
        for i in interactions:
            f.write(i[0].replace("_", "-") + '\t' + i[1].replace("_", "-") + '\n')
    return interactions


def read_all_interactions(path, unique=True):
    if unique:
        interactions = set()
    else:
        interactions = []
    print(f"Reading all interactions in {path}")
    for file in os.listdir(path):
        if not file.startswith('source'):
            continue
        with open(path + file, 'r') as f:
            for line in f:
                split = line.split('\t')
                if split[0] > split[1]:
                    if unique:
                        interactions.add((split[1], split[0]))
                    else:
                        interactions.append((split[1], split[0]))
                else:
                    if unique:
                        interactions.add((split[0], split[1]))
                    else:
                        interactions.append((split[0], split[1]))
    print(f"Total interactions: {len(interactions):,}")


def main(tasks: list):
    for task in tasks:
        function_name = task[0]
        params = task[1:]
        globals()[function_name](*params)

    read_all_interactions("../sourcedata/", True)

