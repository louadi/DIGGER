# This script cleans protein interaction data from various sources that can then be used in the ppidm-check project
# found at https://github.com/Eeeeelias/ppidm-check. This enables adding predicted domain-domain interactions
import os
import mygene
import pandas as pd


def mygene_query(ensembl_prot_ids, field='uniprot'):
    mg = mygene.MyGeneInfo()
    out = mg.querymany(ensembl_prot_ids, scopes='ensembl.protein', fields=field, species='mouse')
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


def clean_mitab(source, target, other_spelling=None):
    interactions = set()
    if other_spelling is 'dip':
        df = pd.read_csv(source, sep='\t', index_col=False)
    else:
        df = pd.read_csv(source, sep='\t')
    if other_spelling == 'intact':
        df = df.rename(columns={'Alt. ID(s) interactor A': 'Alt IDs Interactor A',
                                'Alt. ID(s) interactor B': 'Alt IDs Interactor B'})
    if other_spelling == 'dip':
        df = df.rename(columns={'ID interactor A': 'Alt IDs Interactor A',
                                'ID interactor B': 'Alt IDs Interactor B'})
    print(df.columns)
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
    print(len(interactions))
    with open(target, 'w') as f:
        for i in interactions:
            f.write(i[0] + '\t' + i[1] + '\n')
    return interactions


def clean_mint(source, target):
    interactions = set()
    df = pd.read_csv(source, sep='\t', header=None)

    for index, row in df.iterrows():
        try:
            interactor_a = row[0].split(':')[1]
            interactor_b = row[1].split(':')[1]
        except:
            print(row)
            continue
        interactions.add((interactor_a, interactor_b))

    print(len(interactions))
    with open(target + 'source_mint', 'w') as f:
        for i in interactions:
            f.write(i[0] + '\t' + i[1] + '\n')
    return interactions


def clean_string(source, target, mapping):
    interactions = set()
    df = pd.read_csv(source, sep=' ')
    df['protein1'] = df['protein1'].str[6:]
    df['protein2'] = df['protein2'].str[6:]
    print(len(set(df['protein1'])))

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
    print(len(interactions))

    with open(target + 'source_string', 'w') as f:
        for i in interactions:
            f.write(i[0] + '\t' + i[1] + '\n')
    return interactions


def clean_mippie(source, target, mapping):
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
    print(len(interactions))

    with open(target + 'source_mippie', 'w') as f:
        for i in interactions:
            f.write(i[0] + '\t' + i[1] + '\n')
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


if __name__ == '__main__':
    # target path
    target_path = "/home/elias/hamburg/ppidm-check/sourcedata/"
    target_path_human = "/mnt/d/programming/bachelor_projects/ppidm-check/sourcedata/"

    # source paths
    biogrid_source = "/mnt/d/Downloads/mouse_datasets/BIOGRID-ORGANISM-Mus_musculus-4.4.230.mitab.txt"
    mint_source = "/mnt/d/Downloads/mouse_datasets/mint.txt"
    string_source = "/mnt/d/Downloads/mouse_datasets/string_10090.protein.links.v12.0.txt"
    intact_source = "/mnt/d/Downloads/mouse_datasets/intact.txt"
    dip_source = "/mnt/d/Downloads/mouse_datasets/Mmusc20170205.txt"

    mippie_source = "/mnt/d/Downloads/mouse_datasets/mippie_ppi_v1_0.tsv"

    # mapping
    biomart_mapping = "/mnt/d/Downloads/mouse_datasets/mart_export.txt"
    # clean_mitab(biogrid_source, target_path + "source_biogrid")
    # clean_mitab(intact_source, target_path + "source_intact", 'intact')
    clean_mitab(dip_source, target_path + "source_dip", 'dip')
    # clean_mint(mint_source, target_path)
    # clean_string(string_source, target_path, biomart_mapping)
    # clean_mippie(mippie_source, target_path, biomart_mapping)


    read_all_interactions(target_path, True)
    # read_all_interactions(target_path_human, False)