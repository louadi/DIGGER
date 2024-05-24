# this aims to recreate the database mapping (Human, Mouse) used in nease
import sys

import tqdm
import pandas as pd
import pickle5 as pickle


def load_df(path):
    with open(path, 'rb') as f:
        return pickle.load(f)


human_ground_truth = load_df('db_data/Human')
mouse_ppi_interface = pd.read_csv('db_data/PPI_interface_mapped_to_exon_mouse.csv')
biomart_exons = pd.read_csv('db_data/exon_info_mouse.txt', sep='\t')
biomart_domains = pd.read_csv('db_data/domain_info_mouse.txt', sep='\t')

old_pdb_file = load_df('db_data/pdb')
new_pdb_file = load_df('db_data/pdb_mouse.pkl')


# for PPI interface to exon mapping
def create_ppi_interface_to_exon(exon_info):
    transcript_to_uniprot = pd.read_csv('db_data/transcript_to_uniprot.txt', sep='\t')
    mouse_proteome_pos_pos = pd.read_csv('db_data/mouse_proteome.position_positions_interactions.tsv', sep='\t')

    dict1 = transcript_to_uniprot.set_index('UniProtKB/Swiss-Prot ID')['Transcript stable ID'].to_dict()
    dict2 = transcript_to_uniprot.set_index('UniProtKB/TrEMBL ID')['Transcript stable ID'].to_dict()

    trans_to_uni = [dict1, dict2]

    biomart_exons['CDS start'] = biomart_exons['CDS start'].fillna(-1).astype(int)
    biomart_exons['CDS end'] = biomart_exons['CDS end'].fillna(-1).astype(int)
    exon_info = exon_info.set_index('Transcript stable ID').sort_index()
    total_len = len(mouse_proteome_pos_pos)
    interaction_list = [['Transcript stable ID_x', 'u_ac_1', 'Exon stable ID_x', 'Transcript stable ID_y',
                         'u_ac_2', 'Exon stable ID_y']]
    for _, row in tqdm.tqdm(mouse_proteome_pos_pos.iterrows(), total=total_len):
        interact_a = row['Input Protein ID A']
        interact_b = row['Input Protein ID B']
        pos_a = row['Position A']
        pos_b = row['Position B']

        # check if interact_a is in dict1 or dict2
        transcript_a = trans_to_uni[0].get(interact_a, None)
        if not transcript_a:
            transcript_a = trans_to_uni[1].get(interact_a, None)
            if not transcript_a:
                continue
        transcript_b = trans_to_uni[0].get(interact_b, None)
        if not transcript_b:
            transcript_b = trans_to_uni[1].get(interact_b, None)
            if not transcript_b:
                continue
        # check if transcript_a and transcript_b are in exon_info
        if transcript_a not in exon_info.index or transcript_b not in exon_info.index:
            continue

        # get the exon info for transcript_a and transcript_b
        exon_info_a = exon_info.loc[transcript_a]
        exon_info_b = exon_info.loc[transcript_b]

        if isinstance(exon_info_a, pd.Series):
            exon_info_a = exon_info_a.to_frame().T
        if isinstance(exon_info_b, pd.Series):
            exon_info_b = exon_info_b.to_frame().T

        exon_a, exon_b = None, None

        for _, exon_row in exon_info_a.iterrows():
            # check if CDS start and CDS end are not -1
            if exon_row['CDS start'] == -1 or exon_row['CDS end'] == -1:
                continue

            # check if pos_a is in the range of CDS start and CDS end
            if exon_row['CDS start'] <= pos_a <= exon_row['CDS end']:
                exon_a = exon_row['Exon stable ID']
                break
        for _, exon_row in exon_info_b.iterrows():
            if exon_row['CDS start'] == -1 or exon_row['CDS end'] == -1:
                continue
            if exon_row['CDS start'] <= pos_b <= exon_row['CDS end']:
                exon_b = exon_row['Exon stable ID']
                break

        if exon_a and exon_b:
            interaction_list.append([transcript_a, interact_a, exon_a, transcript_b, interact_b, exon_b])
    return pd.DataFrame(interaction_list[1:], columns=interaction_list[0])


#interaction_list = create_ppi_interface_to_exon(transcript_to_uniprot, biomart_exons, mouse_proteome_pos_pos)
#print(len(interaction_list))
#interaction_list.to_csv('db_data/PPI_interface_mapped_to_exon_mouse.csv', index=False)

# remove the first column in human_ground_truth
human_ground_truth = human_ground_truth.drop(columns=['Unnamed: 0.1'])


# remove 'Exon region start (bp)', 'Exon region end (bp)' from biomart_exons
biomart_exons = biomart_exons.drop(columns=['Exon region start (bp)', 'Exon region end (bp)'])

# join exons and domains on 'Transcript stable ID'
biomart_table = pd.merge(biomart_exons, biomart_domains, on='Transcript stable ID', how='inner')
biomart_table = biomart_table.drop_duplicates().dropna()

# rename NCBI gene (formerly Entrezgene) to NCBI gene ID
biomart_table = biomart_table.rename(columns={'NCBI gene (formerly Entrezgene) ID': 'NCBI gene ID'})

# remove all rows where 'Exon stable ID' hast multiple NCBI gene IDs
grouped = biomart_table.groupby('Exon stable ID')['NCBI gene ID'].nunique()
remove_exon_ids = set(grouped[grouped > 1].index)

# Filter out the rows where 'Exon stable ID' is in remove_exon_ids
biomart_table = biomart_table[~biomart_table['Exon stable ID'].isin(remove_exon_ids)]

# reorder columns to be consistent with the original table
biomart_table = biomart_table[list(human_ground_truth.columns)]

# convert the dtypes to be consistent with the original table
biomart_table = biomart_table.astype(human_ground_truth.dtypes)


# Filter biomart_table down to only the needed exon ids
needed_exon_ids = set(mouse_ppi_interface['Exon stable ID_x']) | set(mouse_ppi_interface['Exon stable ID_y'])
biomart_table_filtered = biomart_table[biomart_table['Exon stable ID'].isin(needed_exon_ids)]

# Set and sort the index for efficient lookups
exon_db_map = biomart_table_filtered.set_index('Exon stable ID').sort_index()

# Filter rows where 'Exon stable ID_x' and 'Exon stable ID_y' are in the filtered biomart_table
existing_exon_ids = set(exon_db_map.index)
filtered_human_ppi_interface = mouse_ppi_interface[
    mouse_ppi_interface['Exon stable ID_x'].isin(existing_exon_ids) &
    mouse_ppi_interface['Exon stable ID_y'].isin(existing_exon_ids)
]

# Initialize new_pdb list with header
new_pdb = [['symbol', 'Gene stable ID', 'Exon with residue ID', 'Genomic coding start', 'Genomic coding end', 'entrezgene']]

# Iterate through each row of filtered_human_ppi_interface
for _, row in tqdm.tqdm(filtered_human_ppi_interface.iterrows(), total=filtered_human_ppi_interface.shape[0]):
    exon_a = row['Exon stable ID_x']
    exon_b = row['Exon stable ID_y']

    # Retrieve the necessary data from exon_db_map and save it
    try:
        db_info_a = exon_db_map.loc[exon_a]
        if len(db_info_a.shape) > 1:
            db_info_a = db_info_a.iloc[0]

        db_info_b = exon_db_map.loc[exon_b]
        if len(db_info_b.shape) > 1:
            db_info_b = db_info_b.iloc[0]

        new_pdb.append([
            db_info_a['Gene name'],
            db_info_a['Gene stable ID'],
            exon_a,
            db_info_a['Genomic coding start'],
            db_info_a['Genomic coding end'],
            db_info_b['NCBI gene ID']
        ])
    except KeyError:
        # Handle missing exon IDs
        continue

# Convert to DataFrame
new_pdb_df = pd.DataFrame(new_pdb[1:], columns=new_pdb[0])
new_pdb_df = new_pdb_df.drop_duplicates()
pickle.dump(new_pdb_df, open('db_data/pdb_mouse.pkl', 'wb'))
print(new_pdb_df.shape)
