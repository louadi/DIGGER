# this aims to recreate the database mapping (Human, Mouse) used in nease
import tqdm
import pandas as pd
import pickle5 as pickle


def load_df(path):
    with open(path, 'rb') as f:
        return pickle.load(f)


human_ground_truth = load_df('db_data/Human')
mouse_ppi_interface = pd.read_csv('../container/domain/data/Mus musculus[mouse]/PPI_interface_mapped_to_exon.csv')
biomart_exons = pd.read_csv('db_data/exon_info_mouse.txt', sep='\t')
biomart_domains = pd.read_csv('db_data/domain_info_mouse.txt', sep='\t')

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
print(new_pdb_df.shape)
