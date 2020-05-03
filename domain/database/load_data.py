import pandas as pd

from os import path
from sqlalchemy import create_engine

# --- Create database connection aka 'SQLAlchemie engine'
# sqlite://<no_hostname>/<path>
# where <path> is relative:
engine = create_engine('sqlite:///datasets.db')

# --- Import the datasets into dataframes
data_base_path = '../data/'

# PPI_interface_mapped_to_exon
# final
# ToDo: all_Proteins
# ToDo: gene_info.csv
# ToDo: gene_to_transcripts

datasets = {
    'ppi_data': pd.read_csv(path.join(data_base_path, 'PPI_interface_mapped_to_exon.csv')),
    'exons_to_domains_data': pd.read_csv(path.join(data_base_path, 'final.csv')),
    'all_Proteins': pd.read_csv(path.join(data_base_path, 'all_Proteins.csv')),
}

# --- Write dataframes to tables in database
for table_name, data_df in datasets.items():
    data_df.to_sql(table_name, engine, if_exists='replace', index=False)
