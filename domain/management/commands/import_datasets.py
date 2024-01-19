import os
from os import path

import django
import pandas as pd
from django.conf import settings
from django.core.management import BaseCommand
from tqdm import tqdm

from domain.models import Gene, Domain

# --- Get database connection aka 'SQLAlchemie engine'

engine = settings.DATABASE_ENGINE

# --- Import the datasets into dataframes
data_base_path = path.join(settings.PROJECT_ROOT, 'domain/data')

datasets = {}
for organism in os.listdir('domain/data'):
    if not os.path.isdir('domain/data/' + organism):
        continue
    trivial_name = organism.split("[")[1][:-1]
    datasets['ppi_data_' + trivial_name] = (f'{organism}/PPI_interface_mapped_to_exon.csv', ',')
    datasets['exons_to_domains_data_' + trivial_name] = (f'{organism}/final.csv', ',')
    datasets['gene_info_' + trivial_name] = (f'{organism}/gene_info.csv', ',')
    datasets['domain_gene_' + trivial_name] = (f'{organism}/gene_name2entrez_id.csv', ',')

    if os.path.isfile('domain/data/' + organism + '/Pfam-A.clans.tsv'):
        datasets['domain_domain_' + trivial_name] = (f'{organism}/Pfam-A.clans.tsv', '\t')


def load_datasets(datasets):
    print('=' * 80)
    print(f'Importing data sets from "{data_base_path}"\ninto {engine.name} database')
    print('=' * 80)
    # Dictionary of table_name : (file_name, separator)

    # --- Write dataframes to tables in database
    for table_name, data_file in datasets.items():
        print(f'Adding table {table_name}:')
        print(f'\tParsing file "{data_file[0]}"')
        data = pd.read_csv(path.join(data_base_path, data_file[0]), sep=data_file[1], low_memory=False)
        print(f'\tWriting {data.shape[0]} rows and {data.shape[1]} columns to database')

        if table_name.startswith('domain_gene'):
            # reduce data columns to only the ones we need
            data = data[['Gene name', 'Gene stable ID']]
            data.columns = ['gene_symbol', 'ensembl_id']
            data['id'] = data.index
            data.set_index('id', inplace=True)
            data.to_sql(table_name, engine, if_exists='replace', index=True)

        if table_name.startswith('domain_domain'):
            # reduce data columns to only the ones we need
            data = data[['PfamId', 'Symbol3', 'Description']]
            data.columns = ['pfam_id', 'symbol', 'description']
            data['id'] = data.index
            data.set_index('id', inplace=True)
            data.to_sql(table_name, engine, if_exists='replace', index=True)

        # Write to database directly using Pandas to_sql method
        else:
            data.to_sql(table_name, engine, if_exists='replace', index=False)
        print(f'\tDone')


class Command(BaseCommand):
    help = 'Import the datasets (CSVs) into the specified database'

    def handle(self, *args, **kwargs):
        load_datasets(datasets)


if __name__ == '__main__':
    load_datasets(datasets)
