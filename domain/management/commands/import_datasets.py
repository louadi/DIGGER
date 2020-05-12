import pandas as pd
from django.core.management import BaseCommand

from django.conf import settings
from os import path

# --- Get database connection aka 'SQLAlchemie engine'
engine = settings.DATABASE_ENGINE

# --- Import the datasets into dataframes
data_base_path = path.join(settings.PROJECT_ROOT, 'domain/data')


def load_datasets():
    print('='*80)
    print(f'Importing data sets from "{data_base_path}"\ninto {engine.name} database')
    print('='*80)
    datasets = {
        'ppi_data': path.join(data_base_path, 'PPI_interface_mapped_to_exon.csv'),
        'exons_to_domains_data': path.join(data_base_path, 'final.csv'),
        'gene_info': path.join(data_base_path, 'gene_info.csv'),
    }

    # --- Write dataframes to tables in database
    for table_name, data_path in datasets.items():
        print(f'Adding table {table_name}:')
        data = pd.read_csv(data_path)
        print(f'\tContains {data.shape[0]} rows and {data.shape[1]} columns')
        data.to_sql(table_name, engine, if_exists='replace', index=False)
        print(f'\tDone')


class Command(BaseCommand):
    help = 'Import the datasets (CSVs) into the specified database'

    def handle(self, *args, **kwargs):
        load_datasets()


if __name__ == '__main__':
    load_datasets()
