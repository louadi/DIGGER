import pandas as pd
from django.core.management import BaseCommand

from django.conf import settings
from os import path

# --- Get database connection aka 'SQLAlchemie engine'
engine = settings.DATABASE_ENGINE

# --- Import the datasets into dataframes
data_base_path = path.join(settings.BASE_DIR, 'domain/data')


def load_datasets():
    datasets = {
        'ppi_data': pd.read_csv(path.join(data_base_path, 'PPI_interface_mapped_to_exon.csv')),
        'exons_to_domains_data': pd.read_csv(path.join(data_base_path, 'final.csv')),
        'gene_info': pd.read_csv(path.join(data_base_path, 'gene_info.csv')),
    }

    # --- Write dataframes to tables in database
    for table_name, data_df in datasets.items():
        data_df.to_sql(table_name, engine, if_exists='replace', index=False)


class Command(BaseCommand):
    help = 'Import the datasets (CSVs) into the specified database'

    def handle(self, *args, **kwargs):
        load_datasets()


if __name__ == '__main__':
    load_datasets()
