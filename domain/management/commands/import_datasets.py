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


def load_datasets():
    print('=' * 80)
    print(f'Importing data sets from "{data_base_path}"\ninto {engine.name} database')
    print('=' * 80)
    # Dictionary of table_name : (file_name, separator)
    datasets = {
        'ppi_data': ('PPI_interface_mapped_to_exon.csv', ','),
        'exons_to_domains_data': ('final.csv', ','),
        'gene_info': ('gene_info.csv', ','),
        Gene: ('gene_name2entrez_id.csv', ','),
        Domain: ('Pfam-A.clans.tsv', '\t'),

    }

    # --- Write dataframes to tables in database
    for table_name, data_file in datasets.items():
        print(f'Adding table {table_name}:')
        print(f'\tParsing file "{data_file[0]}"')
        data = pd.read_csv(path.join(data_base_path, data_file[0]), sep=data_file[1])
        print(f'\tWriting {data.shape[0]} rows and {data.shape[1]} columns to database')

        # Write to database using Django models
        if type(table_name) is django.db.models.base.ModelBase:
            table_name.objects.all().delete()

            # Gene
            if table_name is Gene:
                for _, row in tqdm(data.iterrows(), total=data.shape[0]):
                    entry = Gene()
                    entry.gene_symbol = row['Gene name']
                    entry.ensembl_id = row['Gene stable ID']
                    entry.save()

            # Domain
            elif table_name is Domain:
                for _, row in tqdm(data.iterrows(), total=data.shape[0]):   	
                    entry = Domain()
                    entry.pfam_id = row['PfamId']
                    entry.symbol = row['Symbol3']
                    entry.description = row['Description']
                    entry.save()

        # Write to database directly using Pandas to_sql method
        else:
            data.to_sql(table_name, engine, if_exists='replace', index=False)
        print(f'\tDone')


class Command(BaseCommand):
    help = 'Import the datasets (CSVs) into the specified database'

    def handle(self, *args, **kwargs):
        load_datasets()


if __name__ == '__main__':
    load_datasets()
