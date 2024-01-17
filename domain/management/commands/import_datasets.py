import os
from os import path

import django
import pandas as pd
from django.conf import settings
from django.core.management import BaseCommand
from tqdm import tqdm

from domain.models import Gene, Domain, GeneMouse

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
    if trivial_name == 'human':
        datasets[Gene] = (f'{organism}/gene_name2entrez_id.csv', ',')
        datasets[Domain] = (f'{organism}/Pfam-A.clans.tsv', '\t')
    elif trivial_name == 'mouse':
        datasets[GeneMouse] = (f'{organism}/gene_name2entrez_id.csv', ',')

# datasets = {
#     'ppi_data_human': ('PPI_interface_mapped_to_exon_human.csv', ','),
#     'exons_to_domains_data_human': ('final_human.csv', ','),
#     'gene_info_human': ('gene_info_human.csv', ','),
#     Gene: ('gene_name2entrez_id_human.csv', ','),
#     Domain: ('Pfam-A.clans.tsv', '\t'),
#
#     # --- Adding mouse data to database
#     GeneMouse: ('gene_name2entrez_id_mouse.csv', ','),
#     'ppi_data_mouse': ('PPI_interface_mapped_to_exon_mouse.csv', ','),
#     'exons_to_domains_data_mouse': ('final_mouse.csv', ','),
#     'gene_info_mouse': ('gene_info_mouse.csv', ','),
# }


def load_datasets(datasets):
    print('=' * 80)
    print(f'Importing data sets from "{data_base_path}"\ninto {engine.name} database')
    print('=' * 80)
    # Dictionary of table_name : (file_name, separator)

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

            # Mouse gene
            elif table_name is GeneMouse:
                for _, row in tqdm(data.iterrows(), total=data.shape[0]):
                    entry = GeneMouse()
                    entry.gene_symbol = row['Gene name']
                    entry.ensembl_id = row['Gene stable ID']
                    entry.save()

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
