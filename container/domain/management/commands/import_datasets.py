import os
from os import path

import django
import pandas as pd
from django.conf import settings
from django.core.management import BaseCommand
from sqlalchemy import text
from tqdm import tqdm

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


# Columns each table family is filtered on in domain/Process/*.py. Indexing them
# turns the per-lookup full table scans into index seeks. domain_gene_* is never
# queried by the app, so it deliberately gets no index.
INDEX_COLUMNS = {
    'ppi_data': ['Exon stable ID_x', 'Exon stable ID_y',
                 'Transcript stable ID_x', 'Transcript stable ID_y'],
    'exons_to_domains_data': ['Transcript stable ID', 'Exon stable ID', 'Pfam ID'],
    'gene_info': ['Transcript stable ID', 'Gene stable ID'],
    'domain_domain': ['pfam_id'],
}


def create_indexes(table_name):
    """Create indexes on the columns `table_name`'s family is filtered on.

    Uses CREATE INDEX IF NOT EXISTS (supported by both Postgres and SQLite) so it is
    safe to re-run. Identifiers are double-quoted because several column names contain
    spaces (e.g. "Transcript stable ID").
    """
    for prefix, columns in INDEX_COLUMNS.items():
        if table_name.startswith(prefix + '_'):
            with engine.begin() as conn:
                for column in columns:
                    index_name = f'ix_{table_name}_' + column.replace(' ', '_')
                    conn.execute(text(
                        f'CREATE INDEX IF NOT EXISTS "{index_name}" '
                        f'ON "{table_name}" ("{column}")'
                    ))
            print(f'\tIndexed columns: {", ".join(columns)}')
            break


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

        # Recreate indexes on the columns the app filters on. Tables are dropped and
        # recreated on every import (if_exists='replace'), so their indexes must be
        # rebuilt each time too.
        create_indexes(table_name)
        print(f'\tDone')


class Command(BaseCommand):
    help = 'Import the datasets (CSVs) into the specified database'

    def handle(self, *args, **kwargs):
        load_datasets(datasets)


if __name__ == '__main__':
    load_datasets(datasets)
