# from domain.models import Gene
from pprint import pprint

from pybiomart import Dataset

dataset = Dataset(name='hsapiens_gene_ensembl',
                  host='http://www.ensembl.org')

mapping_df = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
                           filters={'chromosome_name': list(range(1,23)) + ['X', 'Y']})

print(mapping_df)
pprint(dataset.filters)

# for _, row in mapping_df.iterrows():
#     gene = Gene()
#     gene.ensembl_id = row['Gene stable ID']
#     gene.gene_symbol = row['Gene name']
#     gene.save()
