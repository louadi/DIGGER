from django.db import models


class Gene(models.Model):
    ensembl_id = models.CharField(max_length=20)
    gene_symbol = models.CharField(max_length=20, db_index=True)

    def __str__(self):
        return self.gene_symbol
