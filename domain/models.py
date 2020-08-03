from django.db import models


class Gene(models.Model):
    ensembl_id = models.CharField(max_length=20)
    gene_symbol = models.CharField(max_length=20, db_index=True)

    def __str__(self):
        return self.gene_symbol


class Domain(models.Model):
    pfam_id = models.CharField(max_length=10, db_index=True, primary_key=True)
    symbol = models.CharField(max_length=10)
    description = models.CharField(max_length=150)

    def __str__(self):
        return f"Pfam: {self.pfam_id} - Symbol: {self.symbol}"
