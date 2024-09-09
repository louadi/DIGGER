from django.db import models
from django.utils import timezone


class Gene(models.Model):
    ensembl_id = models.CharField(max_length=20)
    gene_symbol = models.CharField(max_length=20, db_index=True)

    def __str__(self):
        return self.gene_symbol


class Domain(models.Model):
    pfam_id = models.CharField(max_length=10, db_index=True, primary_key=True)
    symbol = models.CharField(max_length=20)
    description = models.CharField(max_length=150)

    #def __str__(self):
    #   return f"Pfam: {self.pfam_id} - Symbol: {self.symbol}"

class NeaseSaveLocationMapping(models.Model):
    run_id = models.CharField(max_length=36, primary_key=True, db_index=True)
    saved_for_days = models.IntegerField()
    date_of_creation = models.DateTimeField(auto_now_add=True)
    name = models.CharField(max_length=255, default='')

    # Method to query the database for the run_id and return the saved_for_days
    @staticmethod
    def get_saved_for_days(run_id):
        return str(NeaseSaveLocationMapping.objects.get(run_id=run_id).saved_for_days)

    def get_number_of_saved_for_days(self):
        return str(self.saved_for_days)

    # Calculate how many days are left until deletion
    def days_left(self):
        return self.saved_for_days - (timezone.now() - self.date_of_creation).days