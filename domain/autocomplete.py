from django.http import JsonResponse

from domain.models import Gene
from django.views.decorators.http import require_http_methods


# Derived from here https://stackoverflow.com/questions/32465052/using-typeahead-js-in-django-project
@require_http_methods(['GET'])
def gene_symbol_autocomplete(request):
    q = request.GET.get('q')
    data = {}

    if q:
        qs = Gene.objects.filter(gene_symbol__istartswith=q)[:20]
        data = [{'gene_symbol': gene.gene_symbol, 'ensembl_id': gene.ensembl_id} for gene in qs]

    return JsonResponse(data, safe=False)
