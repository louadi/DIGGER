from django.urls import path
from django.views.generic import TemplateView

from domain import views, autocomplete

urlpatterns = [
    path('', views.home, name="home"),
    path('Exon_level', views.Exon_level, name="Exon_level"),
    path('isoform_level', views.isoform_level, name="isoform_level"),
    path('Network_analysis/', views.network, name="Network"),
    path('vis_network/job/<str:job>', views.Multi_proteins, name="Network_vis"),
    path('Network_analysis/example1', TemplateView.as_view(template_name='domain/Network_example_1.html'), name="Network_example1"),
    path('Network_analysis/example2', TemplateView.as_view(template_name='domain/Network_example_2.html'), name="Network_example2"),
    path('Network_analysis/example3', TemplateView.as_view(template_name='domain/Network_example_3.html'), name="Network_example3"),
    path('about/', views.about, name="about_page"),
    path('documentation/', views.doc, name="doc_page"),
    path('download/', views.download, name="download_page"),
    path('graph/', views.graph, name="graph"),
    path('graph/<str:Pfam_id>', views.display, name="Node_vis"),
    path('ID/<str:P_id>', views.transcript, name="transcript"),
    path('ID/<str:P_id>/InteractionView/<str:P2_id>', views.InteractionView, name="InteractionView"),
    path('ID/gene/<str:gene_ID>/', views.gene, name="InteractionView"),
    path('ID/exon/<str:exon_ID>/', views.exon, name="Exon"),

    # Autocomplete views
    path('gene_symbol-autocomplete/', autocomplete.gene_symbol_autocomplete, name='gene_symbol-autocomplete',)
]
