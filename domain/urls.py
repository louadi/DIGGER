from django.urls import path
from django.views.generic import TemplateView

from domain import views, autocomplete

urlpatterns = [
    # --- Analysis Setup Pages ---
    # Isoform level
    path('isoform_level', views.isoform_level, name="isoform-level"),
    # Exon level
    path('exon_level', views.exon_level, name="exon-level"),
    # Network level
    path('network_analysis/', views.network, name="network-analysis"),
    path('network_analysis/example1', TemplateView.as_view(template_name='setup/network_example_1.html'),
         name="network-analysis-example-1"),
    path('network_analysis/example2', TemplateView.as_view(template_name='setup/network_example_2.html'),
         name="network-analysis-example-2"),
    path('network_analysis/example3', TemplateView.as_view(template_name='setup/network_example_3.html'),
         name="network-analysis-example-3"),

    # --- Selection Page
    # Interaction view page (select one transcript of gene)
    path('ID/gene/<str:gene_ID>/', views.gene, name="gene-overview"),

    # --- Detailed Result Pages ---
    # Network interaction view
    path('vis_network/job/<str:job>', views.Multi_proteins, name="network-visualization"),
    # Transcript view
    path('ID/<str:P_id>', views.transcript, name="transcript"),
    # Exon view
    path('ID/exon/<str:exon_ID>/', views.exon, name="exon"),

    # --- Autocomplete views ---
    path('gene_symbol-autocomplete/', autocomplete.gene_symbol_autocomplete, name='gene-symbol-autocomplete', ),

    # --- MISC ---
    # path('graph/', views.graph, name="graph"),
    # <<< Not used anymore?
    # Domain view
    path('graph/<str:Pfam_id>', views.display, name="Node_vis"),
    # Interaction view
    path('ID/<str:P_id>/InteractionView/<str:P2_id>', views.InteractionView, name="InteractionView"),
    # >>>
]
