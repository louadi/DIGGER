from django.contrib import admin
from django.urls import path, re_path
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

    # NEASE analysis
    path('nease_analysis', views.setup_nease, name="nease-analysis"),
    path('nease_analysis/extra_functions', views.nease_extra_functions, name="nease-image-test"),

    # --- Selection Page
    # Interaction view page (select one transcript of gene)
    path('ID/gene/<str:organism>/<str:gene_ID>/', views.gene, name="gene-overview"),
    path('ID/gene/<str:organism>/multiple/<str:inputs>/', views.multiple_queries, name="multiple-queries"),

    #exon selection page:
    



    # --- Detailed Result Pages ---
    # Network interaction view
    path('vis_network/job/<str:organism>/<str:job>', views.Multi_proteins, name="network-visualization"),
    # Transcript view
    path('ID/<str:organism>/<str:P_id>', views.transcript, name="transcript"),
    # Exon view
    path('ID/exon/<str:organism>/<str:exon_ID>/', views.exon, name="exon"),

    # --- Autocomplete views ---
    path('gene_symbol-autocomplete/', autocomplete.gene_symbol_autocomplete, name='gene-symbol-autocomplete', ),

    # --- MISC ---
    re_path(r'^.*?/get_organisms/', views.get_organisms, name='get_organisms'),
    path('get_organisms/', views.get_organisms, name='get_organisms'),
    # path('graph/', views.graph, name="graph"),
    # <<< Not used anymore?
    # Domain view
    # path('graph/<str:Pfam_id>', views.display, name="Node_vis"),
    # Interaction view
    # path('ID/<str:P_id>/InteractionView/<str:P2_id>', views.InteractionView, name="InteractionView"),
    # >>>
    path("admin/", admin.site.urls),
]
