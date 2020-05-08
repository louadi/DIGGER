from django.urls import path
from . import views 



urlpatterns = [
    path('', views.home, name="home"),
    path('Exon_level', views.Exon_level, name="Exon_level"),
    path('isoform_level', views.isoform_level, name="isoform_level"),
    path('Network_analysis/', views.network, name="Network"),
    path('vis_network/job/<str:job>', views.Multi_proteins, name="Network_vis"),
    path('Network_analysis/example', views.example1, name="Network_example1"),
    path('Network_analysis/example2', views.example2, name="Network_example2"),
    path('about/', views.about, name="about_page"),
    path('graph/', views.graph, name="graph"),
    path('graph/<str:Pfam_id>', views.display, name="Node_vis"),
    path('ID/<str:P_id>', views.transcript, name="transcript"),
    path('ID/<str:P_id>/InteractionView/<str:P2_id>', views.InteractionView, name="InteractionView"),
    path('ID/gene/<str:gene_ID>/', views.gene, name="InteractionView"),
    path('ID/exon/<str:exon_ID>/', views.exon, name="Exon"),
     
    
    
]
