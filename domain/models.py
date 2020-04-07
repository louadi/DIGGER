from django.db import models
from .Process import process_data as pr
from .Process import exonstodomain as exd 
from .Process import proteininfo as  info
# Create your models here.
import networkx as nx
import numpy as np

import pickle


def get_protein_info(ID):

  #Get data about the input transcript or the protein 
  domains,exons,text1,domains.to_html,Text_nodes,text_edges,tran_name,gene_name,Ensemble_geneID,entrezID,gene_description,exons,trID=info.ID_mapper(ID)   

  #save Image of protein Structure
  Protein_structure(ID,exons,domains)

  return domains,exons,text1,domains.to_html,Text_nodes,text_edges,tran_name,gene_name,Ensemble_geneID,entrezID,gene_description,exons,trID







def Protein_structure(ID,exons,domains):
    #save Image of protein Structure
    features1,features2,fend=Visualize_transciript(exons,domains)
    
    fig, (ax1, ax2) = plt.subplots(
        2, 1,figsize=(22, 3.5))
    
    record = GraphicRecord(sequence_length=fend, features=features1,)
    record.plot(ax=ax1,figure_width=23,with_ruler=False)
    
    record = GraphicRecord(sequence_length=fend, features=features2,)
    p=record.plot(ax=ax2,figure_width=23,with_ruler=True,annotate_inline=True)
    
    ax1.title.set_text('Coding Exons')
    ax1.title.set_position([.5, -0.4])
    ax2.title.set_text('Domains')
    ax2.title.set_position([.5, -0.4])
    p.figure.savefig('domain/static/images/transcript'+trID, bbox_inches='tight')


    