#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
import os
import mygene
import requests, sys

from sqlalchemy import text

from django_project import settings

server = "http://rest.ensembl.org"
cwd = os.getcwd()

# --- Get database connection aka 'SQLAlchemie engine'
engine = settings.DATABASE_ENGINE


#The data is generated by "map_exons_to_domains.py":
data=pd.read_csv( 'domain/data/final.csv')

#List of transcripts from genes from Biomart":
genes=pd.read_csv( 'domain/data/gene_to transcripts.csv')

#List of transcripts name and discriptions from Biomart":
gene_info=pd.read_csv( 'domain/data/gene_info.csv')



#load coverter
def load_obj(name):
    with open('domain/data/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)




#from Biomart":
gid2name=load_obj("gid2name")


def entrez_to_name_online(entrezID):
    mg = mygene.MyGeneInfo()
    out = mg.querymany(entrezID, scopes='entrezgene', fields='symbol', species='human')
    return out[0]['symbol']


def entrez_to_name(entrezID):

    return gid2name[entrezID]






def protein_to_transcript(Ensemble_prID):
    ext = "/lookup/id/"+Ensemble_prID +"?"
    info=req(ext)
    return info["Parent"]





def transcript(transcript_ID):
    # The input transcript ID
    # Outputs
        # Exons: all exons in the transcript
        # D: mapping coding exons to their respective protein domains 
        # uniaue domains in the protein

    query = """
            SELECT * 
            FROM exons_to_domains_data 
            WHERE "Transcript stable ID"=:transcript_id 
            ORDER BY "Exon rank in transcript"
            """
    tdata = pd.read_sql_query(sql=text(query), con=engine, params={'transcript_id': transcript_ID})

    # df_filter = data['Transcript stable ID'].isin([transcript_ID])
    # tdata=data[df_filter].sort_values(by=['Exon rank in transcript'])
    exons=tdata.drop(columns=["Pfam ID","Pfam start","Pfam end"]).drop_duplicates()
    D=tdata[tdata["Pfam ID"].notna()].drop_duplicates()
    p=D["Pfam ID"].unique()
    
    return exons,D,p




# cool function to search for all transcripts that contains a specific domain in the genome

def domain_search(Pfam_ID):
    return data[data['Pfam ID'].isin([Pfam_ID])].sort_values(by=['Exon rank in transcript'])









def req(ext):
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()
    return  r.json()




def gene_to_all_transcripts(gene_ID):



    df_filter = genes['Gene stable ID'].isin([gene_ID])
    tdata=genes[df_filter]
    
    
    
    #tdata=tdata["Transcript stable ID"].drop_duplicates()
    tdata=tdata["Transcript stable ID"].unique()
    
    
    
    return tdata





#not used now
def gene_to_all_transcripts_online(gene_ID):
    mg = mygene.MyGeneInfo()
    out = mg.querymany(gene_ID, scopes='ensembl.gene', fields='ensembl', species='human')[0]['ensembl']
    try:
        out=out['transcript']
    except TypeError:
        #print('error')
        for output in out:
            if output['gene']==gene_ID:
                out=output['transcript']
    return out
    
    
    
#not used    
def tranID_convert(Ensemble_transID):
    df_filter = gene_info['Transcript stable ID'].isin([Ensemble_transID])
    tdata=gene_info[df_filter]
    
 
            
    Ensemble_geneID=tdata["Gene stable ID"].unique()[0]
    tran_name=tdata["Transcript name"].unique()[0]
    gene_name=tran_name.split('-')[0]
    entrezID=tdata["NCBI gene ID"].unique()[0]
    if np.isnan(entrezID):
            mg = mygene.MyGeneInfo()
            entrezID = mg.querymany(Ensemble_geneID,  fields='entrezgene', species='human')[0]
            
            
            
    else: entrezID=str(int(tdata["NCBI gene ID"].unique()[0]))
    print('ID ========',entrezID)
    
    gene_description=tdata["Gene description"].unique()[0]
    
    return tran_name,gene_name,Ensemble_geneID,entrezID,gene_description
    
    
    

def tranID_convert_online(Ensemble_transID):
    ext = "/lookup/id/"+Ensemble_transID +"?"

    info=req(ext)
    Ensemble_geneID=info["Parent"]
    tran_name=info["display_name"]

    ext = "/xrefs/id/"+Ensemble_geneID +"?" 
    info=req(ext)
    
    
    for db in info :
        if db['db_display_name']=='NCBI gene': break
            
        
    entrezID=db["primary_id"]
    gene_name=db["display_id"]
    gene_description=db["description"]
    #print('ID ========',entrezID)
    return tran_name,gene_name,Ensemble_geneID,entrezID,gene_description





def unitpr_to_Ensemble(Protein_ID):
    Protein_ID=Protein_ID.split('%')[0]
    ext = "/xrefs/symbol/homo_sapiens/"+Protein_ID +"?"

    info=req(ext)
    for item in info:
        if item["type"]== "transcript":
            Ensemble_trID=item["id"]
            break
    for item in info:
        if item["type"]== "translation":
            Ensemble_prID=item["id"]
            break
        


    return Ensemble_trID,Ensemble_prID


def coordinate_to_exonID(gene_ID,s1,e1):
            
            exon_ID=[]
            trs=list(gene_to_all_transcripts(gene_ID))
            df_filter =data['Transcript stable ID'].isin(trs)
            df=data[df_filter]
            
            for index, row in df.iterrows()    :
            
                
                s2, e2= row["Genomic coding start"],row["Genomic coding end"]
                if pd.notna(s2) and   pd.notna(e2)  and is_overlapping(s1,e1,s2,e2):
                    exon_ID=row["Exon stable ID"]
                    return  exon_ID
    
            return  exon_ID

def is_overlapping(x1,x2,y1,y2):
    return max(x1,y1) <= min(x2,y2)