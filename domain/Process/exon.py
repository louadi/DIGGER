import os

import pandas as pd
import numpy as np

from django.conf import settings
from domain.Process import process_data as pr
from domain.Process import exonstodomain as exd
from domain.Process import proteininfo as  info
from domain.Process import transcript as  tr
from domain.Process import gene as  g
from sqlalchemy import text
from django.urls import reverse

# --- Get database connection aka 'SQLAlchemie engine'
engine = settings.DATABASE_ENGINE


#PPI_old= pd.read_csv("domain/data/PPI_interface_mapped_to_exon.csv")
#tr_to_name_old = pd.read_csv( "domain/data/gene_info.csv")

# --- Create folder
# Global table path
table_path = os.path.join(settings.MEDIA_ROOT, 'table')
if not os.path.exists(table_path):
    os.makedirs(table_path)



def input_exon(exon_ID):
    
            query = """
            SELECT * 
            FROM exons_to_domains_data 
            WHERE "Exon stable ID"=:exon_ID 
            ORDER BY "Exon rank in transcript"
            """
            tb = pd.read_sql_query(sql=text(query), con=engine, params={'exon_ID': exon_ID})
            
            # tb=tb.drop(columns=["Unnamed: 0"]).drop_duplicates()
    
       
            #tb=pr.data[pr.data['Exon stable ID'].isin([exon_ID])]  
            
            transcripts=tb['Transcript stable ID'].unique()
            
            if len(transcripts) == 0:
            
              # ID not found in the database
              return True 
              
            else:
                
                    #get info about the gene
                    _,gene_name,Ensemble_geneID,entrezID,gene_description=pr.tranID_convert(transcripts[0])
                    
                    #print(gene_name,Ensemble_geneID,entrezID,gene_description)
                    
                    ## Table of Transcripts:
                    ## a function from gene page that give list of transcripts with links
                    ## The table is in HTML already formated.
                    
                    tb_transc,_=g.TranscriptsID_to_table(transcripts)
                    
                    
                    
                    
                    
                    
                    
                    #Table of domains
                    tb=tb[tb["Pfam ID"].notna()].drop_duplicates()
                    domains=tb["Pfam ID"].unique()
                    table_domains=[]
                    
                    if len(domains)==0:
                    
                          #No domains found for the exon
                          return transcripts,domains, gene_name,Ensemble_geneID,entrezID,tb_transc,table_domains,-1
                    else:     
                          #function to count number of interactions of the domain:
                          table_domains,_,_=exd.expand_table(tb,domains,entrezID)
                          
                          
                          h=reverse('home')+"graph/"
                          
                          table_domains['Interactions mediated by the domain']=table_domains['Interactions mediated by the domain'].astype(int)
                          
                          
                         
                          
                          
                          
                          
                          table_domains["Link to other Databases"]='<a href="http://pfam.xfam.org/family/'+table_domains['Pfam ID']+'  "target="_blank">Pfam  </a>   &nbsp;&nbsp;&nbsp; <a href="https://3did.irbbarcelona.org/dispatch.php?type=domain&value='+table_domains['Pfam ID']+'"target="_blank">3did  </a>   </h5 class> '
                         
                         
                         

                          
                          table_domains=table_domains.drop(columns=['Exon rank in transcript','Exon stable ID','CDS start','CDS end','Pfam start','Pfam end','Transcript stable ID', 'Chromosome/scaffold name',"Strand","Genomic coding start","Genomic coding end"])
                          table_domains=table_domains.drop_duplicates()

                          table_domains["Symbol"],table_domains["Summary"]=zip(*table_domains['Pfam ID'].map(pr.Domain_name))
                          table_domains=table_domains[['Pfam ID','Interactions mediated by the domain','Symbol','Summary','Link to other Databases']]

                          pd.set_option('display.max_colwidth',1000)
                          Total = table_domains['Interactions mediated by the domain'].sum()
                          
                          
                          
                          
                          
                          
                          
                          table_domains=table_domains.to_html(escape=False, index=False)
                          
                    
            
                          return transcripts,domains, gene_name,Ensemble_geneID,entrezID,tb_transc,table_domains,Total
                          
                          
                          
                    
def vis_exon(missing_domain,entrezID,gene_name,ExonID):

   
    #missing domain: are the domains coded by the exon
    missing_domain = [entrezID+'/'+ e for e in missing_domain]
    
    
    
   
   
    #only if the exon code for domains with known interactions
    
    if tr.PPI.has_node(entrezID):
        g = exd.nx.Graph()
        g.add_edges_from(tr.PPI.edges(entrezID))
        

        
    #search for the protein domains interactions:
    # before coming here need to check if the protein have known domains
    protein_domain=tr.g2d[entrezID]
    

    
    for domain in protein_domain:
        node=entrezID+'/'+domain
        
        

        #add domain interactions to the graph:       
        if tr.DomainG.has_node(node):
            g.add_edges_from(tr.DomainG.edges(node))
            
            
        edges=[]
    #list contains protein confirmed by both PPI and DDI (for visualization)
    protein_with_DDI=[]
    
    #link domains with their proteins
    for node in g.nodes():
            #a node is domain node:
            if len(node.split('/'))==2 :
                edges.append((node.split('/')[0],node))
                protein_with_DDI.append((node.split('/')[0]))
                
    g.add_edges_from(edges)
            
            
        
    protein_with_DDI = list(set(protein_with_DDI))
    nodes,edges,_=tr.vis_node_(g,entrezID,protein_with_DDI,gene_name,missing_domain,[])
    
    pd_interaction=tr.table_interaction(gene_name,entrezID,entrezID,g,protein_with_DDI,missing_domain)
    
    pd_interaction.insert(0,'Affected Protein','')
    pd_interaction["Affected Protein"]=gene_name
    
    pd_interaction=pd_interaction.rename(columns={
    "Protein name": "Partner Protein", })
    pd_interaction.to_csv(f'{table_path}/{ExonID}.csv', index=False,)
    
    
    
    
    #pd_interaction["Retained DDIs"]='<center>&emsp;'+pd_interaction["Retained DDIs"]+'&emsp;</center>'
   #pd_interaction["Lost DDIs"]='<center>&emsp;'+pd_interaction["Lost DDIs"]+'&emsp;</center>'
    
    pd_interaction['Residue evidence']='<center>-<center>'
    
    pd_interaction=pd_interaction.sort_values(by=['Percentage of lost domain-domain interactions'])
    h=reverse('home')+"ID/"+entrezID+'.'+ExonID+'/InteractionView/'
    pd_interaction["Percentage of lost domain-domain interactions"]=pd_interaction["Percentage of lost domain-domain interactions"].astype(int)
    pd_interaction["Score"]=(1-((pd_interaction["Percentage of lost domain-domain interactions"]/100)))
    
    #pd_interaction["Percentage of lost domain-domain interactions"]='<center>'+pd_interaction["Percentage of lost domain-domain interactions"].astype(int).astype(str)+' % '+'</center>'
    #pd_interaction["Affected Protein"]='<center>'+pd_interaction["Affected Protein"]+'</center>'
    #pd_interaction["Protein-protein interaction"]='<center>'+pd_interaction["Protein-protein interaction"]+'<a target="'+'_blank"href="'+h+pd_interaction["NCBI gene ID"]+'">'+" (Visualize) "+'</a>'+'</center>'

        #dont move it from here
    #pd_interaction["NCBI gene ID"]='<center>'+pd_interaction["NCBI gene ID"]+'</center>'
    
    
    
    #pd_interaction['Residue evidence']='<center>-<center>'
    
    pd.set_option('display.max_colwidth',1000)
    
    #pd_interaction=pd_interaction.to_html(escape=False, index=False)
    return nodes,edges,pd_interaction
    
    
    
def PPI_inter(exon_ID,gene_name):
    # --- Get tables from database
    query = """
            SELECT DISTINCT "Transcript stable ID_x", "u_ac_1", "Transcript stable ID_y", "u_ac_2"
            FROM ppi_data 
            WHERE "Exon stable ID_x"=:exon_id
            """
    p1 = pd.read_sql_query(sql=text(query), con=engine, params={'exon_id': exon_ID})

    query = """
            SELECT DISTINCT "Transcript stable ID_x", "u_ac_1", "Transcript stable ID_y", "u_ac_2"
            FROM ppi_data 
            WHERE "Exon stable ID_y"=:exon_id
            """
    p2 = pd.read_sql_query(sql=text(query), con=engine, params={'exon_id': exon_ID})

    # Compare the new and old dataframes
    #p1_old=PPI_old[PPI_old['Exon stable ID_x'] == exon_ID].drop(columns=['Exon stable ID_x', 'Exon stable ID_y']).drop_duplicates()
    #p2_old=PPI_old[PPI_old['Exon stable ID_y'] == exon_ID].drop(columns=['Exon stable ID_y', 'Exon stable ID_x']).drop_duplicates()
    #assert (np.array_equal(p1.values, p1_old.values))
    #assert (np.array_equal(p2.values, p2_old.values))


    p2=p2[['Transcript stable ID_y','u_ac_2','Transcript stable ID_x','u_ac_1']]
    p2=p2.rename(columns= {
        'Transcript stable ID_y':'Transcript stable ID_x',
        'u_ac_2':'u_ac_1',
        'u_ac_1':'u_ac_2',
        'Transcript stable ID_x':'Transcript stable ID_y',   
    })
      
      
    p1=p1.append(p2, ignore_index = True).drop_duplicates()
    
    n=len(p1)
    if n!=0:
   
    
            
        
            p1['Protein with selected exonic region']=gene_name
            
            #Get the partner protein name
            transcripts=p1['Transcript stable ID_y'].tolist()
            p1['Partner Protein']=tr_to_names(transcripts)
            
            
            #DIspaly table in HTML
            p1=p1.rename(columns= {
                
                'u_ac_2':'Uniprot ID of Protein 2',
                'u_ac_1':'Uniprot ID of Protein 1',
                
            })
            
            p1[['Protein with selected exonic region','Partner Protein','Uniprot ID of Protein 1','Uniprot ID of Protein 2']].drop_duplicates().to_csv(f'{table_path}/{exon_ID}_interface.csv', index=False,)
           

               
            
            p1=p1[['Protein with selected exonic region','Partner Protein']].drop_duplicates()
            
            #max display in web page
            max=25
            p1=p1[:max]
            #p1=p1[:10]
    pd.set_option('display.max_colwidth',1000)
    #p_html=p1.to_html(escape=False, index=False)
    
    
    return p1,n


def tr_to_names(list_tr):
    names = []
    for tr in list_tr:
        # SQL
        query = """
                SELECT "Transcript name" 
                FROM gene_info 
                WHERE "Transcript stable ID"=:transcript_id
                LIMIT 1
                """
        name = \
            pd.read_sql_query(sql=text(query), con=engine, params={'transcript_id': tr}).iloc[0, 0].split('-')[0]
        #name_old = tr_to_name_old[tr_to_name_old['Transcript stable ID'] == tr]['Transcript name'].tolist()[0].split('-')[0]

        #assert (np.array_equal(name, name_old))

        names.append(name)
    return names
