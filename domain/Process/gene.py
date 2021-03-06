from domain.Process import process_data as pr
from domain.Process import exonstodomain as exd
from domain.Process import proteininfo as  info
import pandas as pd
from django.urls import reverse
from sqlalchemy import text

from django.conf import settings

from domain.Process import network_analysis as nt

# --- Get database connection aka 'SQLAlchemie engine'
engine = settings.DATABASE_ENGINE  

#load DIGGER Join Graph
DomainG=exd.load_obj("DomainG")

    
def TranscriptsID_to_table(transcripts,entrez='0'):
    if len(transcripts)>=1:
                #print('1111111111') 
                ID=[]
                name=[]
                pfams=[]
                missing_PPI=[]

                all_pfams=nt.g2d[entrez]
                #print(transcripts)
                for tr in transcripts :
                                           
                          query = """
                          SELECT * 
                          FROM exons_to_domains_data 
                          WHERE "Transcript stable ID"=:transcript_id 
                          """
                          tdata = pd.read_sql_query(sql=text(query), con=engine, params={'transcript_id': tr})
                          
                          # tdata=tdata.drop(columns=["Unnamed: 0"]).drop_duplicates()
                          
                          
                          
                          #df_filter = pr.data['Transcript stable ID'].isin([tr])
                          #tdata=pr.data[df_filter]
                          
                          
                          
                          #print(tdata)
                          if len(tdata)!=0  :
                              
                              tmp=pr.tranID_convert(tr)
                              if tmp==0: continue
                              n=tmp[0]
                              name.append(n)
                              ID.append(tr)
                              p=tdata["Pfam ID"].unique()
                              p = p[~pd.isnull(p)]
                              p=sorted(p)

                              # look for interesting isoforms with missing PPI using the join graph
                              missing=[ x for x in all_pfams if x not in p]
                              missing=[ entrez + '/' + x  for x in missing]
                              missing_interaction='-'
                              if any(DomainG.has_node(x) for x in missing):
                                      missing_interaction='&#9989;'
                              missing_PPI.append(missing_interaction)

                              # add hyperlink
                              p = [nt.link(x) for x in p]
                              pfams.append(' ; '.join(p))

                
                
                
                if ID!=[]:
                                
                          pd_isoforms=pd.DataFrame(list(zip(name, ID,pfams,missing_PPI)), columns =['Transcript name', 'Transcript ID','Pfam domains','Interacting domains are missing in the isoform'])
                          pd_isoforms['length'] = pd_isoforms['Pfam domains'].str.len()
                          pd_isoforms.sort_values('length', ascending=False, inplace=True)
                          pd_isoforms=pd_isoforms.drop(columns=['length'])
                          
                          h=reverse('home')+"ID/"
                          pd_isoforms["Link"]='<a href="'+h+pd_isoforms["Transcript ID"]+'">'+" (Visualize) "+'</a>'
    
                          pd.set_option('display.max_colwidth',1000)
                          
                          
                          pd_isoforms=pd_isoforms.to_html(**settings.TO_HTML_PARAMETERS)
              
            
    
    
    
    return pd_isoforms,n.split('-')[0]
    
    
    
    
    #changed
def input_gene(gene_ID):   
      #get a list of all transcripts of the selected gene
          pd_isoforms=[]
          n=''
         
          transcripts=pr.gene_to_all_transcripts(gene_ID)
          entrez=str(nt.ensembl_to_entrez(gene_ID))

          if len(transcripts)==0:
              return [],[]

          pd_isoforms,n=TranscriptsID_to_table(transcripts,entrez)
                 
         
          
          return pd_isoforms,n