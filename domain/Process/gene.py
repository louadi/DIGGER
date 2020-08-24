from domain.Process import process_data as pr
from domain.Process import exonstodomain as exd 
from domain.Process import proteininfo as  info
import pandas as pd
from django.urls import reverse
from sqlalchemy import text

from django.conf import settings
# --- Get database connection aka 'SQLAlchemie engine'
engine = settings.DATABASE_ENGINE  
    
    
    
    
def TranscriptsID_to_table(transcripts):
    if len(transcripts)>=1:
                #print('1111111111') 
                ID=[]
                name=[]
                pfams=[]
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
                              pfams.append(' ; '.join(p))
                
                
                
                if ID!=[]:
                                
                          pd_isoforms=pd.DataFrame(list(zip(name, ID,pfams)), columns =['Transcript name', 'Transcript ID','Pfam domains'])
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
          

          if len(transcripts)==0:
              return [],[]

          pd_isoforms,n=TranscriptsID_to_table(transcripts)
                 
         
          
          return pd_isoforms,n