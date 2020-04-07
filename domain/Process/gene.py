from domain.Process import process_data as pr
from domain.Process import exonstodomain as exd 
from domain.Process import proteininfo as  info
import pandas as pd
    
    
    
    
    
    
def TranscriptsID_to_table(transcripts):
    if len(transcripts)>=1:
                #print('1111111111') 
                ID=[]
                name=[]
                pfams=[]
               
                for tr in transcripts :
                           
                          df_filter = pr.data['Transcript stable ID'].isin([tr])
                          tdata=pr.data[df_filter]
                          #print(tdata)
                          if len(tdata)!=0  :
                              ID.append(tr)
                              n=pr.tranID_convert(tr)[0]
                              name.append(n)
                              
                              p=tdata["Pfam ID"].unique()
                              p = p[~pd.isnull(p)]
                              
                              pfams.append('&emsp;<center>'+' ; '.join(p)+'</center>&emsp;')
                
                
                if ID!=[]:
                                
                          pd_isoforms=pd.DataFrame(list(zip(name, ID,pfams)), columns =['<center>Transcript name<center>', '<center>Transcript ID</center>','<center>Pfam domains<center>'])
                          pd_isoforms['length'] = pd_isoforms['<center>Pfam domains<center>'].str.len()
                          pd_isoforms.sort_values('length', ascending=False, inplace=True)
                          pd_isoforms=pd_isoforms.drop(columns=['length'])
                          
                          h="/ID/"
                          pd_isoforms["<center>Link</center>"]='<center>&emsp;'+'<a target="'+'_blank"href="'+h+pd_isoforms["<center>Transcript ID</center>"]+'">'+" (Visualize) "+'</a>'+'&emsp;</center>'
                          pd_isoforms["<center>Transcript ID</center>"]='<center>&emsp;'+pd_isoforms["<center>Transcript ID</center>"]+'&emsp;</center>'
                          pd.set_option('display.max_colwidth',1000)
                          
                          
                          pd_isoforms=pd_isoforms.to_html(escape=False, index=False)
              
            
    
    
    
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