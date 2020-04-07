from domain.Process import exonstodomain as exd 
from domain.Process import process_data as pr
import pandas as pd
import numpy as np






PPI=exd.load_obj("PPI")
DDI=exd.load_obj("DDI")
g2d=exd.load_obj("g2d")


data = pd.read_csv( "domain/data/all_Proteins.csv")






def Construct_network(proteins_id, missing,job_ID):
    
      E=[]
      N=[]
      affected_nodes=[]
      DDI_nodes=[]
      p1=[]
      p2=[]
      p1_name=[]
      p2_name=[]
      DDIs=[]
      missing_DDIs=[]
      #Score = (retained DDI/All interactions)
      score=[]
      source=[]
      
      # Get the subgraph:
      H = PPI.subgraph(proteins_id)
      
      
      for e in H.edges():
          if e[0] not in N: N.append(e[0])
          if e[1] not in N: N.append(e[1])
      
          gene1=entrez_to_ensembl(e[0])
          gene2=entrez_to_ensembl(e[1])
          
          
          #check if gene is in missing dictionary
          #missing{} have all proteins with domain
          
          
          #Draw The PPI edge for every case
         
          #one of the gene have a domain
          if (gene1 in missing) and (gene2 in missing):
              inter=False
              
              if e[0] in g2d:
                    domains1=g2d[e[0]]
              else: 
                  domains1=[]
              if e[1] in g2d:
                  domains2=g2d[e[1]]
              else:
                  domains2=[]
                  
              DDIs_tmp=[]
              missing_DDIs_tmp=[]

              for d1 in domains1:
                  for d2 in domains2:
                      if DDI.has_edge(d1,d2):
                           inter=True
                              
                          
                           #Interacted domain is missing
                           if (d1 in missing[gene1]) or (d2 in missing[gene2]):
                                  E.append("{from: '"+e[0]+"', to: '"+e[1]+"', dashes:  true,title:' Lost "+d1+'-'+d2+"', color: 'red'},") 
                                  #print('yeaaaaah')
                                  print(d1,d2)
                                  print(gene1,gene2)
                                  affected_nodes.append(e[0])
                                  affected_nodes.append(e[1])
                                  
                                  missing_DDIs_tmp.append(d1+'-'+d2)
                           else:
                                  #Interaction retained
                                  E.append("{from: '"+e[0]+"', to: '"+e[1]+"',title:'"+d1+'-'+d2+"',  color: 'red'},") 
                                  DDI_nodes.append(e[0])
                                  DDI_nodes.append(e[1])  
                                  DDIs_tmp.append(d1+'-'+d2)
                                  
                                  
              if inter:
                   # inter=True: a Domain-Domain Interaction was found:
                   #for the table
                   p1.append(e[0])
                   p1_name.append(pr.entrez_to_name(e[0]))  
                   p2.append(e[1])
                   p2_name.append(pr.entrez_to_name(e[1]))
                   if DDIs_tmp!=[]: DDIs.append(' ; '.join(DDIs_tmp))
                   else: DDIs.append('-')
                   if missing_DDIs_tmp!=[]: missing_DDIs.append(' ; '.join(missing_DDIs_tmp))
                   else: missing_DDIs.append('-')
                   score.append(str(float("{0:.2f}".format(len(DDIs_tmp)/(len(DDIs_tmp)+len(missing_DDIs_tmp))))))         
                   source.append('PPI-DDI')           
                              
              if not inter:    
                   E.append("{from: '"+e[0]+"', to: '"+e[1]+"', title:'PPI',  color: CHOOSEN,       smooth: {type: 'continuous'}},") 
                   p1.append(e[0])
                   p1_name.append(pr.entrez_to_name(e[0]))  
                   p2.append(e[1])
                   p2_name.append(pr.entrez_to_name(e[1]))
                   DDIs.append('-')
                   missing_DDIs.append('-')
                   score.append('1.0')
                   source.append('PPI') 
          else : #Both genes have no domains
                   E.append("{from: '"+e[0]+"', to: '"+e[1]+"', title:'PPI evidence',  color: CHOOSEN,       smooth: {type: 'continuous'}},") 
                   E.append("{from: '"+e[0]+"', to: '"+e[1]+"', title:'PPI',  color: CHOOSEN,       smooth: {type: 'continuous'}},") 
                   p1.append(e[0])
                   p1_name.append(pr.entrez_to_name(e[0]))  
                   p2.append(e[1])
                   p2_name.append(pr.entrez_to_name(e[1]))
                   DDIs.append('-')
                   missing_DDIs.append('-')
                   score.append('1.0')
                   source.append('PPI')
      nodes=[]
      for n in N:
          ensembl=entrez_to_ensembl(n)
          if (n in affected_nodes) and (missing[ensembl]!=[]):
              nodes.append("{id: \""+n+"\",url:  '/ID/gene/"+ensembl+"' , color: CHOOSEN2, label:  '"+pr.entrez_to_name(n)+"'}, ")
              
          elif n in DDI_nodes:
              nodes.append("{id: \""+n+"\", url:  '/ID/gene/"+ensembl+"' , color: CHOOSEN3, label:  '"+pr.entrez_to_name(n)+"'},")
      
          else:    
              nodes.append("{id: \""+n+"\", url:  '/ID/gene/"+ensembl+"' , label:  '"+pr.entrez_to_name(n)+"'},")    
      
      
      # Table of interactions
      pd_interaction=pd.DataFrame( data=[p1,p1_name,p2,p2_name,DDIs,  missing_DDIs,score,source ], index=['Protein 1', 'Protein 1 name','Protein 2', 'Protein 2 name', 'Retained DDIs', 'Lost DDIs','Score','Source of the interaction'] ) 
      pd_interaction=pd_interaction.transpose()
      pd.set_option('display.max_colwidth',1000)
      
      
      pd_interaction.to_csv('domain/static/jobs/tables/'+job_ID+'.csv', index=False)
      
      
      pd_html=pd_interaction.drop(columns=['Protein 1', 'Protein 2'])
      
      
      #If the table is big, only show important Interactions
      if len(pd_html)>10:
          pd_html=pd_html[pd_html['Source of the interaction']=='PPI-DDI']
      
      #Max interaction to display in web page
      max_disp=150
      
      if len(pd_html)>max_disp:
            pd_html=pd_html[:max_disp]
      
      pd_html.sort_values(by=['Source of the interaction'], ascending=[True])
      
      
      pd_html['Protein 1 name']='<center>&emsp;'+pd_html['Protein 1 name']+'&emsp;</center>'
      pd_html['Protein 2 name']='<center>&emsp;'+pd_html['Protein 2 name']+'&emsp;</center>'
      pd_html['Retained DDIs']='<center>&emsp;'+pd_html['Retained DDIs']+'&emsp;</center>'
      pd_html['Lost DDIs']='<center>&emsp;'+pd_html['Lost DDIs']+'&emsp;</center>'
      pd_html['Score']='<center>&emsp;'+pd_html['Score']+'&emsp;</center>'
      pd_html['Source of the interaction']='<center>&emsp;'+pd_html['Source of the interaction']+'&emsp;</center>'
      
      
      pd_html=pd_html.rename(columns={
          "Protein 1 name": "<center>Protein 1 name</center>",
          "Protein 2 name": "<center>Protein 2 name</center>", 
          "Source of the interaction": "<center>&emsp;Source of the interaction&emsp;</center>", 
          "Score": "<center> Score</center>",
          "Retained DDIs": "<center>Retained Domain-Domain interactions</center>", 
          "Lost DDIs": "<center>Lost Domain-Domain interactions</center>",
          
          })
      
      
      
      # Convert Table to HTML
      pd_html=pd_html.to_html(escape=False, index=False)
      
      
      
      #Generate the network file (only edges and scores)
      pd_interaction=pd_interaction.rename(columns={
          "Score": "Weight",
          })
          
          
      pd_interaction['Weight']=pd_interaction['Weight'].astype(float)     
      pd_interaction.drop(columns=['Protein 1 name', 'Protein 2 name','Retained DDIs','Lost DDIs','Source of the interaction']).to_csv('domain/static/jobs/networks/'+job_ID+'.sif', index=False,sep='\t')
      
      
      
      return nodes,E,pd_interaction,pd_html







        
def analysis_input_isoforms(Inputs):
        print(Inputs)
        filtred=filter_proteins_list(Inputs)
        print('-----------filtred--------------')
        


        n=len(filtred)
        #Inputs are so many
        if n>2000:
              return False

        else:
                print('-----------yeaaaaaaah--------------')
                gene_domains={}
                #all genes and missing domains
                missing={}
                # all entrex ID of the genes
                proteins_id=[]

                for tr in filtred:
                          #print(tr)
                          gene=tr_to_gene(tr)
                          domains=tr_to_domain(tr)
                          
                          if len(domains)==0:
                              print(tr,domains)
                              gene_domains[gene]=[]

                          if len(domains)>0:
                              #Check for every domain
                              for domain in domains:

                                      if gene not in gene_domains:

                                          gene_domains[gene]=[domain]

                                      else :
                                          gene_domains[gene].append(domain)
                                          gene_domains[gene]=Remove(gene_domains[gene])

      
                for tr in filtred:
                    ensembl=tr_to_gene(tr)
                    entrez=str(ensembl_to_entrez(ensembl))
                    proteins_id.append(entrez)
                    if entrez in g2d:
                        domains=g2d[entrez]


                        #all domains are lost for this gene
                        if ensembl not in gene_domains :
                          missing[ensembl]=domains

                        #check if all domains are there   
                        else:
                          missing[ensembl]=[]
                          for d in domains:
                              if d not in gene_domains[ensembl]:
                                  missing[ensembl].append(d)


                proteins_id=Remove(proteins_id)       

                return proteins_id, missing,n


def analysis_input_genes(Inputs):
      
    protein_id=[]
    missing={}
    for i in Inputs:
        f=i.replace(" ", "")
        f=f.split(".")[0]
        
        if len(f)==15 and f[0:4]=='ENSG' :
        
              #check if gene is coding:
              if  len(data[data['Gene stable ID'].isin([f])])!=0:
              
                  print(f)
                  #check PPI status:
                  
                  entrez_id=str(ensembl_to_entrez(f))
                  print(entrez_id)
                  if PPI.has_node(entrez_id):
                      print('yeaaaaaah')
                      protein_id.append(entrez_id)
                      missing[f]=[]
    print(len(protein_id))
    if len(protein_id)>2000 or len(protein_id)<1:
            return False 
    else : return protein_id,missing,0
                  
            







def pr_to_tr(pr):
    df_filter = data['Protein stable ID'].isin([pr])
    try: return data[df_filter]['Transcript stable ID'].unique()[0]
    except IndexError:
        return False

def tr_to_gene(tr):
    df_filter = data['Transcript stable ID'].isin([tr])
    try: return data[df_filter]['Gene stable ID'].unique()[0]
    except IndexError:
        return False

def ensembl_to_entrez(gene):
    df_filter = data['Gene stable ID'].isin([gene])
    try: return data[df_filter]['NCBI gene ID'].unique()[0]
    except IndexError:
        return False
    
def entrez_to_ensembl(gene):
    df_filter = data['NCBI gene ID'].isin([gene])
    try: return data[df_filter]['Gene stable ID'].unique()[0]
    except IndexError:
        return False   
    
def tr_to_domain(tr):
    df_filter = data['Transcript stable ID'].isin([tr])
    tdata=data[df_filter]
    tdata=tdata[tdata["Pfam ID"].notna()].drop_duplicates()
    try: return tdata["Pfam ID"].unique()
    except IndexError:
        return False

#add to view
def check_PPI_status(tr):
    df_filter = data['Transcript stable ID'].isin([tr])
    return PPI.has_node(data[df_filter]['NCBI gene ID'].astype('str').unique()[0]) 


def tr_is_coding(tr):
    return  len(data[data['Transcript stable ID'].isin([tr])])!=0

def filter_proteins_list(List):
    filtred_list=[]
    for tr in List:
        #Correct spaces, dotsof version....
        ftr=tr=tr.replace(" ", "")
        ftr=ftr.split(".")[0]
        ftr=ftr.split("'\'")[0]
        ftr=ftr.split("+")[0]
        #print(ftr)
        #make sure Correct ID
        ftr=ftr[0:15]
        if ftr[0:3]=='ENS':
            
            #Check if protein coding
            
            #input is a protein and coverted successfully 
            if ftr[3]=='P' :
                
                tmp= pr_to_tr(ftr)
                
                # The protein is converted to a transcript
                if tmp!=False:
                        # Check PPI status
                        if check_PPI_status(tmp):
                            filtred_list.append(tmp)
                
            #check if transcript is a coding protein
            elif ftr[3]=='T'and tr_is_coding(ftr) and check_PPI_status(ftr):
               filtred_list.append(ftr)

                
    return filtred_list


def Remove(duplicate): 
    final_list = [] 
    for num in duplicate: 
        if num not in final_list: 
            final_list.append(num) 
    return final_list 







