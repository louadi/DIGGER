from domain.Process import process_data as pr
from domain.Process import exonstodomain as exd 
from domain.Process import proteininfo as  info
from domain.Process import transcript as tr
import pandas as pd
from domain.Process import exon as ex 

PPI=exd.load_obj("PPI")
g2d=exd.load_obj("g2d")
DomainG=exd.load_obj("DomainG")




def int_view(P_id,P2_id):
    
    
    if P_id[0]=='E':
                  
                  domains,unique_domains,exons,text1,domainshtml,Text_nodes,text_edges,tran_name,gene_name,Ensemble_geneID,entrezID,gene_description,exons,droped1,droped2,trID=info.get_protein_info2(P_id)
                  
                                 
                  
                      
                  if PPI.has_node(entrezID):
                      g = exd.nx.Graph()
                      g.add_edges_from(PPI.edges(entrezID))
                      
                  else : 
                        print('no interactions')
                      
                  #search for the protein domains interactions:
                  # before coming here need to check if the protein have known domains
                  protein_domain=g2d[entrezID]
                  missing_domain=[]
                  
                  for domain in protein_domain:
                      node=entrezID+'/'+domain
                      
                      #check for missing domains in the transcript (for visualization):
                      if domain not in unique_domains:
                          missing_domain.append(node)
                             
                      #add domain interactions to the graph:       
                      if DomainG.has_node(node):
                          g.add_edges_from(DomainG.edges(node))
                  
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
                   
                  
                  
                  
                  
                  
                  _,DI_edges,lost_edges=tr.Interacted_domain(P2_id,g,entrezID,missing_domain)
              
                  
                  
                  nodes,edges,_=vis_interaction_(g,entrezID,protein_with_DDI,tran_name,missing_domain,P2_id)
              
                  
                  
                  if lost_edges==[]: lost_edges=''
                  if DI_edges==[]: DI_edges=''
              
                  
                  p_name=pr.entrez_to_name(P2_id)
                  pd_interaction=pd.DataFrame( data=[p_name,P2_id,' ; '.join(DI_edges),' ; '.join(lost_edges)], index=['Protein name','NCBI gene ID', 'Retained DDIs', 'Lost DDIs'] ) 
                  
                  pd_interaction=pd_interaction.transpose()
                  
                  
                  pd_interaction["Protein name"]='<center>'+pd_interaction["Protein name"]+'</center>'
                  pd_interaction["NCBI gene ID"]='<center>'+pd_interaction["NCBI gene ID"]+'</center>'
                  pd_interaction["Retained DDIs"]='<center>'+pd_interaction["Retained DDIs"]+'</center>'
                  pd_interaction["Lost DDIs"]='<center>'+pd_interaction["Lost DDIs"]+'</center>'
                  
                  
                  pd_interaction=pd_interaction.rename(columns={
                  
                  "Protein name": "<center>Protein name</center>", 
                  "NCBI gene ID": "<center>NCBI gene ID</center>", 
                  "Retained DDIs": "<center>Retained DDIs</center>", 
                  "Lost DDIs": "<center>Lost DDIs</center>"
                  
                  })
                  
                 
                  
              
                  pd.set_option('display.max_colwidth',1000)
                  pd_interaction=pd_interaction.to_html(escape=False, index=False)
                  
                  return nodes,edges,pd_interaction,tran_name,p_name
                  
                  
                                 
    else:
                 
      exon_ID=P_id.split('.')[1]        
      _,missing_domain, gene_name,Ensemble_geneID,entrezID,tb_transc,table_domains,number=ex.input_exon(exon_ID)
                
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
                     
                         
      _,DI_edges,lost_edges=tr.Interacted_domain(P2_id  ,g,entrezID,missing_domain)
            
                
                
      nodes,edges,_=vis_interaction_(g,entrezID,protein_with_DDI,gene_name,missing_domain,P2_id)         
      
      if lost_edges==[]: lost_edges=''
      if DI_edges==[]: DI_edges=''
      
      
      p_name=pr.entrez_to_name(P2_id)
      pd_interaction=pd.DataFrame( data=[p_name,P2_id,' ; '.join(DI_edges),' ; '.join(lost_edges)], index=['Protein name','NCBI gene ID', 'Retained DDIs', 'Lost DDIs'] ) 
      
      pd_interaction=pd_interaction.transpose()
      
      
      pd_interaction["Protein name"]='<center>'+pd_interaction["Protein name"]+'</center>'
      pd_interaction["NCBI gene ID"]='<center>'+pd_interaction["NCBI gene ID"]+'</center>'
      pd_interaction["Retained DDIs"]='<center>'+pd_interaction["Retained DDIs"]+'</center>'
      pd_interaction["Lost DDIs"]='<center>'+pd_interaction["Lost DDIs"]+'</center>'
      
      
      pd_interaction=pd_interaction.rename(columns={
      
      "Protein name": "<center>Protein name</center>", 
      "NCBI gene ID": "<center>NCBI gene ID</center>", 
      "Retained DDIs": "<center>Retained DDIs</center>", 
      "Lost DDIs": "<center>Lost DDIs</center>"
      
      })
      
      
      
      
      pd.set_option('display.max_colwidth',1000)
      pd_interaction=pd_interaction.to_html(escape=False, index=False)
      
      return nodes,edges,pd_interaction,gene_name,p_name











def vis_interaction_(g,entrezID,protein_with_DDI,tran_name,missing_domain,p2):

    N=[]
    E=[]
    for node in g.nodes():
    
        #filter other proteins:
        if node.split('/')[0]==entrezID or node.split('/')[0]==p2:
            #node of a protein:
            if node==p2 and len(node.split('/'))==1:
                    try:
                        label=pr.entrez_to_name(node)
                        N.append("{id: \""+node+"\", label:  \""+label+"\", group:  \""+tr.group_node(node,entrezID)+"\", physics:"+tr.physics(node,entrezID)+", source:  \""+tr.source_node(node,entrezID,protein_with_DDI)+"\", value:  \""+tr.value_node(node,entrezID)+"\"},")
                     
                    except KeyError:
                        print(node)
    
                        
            else:
                        label=tr.node_label(node,entrezID,tran_name)
                        N.append("{id: \""+node+"\", label:  \""+label+"\", group:  \""+tr.group_node(node,entrezID)+"\", physics:"+tr.physics(node,entrezID)+
                                 ", source:  \""+tr.source_node(node,entrezID,protein_with_DDI)+"\", value:  \""+tr.value_node(node,entrezID)+"\"},")

    for e in g.edges():
        gene1=e[0].split('/')[0]
        gene2=e[1].split('/')[0]
        
        if all(x in [entrezID,p2] for x in [gene1,gene2]):
            E.append("{from: \""+e[0]+"\", to: \""+e[1]+"\", dashes:  "+tr.edge_dashes(e,entrezID,missing_domain)[0]+","+tr.edge_option(e,entrezID)+"},")      

    return N,E,len(g)-1























