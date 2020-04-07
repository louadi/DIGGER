from domain.Process import process_data as pr
from domain.Process import exonstodomain as exd 
from domain.Process import proteininfo as  info
import pandas as pd


PPI=exd.load_obj("PPI")
g2d=exd.load_obj("g2d")
DomainG=exd.load_obj("DomainG")




def Protein_view(P_id):
    
    domains,unique_domains,exons,text1,domainshtml,Text_nodes,text_edges,tran_name,gene_name,Ensemble_geneID,entrezID,gene_description,exons,droped1,droped2,trID,p=info.get_protein_info(P_id)
    

    
        
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
    nodes,edges,_=vis_node_(g,entrezID,protein_with_DDI,tran_name,missing_domain)
    
    df_missed=[]
    if missing_domain!=[]:
            n=[]
            l=[]
            dom=[]
            h="/graph/"
            for d in missing_domain:
                s=d.split('/')
                dom.append(s[1])
                if DomainG.has_node(d):
                    n.append(len(DomainG[d]))
                    
                    
                    l.append('<a target="'+'_blank" href="'+h+entrezID+"."+s[1]+'">'+gene_name+'-'+s[1]+'</a>')
                else: 
                    n.append(0)
                    l.append('')
                    
                  
                  
            dfff = pd.DataFrame(list(zip(dom,n,l)), columns=['Pfam ID','Pfam known interactions','Visualization of the domain interactions'])
            pd.set_option('display.max_colwidth',1000)
            
            
            
            dfff["Pfam known interactions"]='<center>'+dfff["Pfam known interactions"].astype(str)+'</center>'
            dfff["Visualization of the domain interactions"]='<center>'+dfff["Visualization of the domain interactions"]+'</center>'
            df_missed=dfff.to_html(escape=False, index=False)
    
            
    
    
    
    #interactionView:
    
    if len(protein_with_DDI)>1:
          pd_interaction=table_interaction(tran_name,trID,entrezID,g,protein_with_DDI,missing_domain)
              
              
          pd_interaction.insert(0,'Selected Protein variant','')
          pd_interaction["Selected Protein variant"]=tran_name
          pd_interaction.to_csv('domain/static/table 2/'+trID+'.csv', index=False,)
          
    
          pd_interaction["Retained DDIs"]='<center>&emsp;'+pd_interaction["Retained DDIs"]+'&emsp;</center>'
          pd_interaction["Lost DDIs"]='<center>&emsp;'+pd_interaction["Lost DDIs"]+'&emsp;</center>'
          
          
          
          pd_interaction=pd_interaction.sort_values(by=['Percentage of lost domain-domain interactions'])
          h="/ID/"+trID+'/InteractionView/'
          
          pd_interaction["Percentage of lost domain-domain interactions"]='<center>'+pd_interaction["Percentage of lost domain-domain interactions"].astype(int).astype(str)+' % '+'</center>'
          
          pd_interaction["Protein-protein interaction"]='<center>'+pd_interaction["Protein-protein interaction"]+'<a target="'+'_blank"href="'+h+pd_interaction["NCBI gene ID"]+'">'+" (Visualize) "+'</a>'+'</center>'
          
          
          #dont move it from here
          pd_interaction["NCBI gene ID"]='<center>'+pd_interaction["NCBI gene ID"]+'</center>'
          pd_interaction["Selected Protein variant"]='<center>'+pd_interaction["Selected Protein variant"]+'</center>'
          pd_interaction["Protein name"]='<center>'+pd_interaction["Protein name"]+'</center>'
          
          
          pd_interaction=pd_interaction.rename(columns={
          "Selected Protein variant": "<center>Selected Protein variant</center>",
          "Protein name": "<center>Partner Protein </center>", 
          "NCBI gene ID": "<center>NCBI gene ID</center>", 
          "Percentage of lost domain-domain interactions": "<center> % of lost DDIs</center>",
          "Retained DDIs": "<center>&emsp;Retained Domain-Domain interactions</center>", 
          "Lost DDIs": "<center>Lost Domain-Domain interactions</center>",
          "Protein-protein interaction": "<center>Protein-protein interaction</center>"
          })
          
          pd.set_option('display.max_colwidth',1000)
          pd_interaction=pd_interaction.to_html(escape=False, index=False)
          
    else: pd_interaction=[]
    
    
    
    
    #get a list of all Isoforms of the selected transcript
    gene_ID=pr.tranID_convert(trID)[2]
    transcripts=pr.gene_to_all_transcripts(gene_ID)
    
    if len(transcripts)>=1:
  
          ID=[]
          name=[]
          pfams=[]
          
          for tr in transcripts :
              
              df_filter = pr.data['Transcript stable ID'].isin([tr])
              tdata=pr.data[df_filter]
              
              if len(tdata)!=0  :
                  ID.append(tr)
                  name.append(pr.tranID_convert(tr)[0])
                  
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
            
            
            
            pd_isoforms=pd_isoforms.to_html(escape=False, index=False)
            
    else: pd_isoforms=[]
    
    
    
    
    return nodes,edges,_,domains,unique_domains,exons,text1,domainshtml,Text_nodes,text_edges,tran_name,gene_name,Ensemble_geneID,entrezID,gene_description,exons,droped1,droped2,trID,p,df_missed,pd_interaction,pd_isoforms









def Interacted_domain(p,g,entrezID,missing_domain):
        #Search for interacted domains
        edges=[]
        DDI_edges=[]
        lost_edges=[]
        for e in g.edges():
            gene1=e[0].split('/')[0]
            gene2=e[1].split('/')[0]
            if p in [gene1,gene2]:
                edges.append(e)
                if edge_dashes(e,entrezID,missing_domain)[0]!='false':
                    lost_edges.append(e[0].split('/')[1]+'-'+e[1].split('/')[1])

                elif len(e[0].split('/'))==2 and len(e[1].split('/'))==2: 
                    DDI_edges.append(e[0].split('/')[1]+'-'+e[1].split('/')[1])
        
        return edges,DDI_edges,lost_edges     

def table_interaction(tran_name,trID,entrezID,g,protein_with_DDI,missing_domain):
    Interactions=[]
    IDs=[]
    DDIs=[]
    lost_DDIs=[]
    perc=[]
    status=[]
    for protein in protein_with_DDI:
        #select first protein
        if protein !=entrezID:
            #Search for interacted domains
            edges,DDI_edges,lost_edges=Interacted_domain(protein,g,entrezID,missing_domain)
            
            IDs.append(protein)
            p=len(lost_edges)/(len(DDI_edges)+len(lost_edges))
            p=float("{0:.2f}".format(p))
            perc.append(p*100)
            
            st='Affected'
            if p==0: st='Retained'
            if p==1: st='Lost'
            
            status.append(st)
            
            Interactions.append(pr.entrez_to_name(protein))
            DDIs.append(' ; '.join(DDI_edges))
            lost_DDIs.append(' ; '.join(lost_edges))
    return  pd.DataFrame(list(zip(Interactions, IDs,DDIs,lost_DDIs,perc,status)), columns =['Protein name', 'NCBI gene ID','Retained DDIs','Lost DDIs','Percentage of lost domain-domain interactions','Protein-protein interaction']) 





























def vis_node_(g,entrezID,protein_with_DDI,tran_name,missing_domain):

    N=[]
    E=[]
    for node in g.nodes():
        #node of a protein:
        if node!=entrezID and len(node.split('/'))==1:
                try:
                    label=pr.entrez_to_name(node)
                    N.append("{id: \""+node+"\", label:  \""+label+"\", group:  \""+group_node(node,entrezID)+"\", physics:"+physics(node,entrezID)+
                             ", source:  \""+source_node(node,entrezID,protein_with_DDI)+"\", value:  \""+value_node(node,entrezID)+"\"},")
                except KeyError:
                    print(' ')

                    
        else:
                    label=node_label(node,entrezID,tran_name)
                    N.append("{id: \""+node+"\", label:  \""+label+"\", group:  \""+group_node(node,entrezID)+"\", physics:"+physics(node,entrezID)+
                             ", source:  \""+source_node(node,entrezID,protein_with_DDI)+"\", value:  \""+value_node(node,entrezID)+"\"},")

    for e in g.edges():
         E.append("{from: \""+e[0]+"\", to: \""+e[1]+"\", dashes:  "+edge_dashes(e,entrezID,missing_domain)[0]+","+edge_option(e,entrezID)+"},")      

    return N,E,len(g)-1























def edge_option(edge,entrezID):
    e1=edge[0].split('/')
    e2=edge[1].split('/')
    #protein to protein
    if len(e1)==1 and len(e2)==1:
        option='length: PR_LENGTH, width: WIDTH_SCALE * 4'
    #domain to domain
    elif len(e1)==2 and len(e2)==2:
        option='length: PR_DM, color: GREEN, width: WIDTH_SCALE * 2'
        
    #domain to protein    
    else    :
        option='length: LENGTH_domain, color: YELLOW, width: WIDTH_SCALE * 2'
    return option    

def edge_dashes(edge,entrezID,missing_domain):
    v="false"
    lost_int=[]
    #one of the nodes are missing in this transcript:
    if any(x in missing_domain for x in edge):
      
        v='[2, 2, 10, 10] '
        lost_int.append(edge)
    return v,lost_int










# node group
def group_node(node,entrezID):
    group="protein"
    if len(node.split('/'))==2 : 
        group="domain"
    
    return group

# Seperate node by their source: DDI or PPI
# For visualization
# all domains node will be in mode 2 in additon to proteins with interacted domains
def source_node(node,entrezID,protein_with_DDI):
    source="PPI"
    if len(node.split('/'))==2 : 
        source="DDI"
    elif node in protein_with_DDI :
        source="DDI"
    return source


def physics(node,entrezID):
    p='true'
    if node==entrezID: p='false'
    return p


def value_node(node,entrezID,):
    v='2'
    if group_node(node,entrezID)!="domain"  :v='4'
    return v

def node_label(node,entrezID,tran_name):
    
    #node of selected transcript
    if node==entrezID:  label=tran_name
        
    #node is a domain
    elif  len(node.split('/'))==2 :
        label=node.split('/')[1]
    
    else :
        label=node
    return label