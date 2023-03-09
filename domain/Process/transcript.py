import os

from domain.Process import process_data as pr
from domain.Process import exonstodomain as exd 
from domain.Process import proteininfo as  info
from domain.Process import network_analysis as nt
import pandas as pd
from django.urls import reverse


from sqlalchemy import text

from django.conf import settings
# --- Get database connection aka 'SQLAlchemie engine'
engine = settings.DATABASE_ENGINE

# --- Create folder
# Global table path
table_path_2 = os.path.join(settings.MEDIA_ROOT, 'table 2')
if not os.path.exists(table_path_2):
    os.makedirs(table_path_2)

PPI_human=exd.load_obj("PPI_human")
g2d_human=exd.load_obj("g2d_human")

PPI_mouse=exd.load_obj("PPI_mouse")
g2d_mouse=exd.load_obj("g2d_mouse")

#Join PPI-DDI network
DomainG_human=exd.load_obj("DomainG_human")
DomainG_mouse=exd.load_obj("DomainG_mouse")


# PPI network confirmed by residue evidence
#PPI_3D=exd.load_obj("Residue")


def Protein_view(P_id,organism):
    
    i=info.get_protein_info(P_id,organism)
    if i==0: return 0
    domains,unique_domains,exons,text1,domainshtml,Text_nodes,text_edges,tran_name,gene_name,Ensemble_geneID,entrezID,gene_description,exons,droped1,droped2,trID,p,co_partners=i
    
    if organism == "human":
        PPI = PPI_human
        g2d = g2d_human
        DomainG = DomainG_human
    elif organism == "mouse":
        PPI = PPI_mouse
        g2d = g2d_mouse
        DomainG = DomainG_mouse

    if PPI.has_node(entrezID):
        g = exd.nx.Graph()
        g.add_edges_from(PPI.edges(entrezID))

    else :
          print('no interactions')
          return 1
        
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
    nodes,edges,_=vis_node_(g,entrezID,protein_with_DDI,tran_name,missing_domain,co_partners,organism)
    
    df_missed=[]
    if missing_domain!=[]:
            n=[]
            l=[]
            dom=[]
            
            for d in missing_domain:
                s=d.split('/')
                dom.append(s[1])
                l.append(' <a href="http://pfam.xfam.org/family/'+s[1]+'  "target="_blank">Pfam  </a>   &nbsp;&nbsp;&nbsp; <a href="https://3did.irbbarcelona.org/dispatch.php?type=domain&value='+s[1]+'"target="_blank">3did  </a>      </h5 class> ')
                if DomainG.has_node(d):
                    n.append(len(DomainG[d]))
                    
                    
                    
                    
                 
                    
                    
                    
                else: 
                    n.append(0)
                
                
               
                    
                  
                  
            dfff = pd.DataFrame(list(zip(dom,n,l)), columns=['Pfam ID','Interactions mediated by the domain','Link to other databases'])
            pd.set_option('display.max_colwidth',1000)
            dfff["Symbol"],dfff["Summary"]=zip(*dfff['Pfam ID'].map(pr.Domain_name))
            dfff=dfff[['Pfam ID','Interactions mediated by the domain','Symbol','Summary','Link to other databases']]
            df_missed=dfff

    
            
    
    
    
    #interactionView:
    
    if len(protein_with_DDI)>1:
          pd_interaction=table_interaction(tran_name,trID,entrezID,g,protein_with_DDI,missing_domain,organism)

          pd_interaction["Score"] = (
                      1 - ((pd_interaction["Percentage of lost domain-domain interactions"].astype(int) / 100)))
          
          pd_interaction.insert(0,'Selected Protein variant','')
          pd_interaction["Selected Protein variant"]=tran_name
          pd_interaction.drop(columns=['Lost DDIs','Retained DDIs']).to_csv(f'{table_path_2}/{trID}.csv', index=False,)

          pd_interaction=pd_interaction.drop(columns=['retained DDIs','missing DDIs'])
    
          pd_interaction["Retained DDIs"]='&emsp;'+pd_interaction["Retained DDIs"]+'&emsp;'
          pd_interaction["Lost DDIs"]='&emsp;'+pd_interaction["Lost DDIs"]+'&emsp;'
          
          
          
          pd_interaction=pd_interaction.sort_values(by=['Percentage of lost domain-domain interactions'])
          #h=reverse('home')+"ID/"+trID+'/InteractionView/'

          
          
          pd_interaction["Percentage of lost domain-domain interactions"]=pd_interaction["Percentage of lost domain-domain interactions"].astype(int)
          pd_interaction["Percentage of lost domain-domain interactions"]=pd_interaction["Percentage of lost domain-domain interactions"].astype(str)+' % '
          
          
          

          
    else: pd_interaction=[]
    
    
    
    
    #get a list of all Isoforms of the selected transcript
    gene_ID=Ensemble_geneID
    transcripts=pr.gene_to_all_transcripts(gene_ID,organism)
    
    if len(transcripts)>=1:
  
          ID=[]
          name=[]
          pfams=[]
          
          for tr in transcripts :

              query = """
              SELECT * 
              FROM exons_to_domains_data_"""+organism+"""
              WHERE "Transcript stable ID"=:transcript_id 
              """
              tdata = pd.read_sql_query(sql=text(query), con=engine, params={'transcript_id': tr})
              
              # tdata=tdata.drop(columns=["Unnamed: 0"]).drop_duplicates()
              
              #df_filter = pr.data['Transcript stable ID'].isin([tr])
              #tdata=pr.data[df_filter]
              
              if len(tdata)!=0  :
                  ID.append(tr)
                  name.append(pr.tranID_convert(tr,organism)[0])
                  
                  p=tdata["Pfam ID"].unique()
                  p = p[~pd.isnull(p)]
                  p=sorted(p)
                  pfams.append(' ; '.join(p))
          
          
          if ID!=[]:        
            pd_isoforms=pd.DataFrame(list(zip(name, ID,pfams)), columns =['Transcript name', 'Transcript ID','Pfam domains'])
            pd_isoforms['length'] = pd_isoforms['Pfam domains'].str.len()
            pd_isoforms.sort_values('length', ascending=False, inplace=True)
            pd_isoforms=pd_isoforms.drop(columns=['length'])
            
            h=reverse('home')+"ID/"+organism+"/"
            pd_isoforms["Link"]='<a href="'+h+pd_isoforms["Transcript ID"]+'">'+" (Visualize) "+'</a>'
            
            
            
            
            pd_isoforms=pd_isoforms.to_html(**settings.TO_HTML_PARAMETERS)
            
    else: pd_isoforms=[]
    
    
    
    
    return nodes,edges,_,domains,unique_domains,exons,text1,domainshtml,Text_nodes,text_edges,tran_name,gene_name,Ensemble_geneID,entrezID,gene_description,exons,droped1,droped2,trID,p,df_missed,pd_interaction,pd_isoforms,co_partners









def Interacted_domain(p,g,entrezID,missing_domain):
        #Search for interacted domains
        edges=[]
        DDI_edges=[]
        lost_edges=[]
        DDI_edges2=[]
        lost_edges2=[]
        for e in g.edges():
            gene1=e[0].split('/')[0]
            gene2=e[1].split('/')[0]
            if p in [gene1,gene2]:
                edges.append(e)
                if edge_dashes(e,entrezID,missing_domain)[0]!='false':
                    lost_edges.append(e[0].split('/')[1]+'-'+e[1].split('/')[1])
                    lost_edges2.append(nt.link(e[0].split('/')[1]) + '-' + nt.link(e[1].split('/')[1]))


                elif len(e[0].split('/'))==2 and len(e[1].split('/'))==2: 
                    DDI_edges.append(e[0].split('/')[1]+'-'+e[1].split('/')[1])
                    DDI_edges2.append(nt.link(e[0].split('/')[1]) + '-' +nt.link( e[1].split('/')[1]))
        return edges,DDI_edges,lost_edges,DDI_edges2,lost_edges2

def table_interaction(tran_name,trID,entrezID,g,protein_with_DDI,missing_domain,organism):
    Interactions=[]
    IDs=[]
    DDIs=[]
    lost_DDIs=[]
    DDIs2=[]
    lost_DDIs2=[]
    perc=[]
    status=[]
    for protein in protein_with_DDI:
        #select first protein
        if protein !=entrezID:
            #Search for interacted domains
            edges,DDI_edges,lost_edges,DDI_edges2,lost_edges2=Interacted_domain(protein,g,entrezID,missing_domain)
            
            IDs.append(protein)
            p=len(lost_edges)/(len(DDI_edges)+len(lost_edges))
            p=float("{0:.2f}".format(p))
            perc.append(p*100)
            
            st='Affected'
            if p==0: st='Retained'
            if p==1: st='Missing'
            
            status.append(st)
            
            Interactions.append(pr.entrez_to_name(protein,organism))
            DDIs.append(' ; '.join(DDI_edges))
            lost_DDIs.append(' ; '.join(lost_edges))
            DDIs2.append(' ; '.join(DDI_edges2))
            lost_DDIs2.append(' ; '.join(lost_edges2))
    return  pd.DataFrame(list(zip(Interactions, IDs,DDIs2,lost_DDIs2,DDIs,lost_DDIs,perc,status)), columns =['Protein name', 'NCBI gene ID','Retained DDIs','Lost DDIs','retained DDIs','missing DDIs','Percentage of lost domain-domain interactions','Protein-protein interaction'])





























def vis_node_(g,entrezID,protein_with_DDI,tran_name,missing_domain,co_partners,organism):

    N=[]
    E=[]
    for node in g.nodes():
        #node of a protein:
        if node!=entrezID and len(node.split('/'))==1:
                try:
                    label=pr.entrez_to_name(node,organism)
                    N.append("{id: \""+node+"\", label:  \""+label+"\", group:  \""+group_node(node,entrezID)+"\", physics:"+physics(node,entrezID)+
                             ", source:  \""+source_node(node,entrezID,protein_with_DDI)+"\", value:  \""+value_node(node,entrezID)+"\"},")
                except KeyError:
                    print(' ')
 
                    
        else:
                    color=''
                    if node in missing_domain:  color='color: missing, '
                    label=node_label(node,entrezID,tran_name)
                    N.append("{id: \""+node+"\", label:  \""+label+"\", "+color+"group:  \""+group_node(node,entrezID)+"\", physics:"+physics(node,entrezID)+
                             ", source:  \""+source_node(node,entrezID,protein_with_DDI)+"\", value:  \""+value_node(node,entrezID)+"\"},")

    for e in g.edges():
         E.append("{from: \""+e[0]+"\", to: \""+e[1]+"\", dashes:  "+edge_dashes(e,entrezID,missing_domain)[0]+","+edge_option(e,entrezID,co_partners)+"},")      

    return N,E,len(g)-1























def edge_option(edge,entrezID,co_partners):
    e1=edge[0].split('/')
    e2=edge[1].split('/')
    #protein to protein
    if len(e1)==1 and len(e2)==1:
        edge_color=''
        #print(co_partners)
        if (e1[0] in co_partners and e1[0]!=entrezID) or  (e2[0] in co_partners and e2[0]!=entrezID): edge_color='color: Residue,'
        
        option='length: PR_LENGTH,'+edge_color+' width: WIDTH_SCALE * 4'
        
        
        
        
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
        source='DDI" , origin: "'+node.split('/')[0]
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
