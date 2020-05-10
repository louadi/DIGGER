from django.shortcuts import render
from django.http import HttpResponse
from .load.domaingraph import load_obj
from django_pandas.managers import DataFrameManager
from django.shortcuts import redirect

from .Process import exonstodomain as exd 
from .Process import exon as ex 
from .Process import process_data as pr
from .Process import proteininfo as  info
from .Process import transcript as  tr
from .Process import InteractionView as  iv
from .Process import gene as  g
from .Process import network_analysis as nt 


import pickle
import random 

'''
def get_input(query=None):
    queryset=[]
    queries=query.split('')
'''


# or reference a template from templates folder
def home(request):

    
    if "search" in request.GET :     # If the form is submitted
    
        search_query = request.GET['search']
        search_query=search_query.replace(" ", "")
        search_query=search_query.split("+")[0]
        search_query=search_query.split("%")[0]
        search_query=search_query.split(".")[0]
        #Input search is a protein:
        if len(search_query)==15 and search_query[:4]=='ENST' or search_query[:4]=='ENSP':
          return transcript(request,search_query)
          
          
          
          
        #Input search is an exon: 
        elif search_query[:4]=='ENSG':
          return gene(request,search_query)
     
     
    if "search 2" in request.GET :     # If the form is submitted
    
        search_query = request.GET['search 2']
        search_query=search_query.replace(" ", "")
        search_query=search_query.split("+")[0]
        search_query=search_query.split("%")[0]
        search_query=search_query.split(".")[0]

        if search_query[:4]=='ENSE':
          return exon(request,search_query)
        
    return render(request,'domain/home.html' ) 



#Display transcripts of a gene
def gene(request,gene_ID):
  
   transcript_table,gene_name=g.input_gene(gene_ID)

   if transcript_table==[]:
       return HttpResponse(' wrong entry or protein without any known Pfam domains')
   
   context = {
      'tb':transcript_table,
      'name':gene_name,
       }  
   return render(request,'domain/gene.html',context)  





#Input is an exon:

def exon(request,exon_ID):
    
   v=ex.input_exon(exon_ID)
   
   if v==True:
       return HttpResponse(' wrong entry or exon in a gene without any known Pfam domains')
   else:
       _,domains, gene_name,Ensemble_geneID,entrezID,tb_transc,table_domains,number=v
       
       
    #only if the exon code for domains with known interactions
   
    
   if number >0 :
           
           nodes,edges,pd_interaction=ex.vis_exon(domains,entrezID,gene_name,exon_ID)
   else:
        nodes,edges,pd_interaction=[],[],[]
        
   #PPI res interfaces on the exon:      
   #table: HTML table with all PPIs that have res interface mapped to the exon
   #number_of_PPI: number of interactions
   
   table,number_of_PPI=ex.PPI_inter(exon_ID,gene_name)
   
   
   
   if number >0 :
                # added to combine evidence of DDI and Residue in one final table
                if number_of_PPI>0:
                      
                      ppi_from_res=table['Partner Protein'].unique()
                      f=pd_interaction['Partner Protein'].isin(ppi_from_res)
                      pd_interaction.loc[f, 'Residue evidence'] = '<center>&#128504;</center>'
                 
                
                pd_interaction["Partner Protein"]='<center>'+pd_interaction["Partner Protein"]+'</center>'
                
                pd_interaction=pd_interaction[["Affected Protein",'Partner Protein','NCBI gene ID','Retained DDIs','Lost DDIs','Percentage of lost domain-domain interactions','Residue evidence',"Protein-protein interaction"]]
                
                
                pd_interaction=pd_interaction.rename(columns={
                "Partner Protein": "<center>Partner Protein</center>", 
                "Affected Protein": "<center>Affected Protein</center>",
                "NCBI gene ID": "<center>NCBI gene ID</center>", 
                "Percentage of lost domain-domain interactions": "<center> % of lost DDIs</center>",
                "Retained DDIs": "<center>&emsp;Retained Domain-Domain interactions</center>", 
                "Lost DDIs": "<center>Lost Domain-Domain interactions</center>",
                "Protein-protein interaction": "<center>Protein-protein interaction</center>",
                'Residue evidence':'<center>Residue evidence</center>'
                })


                pd_interaction=pd_interaction.to_html(escape=False, index=False)
             
   table=table.to_html(escape=False, index=False)  
   context = {
   
      'tb1':tb_transc,
      'tb2':table_domains,
      'tb3':pd_interaction,
      'tb4':table,
      'name':gene_name,
         'exon_ID':exon_ID,
    'entrezID':entrezID,
    'gID':Ensemble_geneID,
    'dis':number>0,
    "dis2": number==-1,
    'dis3':number_of_PPI!=0,
    'long_table': number_of_PPI>25,
    'pv_nodes': nodes,
    'pv_edges': edges,
    
       }  
   return render(request,'domain/exon.html',context)  







#Display a single node
def display(request,Pfam_id):
  
   n,e,gene,domain=exd.vis_node_(Pfam_id)
   context = {
     'nodes':n,
     'edges':e,
     'gene':gene,
     'domain':domain
       }  
   return render(request,'domain/display.html',context)  



#Dsiplay information of a transcript or a protein
def transcript(request,P_id):


    print(P_id)
    
    nodes,edges,_,domains,unique,exons,text1,domainshtml,Text_nodes,text_edges,tran_name,gene_name,Ensemble_geneID,entrezID,gene_description,exons,droped1,droped2,trID,p,missed,interaction,isoforms,=tr.Protein_view(P_id)
    
    '''
    dis1=True
    if missed==[]: dis1=False
    
    
    dis2=True
    if interaction==[]: dis2=False
    
    
    dis3=True
    if isoforms==[]: dis3=False
    '''
        
    context={
    'dt':droped1,
    'text1':text1,
    'tran_name':tran_name,
    'gene_description':gene_description,
    'trID':trID,
    'gID':Ensemble_geneID,
    'path':p,
    'pv_nodes': nodes,
    'pv_edges': edges,
    'entrezID' : entrezID,
    'dt2': droped2,
    'dt3' :missed,
    'dt4' :interaction,
     "dt5": isoforms,
     
    'dis1': missed!=[],
    'dis2': interaction!=[],
  
    'dis3': isoforms!=[],
    
    }
 
    return render(request,'domain/transcr.html',context) 







#InteractionView
def InteractionView(request,P_id,P2_id):

    
    
    
    nodes,edges,interaction,tr,p_name=iv.int_view(P_id,P2_id)
    
  
    
    context={

    'pv_nodes': nodes,
    'pv_edges': edges,

    'dt1' :interaction,
    'pr':p_name,
    'tr':tr,
    
    }
 
    return render(request,'domain/InteractionView.html',context) 






def isoform_level(request):
  if "search" in request.GET :     # If the form is submitted
    
        search_query = request.GET['search']
        search_query=search_query.replace(" ", "")
        search_query=search_query.split("+")[0]
        search_query=search_query.split("%")[0]
        search_query=search_query.split(".")[0]
        search_query=search_query[:15]
        #Input search is a protein:
        if  search_query[:4]=='ENST' or search_query[:4]=='ENSP':
          return redirect(transcript, P_id = search_query)
          #return transcript(request,search_query)
          
          
          
          
        #Input search is an exon: 
        elif len(search_query)==15 and search_query[:4]=='ENSG':
          return redirect(gene, gene_ID = search_query)
          #return gene(request,search_query)
     
  return render(request,'domain/isoform_level.html',) 






def Exon_level(request):
    if "search 2" in request.GET :     # If the form is submitted
      #Input and Exon ID
      print('-----------------------------------------------------------')
      search_query = request.GET['search 2']
      search_query=search_query.replace(" ", "")
      search_query=search_query.split("+")[0]
      search_query=search_query.split("%")[0]
      search_query=search_query.split(".")[0]
      search_query=search_query[:15]


      if  search_query[:4]=='ENSE':
          return redirect(exon, exon_ID = search_query)
          #return exon(request,search_query)   
    




    if "search" in request.GET :     # If the form is submitted
        #Input coordinate of the exon
        #Check if coordinate are correct
        #  Example   ' ENSG00000266028  206437964 206437042 '
    
        print('-----------------------------------------------------------')
        search_query = request.GET['search']
        
        
        search_query=search_query.split(" ")
        search_query =[x for x in search_query if x!='']
        #search_query[0]=search_query[0].split(".")[0]
        if len(search_query)==3 and  len(search_query[0])==15 and search_query[0][:4]=='ENSG' and search_query[1].isdigit() and search_query[2].isdigit():

            gene_ID=search_query[0]
            s1=int(search_query[1])
            e1=int(search_query[2])
         
            exonID=pr.coordinate_to_exonID(gene_ID,s1,e1)
            
            if exonID!=[]:
                return redirect(exon, exon_ID = exonID)
                #return exon(request,exonID)   
            else: return HttpResponse("<h1>No match</h1>")
    return render(request,'domain/Exon_level.html',) 
    
    
    

#PPI network analysis
def network(request):
    
    if "input" in request.GET :
          input_query = request.GET['input']
          print(input_query)
          #print(input_query)
          input_query=input_query.split("\r\n")
          #print(input_query)
          input_query[0]=input_query[0].replace(" ", "")
          
          #max input IDs
          if len(input_query)<2000 and len(input_query)>1:
                if input_query[0][0:4]=='ENSG' or input_query[0][0:4]=='ENST' or input_query[0][0:4]=='ENSP':
                      job_num=str(random.randrange(500))
                      with open('domain/static/jobs/'+job_num+".txt", "wb") as fp:   #Pickling
                             pickle.dump(input_query, fp)
                      return redirect(Multi_proteins,job=job_num)
                
    return render(request,'domain/network.html')  



def Multi_proteins(request, job='0'):

    with open('domain/static/jobs/'+job+".txt", "rb") as fp:   # Unpickling
            inputs = pickle.load(fp)
            
    
    if inputs[0][0:4]=='ENSG' :
       info=nt.analysis_input_genes(inputs)
       
    elif inputs[0][0:4]=='ENSG' or inputs[0][0:4]=='ENST' or inputs[0][0:4]=='ENSP':
          info=nt.analysis_input_isoforms(inputs)
    else: 
           return HttpResponse("<h1>wrong entry</h1>")
           
    if info==False:
        return HttpResponse("<h1>Too many inputs (max=2000 genes)</h1>")
    
    else:    
      genes, missing,num_isoforms=info
      
      Net=nt.Construct_network(genes, missing,job)
      
      if Net==0:  
          return HttpResponse("<h1>There is no known interaction between these proteins</h1>")
      
      else: nodes,edges,tab,tb_html=Net

    
    

    context={

    'pv_nodes': nodes,
    'pv_edges': edges,
    "tab":tb_html,
    'ID':job,
    "genes_number": len(missing),
    "isoforms_num": num_isoforms,
    'interacted_nodes':len(nodes),
    
    }
    
    return render(request,'domain/vis_network.html',context) 

def example2(request):

    if "input" in request.GET :
          input_query = request.GET['input']
          print(input_query)
          #print(input_query)
          input_query=input_query.split("\r\n")
          #print(input_query)
          input_query[0]=input_query[0].replace(" ", "")
          
          #max input IDs
          if len(input_query)<2000:
                if input_query[0][0:4]=='ENSG' or input_query[0][0:4]=='ENST' or input_query[0][0:4]=='ENSP':
                      job_num=str(random.randrange(500))
                      with open('domain/static/jobs/'+job_num+".txt", "wb") as fp:   #Pickling
                             pickle.dump(input_query, fp)
                      return redirect(Multi_proteins,job=job_num)
    return render(request,'domain/Network_example2.html')   
    
    
    
        
def example1(request):

    if "input" in request.GET :
          input_query = request.GET['input']
          print(input_query)
          #print(input_query)
          input_query=input_query.split("\r\n")
          #print(input_query)
          input_query[0]=input_query[0].replace(" ", "")
          
          #max input IDs
          if len(input_query)<2000:
                if input_query[0][0:4]=='ENSG' or input_query[0][0:4]=='ENST' or input_query[0][0:4]=='ENSP':
                      job_num=str(random.randrange(500))
                      with open('domain/static/jobs/'+job_num+".txt", "wb") as fp:   #Pickling
                             pickle.dump(input_query, fp)
                      return redirect(Multi_proteins,job=job_num)
    return render(request,'domain/Network_example1.html')   
    
def example3(request):

    if "input" in request.GET :
          input_query = request.GET['input']
          print(input_query)
          #print(input_query)
          input_query=input_query.split("\r\n")
          #print(input_query)
          input_query[0]=input_query[0].replace(" ", "")
          
          #max input IDs
          if len(input_query)<2000:
                if input_query[0][0:4]=='ENSG' or input_query[0][0:4]=='ENST' or input_query[0][0:4]=='ENSP':
                      job_num=str(random.randrange(500))
                      with open('domain/static/jobs/'+job_num+".txt", "wb") as fp:   #Pickling
                             pickle.dump(input_query, fp)
                      return redirect(Multi_proteins,job=job_num)
    return render(request,'domain/Network_example4.html')   
    
    
def about(request):
 
    return render(request,'domain/about.html',) 
 

 
def graph(request):
    jdata = {
     'data':jsonData,
     'data2':jsonData2
       }
    return render(request,'domain/index.html',jdata) 







def homeold(request):
  return HttpResponse("<h1>Blog home</h1>") 
  
  
  
  
  
  
  
  
  
  
  
  



  
  
  
  
  
  
  
  
  
