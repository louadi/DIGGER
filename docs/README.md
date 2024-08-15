<p align="center">
  <img src="https://exbio.wzw.tum.de/digger/static/image/DIGGER.png" height="130">
</p>

Domain Interaction Graph Guided ExploreR (DIGGER) helps to investigate the impact of alternative splicing on protein-protein interactions (PPI) on the level of an exon, event, transcript, or the whole network. DIGGER is ideally suited to investigate the difference in interactions
between the isoforms, analyze the effect of isoform switch, or to explore how alternative splicing events such as exon skipping lead to altered
interactions of protein isoforms.   

# The concept of DIGGER

Experimental evidence suggests that the majority of protein isoforms have different interaction partners. However, conducting experiments for all possible isoforms will be time- and resource-consuming. Domain Interaction Graph Guided ExploreR (DIGGER) integrates protein-protein and domain-domain interactions (DDI) into a joint graph and maps interacting residues to exons (Figure, **A**, upper part). The concept is the following: if an alternative splicing event leads to splicing out an interacting domain or residues, the interaction between proteins might be lost or impaired (Figure, **A**, lower part). 

<p align="center">
<img alt="DIGGER workflow" src="https://exbio.wzw.tum.de/digger/static/image/figure%201.png" width="800"/>
</p>

DIGGER uses the following sources for PPI, DDIs, and interacting residues:

Homo sapiens:
[TODO]
Mus musculus:
[TODO]

DIGGER provides PPI visualizations in a structural context using a dynamic graph visualization that can be toggled between a protein
isoform and a domain-centric view (Figure, **B**). Furthermore, the tool maps the protein features encoded by a selected exon,
to judge the functional role of individual exons in the PPI (Figure, **C**).  DIGGER’s joint PPI and domain-domain interaction network can also be used for subnetwork extraction, providing a basis for network analysis (**D**).

Finally, DIGGER provides a web interface for **NEASE** (Network Enrichment method for Alternative Splicing Events) - a tool for the functional enrichment of alternative splicing exons/events.

# Step-by-step functional analysis of alternative splicing

# Isoform-level analysis

<img align="left" src="https://github.com/OlgaVT/DIGGER/blob/patch-1/docs/Figure%202.png" width="64"/>

The isoform-level analysis investigates the difference in interactions between isoforms or the consequences of isoform switch. In the documentation, we will use an example of BAG1 and NCK2 genes from the DIGGER publication. 

As input, isoform-level analysis accepts the name of a gene, an Ensebml identifier of a gene, and an Ensembl identifier of a transcript. A user can also put multiple identifiers separated by a comma. Then a user should select an organism. Currently, DIGGER contains the data from Homo sapiens and Mus musculus.

*Example: Input genes are BAG1 and NCK2; the chosen organism is Homo sapiens.*
<p align="center">
</p>

If there is information about interactions for the isoforms of interest, the resulting window will show a summary of information about them in a network ('Network') and tabular views ('Queries'). A user can switch between them in the menu in the upper left corner.

A network view provides information about known interactions of query isoforms. Predicted DDIs toggle bars allow the addition of predicted DDIs of different levels of confidence (**A**). Interactions lost or impaired due to splicing are marked as dashed lines. A user could also adjust the view of the network by enabling/disabling physics and adjusting the height of a canvas (**B**).
<p align="center">
</p>

The tabular view presents the information about query isoforms or all annotated protein-coding isoforms if a query consists of a gene name/identifier. The Pfam domains (**A**) column contains information about known PFAM domains of a resulting protein. "Present / Missing interacting domains in the isoform" column (**B**) demonstrates the percentage of interactions lost or impaired for this particular isoform due to alternative splicing. To further investigate the impact of a splicing event, a user can click 'Visualize' (**C**).

*Example: a transcript BAG-207 possibly lost a PPI due to alternative splicing that leads to removing the part of the protein that encodes a domain PF02179. A transcript NCK-202 also lost a PPI due to splicing of domains PF00017 and PF14604*
<p align="center">
</p>

A "Visualize" window, provides five different visualization panels: Overview, Protein View, Interaction View, Domain View, and Other Isoforms.

The Overview window header shows the general information about the transcript.
<p align="center">
</p>

Below is a scheme of transcript exons and mapped Pfam domains, and the table of all exons and Pfam domains mapped on exons. 

*Example: Comparing BAG1-210 and BAG-207 isoforms, one can detect that these isoforms differ in exons that encode a domain PF02179 as well as exons that contain interacting residues.*
<p align="center">
</p>

*Example: Comparing NCK-201 and NCK-202 isoforms, one can detect that these isoforms differ in exons that encode domains PF14604, the second PF00018, PF00017 as well as exons that contain interacting residues.*
<p align="center">
</p>

The Protein view shows all known direct interactors of an isoform. View mode can be changed from PPI to PPI-DDI view. 

*Example: PPI-DDI view for BAG1-207 depicts possibly missing interactions with BAG3, HSPA1A, HSPA4, and HSPA8.*
<p align="center">
</p>

Interaction View allows investigation into details of the PPI between a query isoform and all known and predicted interactors. The interactors can be chosen from the Select Interaction Partner Menu (just tick the box with the interaction of interest there).

*Example: Interaction View depicts the missing interaction between BAG1-207 and HSPA8*
<p align="center">
</p>

Domain View shows all interactors known for a particular domain.

*Example: Splicing out a domain PF02179 in the isoform BAG1-207 leads to missing PPIs with BAG3, HSPA1A, HSPA4, and HSPA8 comparing to BAG1-201.*
<p align="center">
</p>

Other isoforms View again depicts the table of all other known protein-coding isoforms of a gene.
<p align="center">
</p>

# Exon-level analysis

The Exon-level analysis investigates the consequences of exon skipping event to interactions between isoforms. In the documentation, we will use an example of the INSR gene from the DIGGER publication.
<p align="center">
</p>

In addition to input as in isoform-level analysis, exon-level analysis accepts an Ensebml identifier of an exon, an Ensembl identifier of a gene with the coordinates of an exon. A user can also put multiple identifiers separated by a comma. Then a user should select an organism.

*Example: an input exon is exon ENSE00003569638 from BAG1; the organism is Homo sapiens*
<p align="center">
</p>

Results of the exon analysis consists of DomainView, ProteinView, InteractionView, and the table of Proteins that use this exon.

DomainView presents the general information about an exon (**A**), the information about the corresponding Pfam domain with a link to Pfam and 3did databases (**B**) and a network with all known interactions (**C**). Using toggle bars, a user can add predicted DDIs of different level of confidence (**D**).

*Example: ENSE00003569638 encode a domain PF02719. This domain is known to interact with BAG through a DDI with PF02719 and HSPA1A, HSPA4, and HSPA8 through a DDI with PF00012.*
<p align="center">
</p>

ProteinView is similar to Isoform-level analysis.

*Example: Skipping of ENSE00003569638 might lead to loss of interactions with BAG, HSPA1A, HSPA4, and HSPA8.*
<p align="center">
</p>

InteractionView is also similar to Isoform-level analysis.

Example: Skipping of ENSE00003569638 might lead to loss of interactions with HSPA8.
<p align="center">
</p>

Finally, a table "Proteins that use this exon" give use the list of such protein-coding isoforms:

<p align="center">
</p>

# Network-level analysis
Network analysis investigates the systematic impact of alternative splicing events on the PPIs.
<p align="center">
</p>

As input, Network-level analysis accepts the list of genes, transcripts, proteins, or a transcript count file and constructs a subnetwork from this list.
<p align="center">
</p>

If interactions between isoforms are known, the isoforms are connected to the subnetworks and the interactions are marked as lost or present as well as based on the interaction evidence. For further exploration, the network can be downloaded in SIF or GraphML formats.
<p align="center">
</p>


# Running NEASE

**NEASE** (Network Enrichment method for Alternative Splicing Events) is a splicing-based functional network enrichment tool. Currently, NEASE supports the following formats:
- Standard (List of genes with exon coordinates)
- MAJIQ *deltapsi.tsv
- Whippet *.diff
- rMATS SE.MATS.JC.txt

If you have a lot of files to run with NEASE, consider using the Python package: https://github.com/louadi/NEASE.

To start the analysis, upload your file (**A**), define an organism (currently, Homo sapiens or Mus musculus, **B**), the input file type (**C**), and adjust the parameters.

For all input files except MAJIQ:
- Analysis name
- Usage of predicted DDIs (**D**)
- P-value cutoff (**E**)
- Minimum delta PSI (Min delta, **F**)

For MAJIQ:
- Minimum delta PSI will start by default from 0.2 (Min delta, **F**)
- Confidence interval (0.95 corresponds to P-value 0.05)
<p align="center">
</p>

The results of the analysis will be stored for **7 days** and visible from the NEASE starting page (Previous Analysis). For large files, the running time might take up to several minutes.

The resulting 'Overview' page contains: 
- the general information about the number of genes with AS that affect any known protein feature;
- the table of affected protein domains;
- the table with affected interactions;
- the table with affected linear features;
- the table with affected interacting residues.

All tables are downloadable in a csv format.

NEASE contains three further types of analysis: enrichment analysis, pathway analysis, and visualization.

NEASE edge enrichment supports the following ontologies:

Homo sapiens: TODO
Mus musculus: KEGG, Reactome, MouseCyc.
<p align="center">
</p>

Simply type the name of the ontology of interest and click the "Run enrichment" button. The resulting table will show the list of pathways, the list of genes that lost/impaired interactions with these pathways due to alternative splicing, and the significance assessment with p-value and adjusted p-value. NEASE score re-weights the pathways according to the hub bias. For more details, check the publication:  https://doi.org/10.1186/s13059-021-02538-1.

To investigate the enriched pathways in detail, use Pathway analysis. Type in the pathway ID from NEASE edge enrichment analysis, e.g., "path:mmu04010". The resulting table contains the information about each gene with list/impaired interactions with this pathway due to alternative splicing: whether it belongs to the pathway, how many interactions were lost/impaired, and the p-value.

Finally, "Visualize pathways" option provides the visualization of teh resulting network with lost/impaired interactions. Be aware, the visualization might take up to several minutes. 

# Hosting your own instance

If you want, you can also host your own instance of DIGGER. To do so, follow the instructions below.

## Deploying DIGGER
To install DIGGER 2.0 (running it for the first time) follow these steps:
```shell script
# First, follow this link, if you want to install docker-compose: 
# https://docs.docker.com/compose/install/

#  Clone this repository and change into the created directory
git clone https://github.com/louadi/DIGGER.git && cd DIGGER

# Download a copy of all the data files into domain/data/Homo sapiens[human]/
# If you have more organisms, add the files into respective foldes, e.g. domain/data/Mus musculus[mouse]/
wget https://zenodo.org/record/3973368/files/data.zip    # extract the files


# Create a copy of the .env.sample file and edit the .env file
cp .env.sample .env   # now edit the .env file 

# Deploy and build the containers
docker-compose up -d --build

# Apply migrations (make sure the containers are up an running, else you will get an error)
docker-compose exec web python manage.py migrate --noinput 

# Import all the datasets into the database
docker-compose exec web python manage.py import_datasets

# Collect all the static files
docker-compose exec web python manage.py collectstatic --no-input

# Enjoy your instance of DIGGER

```


## Extending DIGGER by DDIs
If you downloaded the data provided by us, the extended DDIs are already included and you will **not** have to do 
this step. Only extend DIGGER like this if you want to add new organisms and have sufficient data to do so.  
To extend DIGGER using generated Domain-Domain Interactions (DDIs) follow these steps:
````bash
# First, create a conda environment to ensure you have all dependencies installed
conda env create -f DIGGER_env.yml

# Activate the environment
conda activate DIGGER

# Make sure you have the necessary files in the sourcedata folder
# more info about what these files should look like can be found in the sourcedata README.md

# copy the example database_sources to have a backup and list of all options. 
# Edit the database_sources.yml file to your needs
cp preprocess/sourcedata/example.database_sources.yml preprocess/sourcedata/database_sources.yml

# Run the prediction script
cd preprocess
python main_ddi_extend.py

# Once this has run, depending on your settings, restarting the DIGGER container will show the new data
cd ..
docker-compose up -d --force-recreate
````


# Cite

If you use DIGGER, please cite:

Zakaria Louadi, Kevin Yuan, Alexander Gress, Olga Tsoy, Olga Kalinina, Jan Baumbach, Tim Kacprowski*, Markus List*. DIGGER: exploring the functional role of alternative splicing in protein interactions, Nucleic Acids Research, https://doi.org/10.1093/nar/gkaa768  (*joint last author)


# Contact us
Elias Albrecht: elias.albrecht@in.tum.de  
Olga Tsoy: olga.tsoy@uni-hamburg.de
