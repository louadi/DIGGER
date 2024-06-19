<p align="center">
  <img src="https://exbio.wzw.tum.de/digger/static/image/DIGGER.png" height="130">
</p>

Protein-protein interaction (PPI) networks are a key resource for systems biology. However, they do not consider the influence of alternative splicing, even though experimental evidence suggests that interaction partners are different for isoforms of the same protein. Domain Interaction Graph Guided ExploreR (DIGGER) integrates protein-protein interactions and domain-domain interactions into a joint graph and maps interacting residues to exons. DIGGER allows the users to query exons or isoforms individually or as a set to visually explore their interactions and it is available at: https://exbio.wzw.tum.de/digger

# Overview of the workflow

<img alt="DIGGER workflow" src="https://exbio.wzw.tum.de/digger/static/image/figure%201.png" width="800"/>

PPI networks do not consider the influence of alternative splicing, even
though experimental evidence suggests that the majority of protein isoforms
have different interaction partners.

Domain Interaction Graph Guided ExploreR (DIGGER) integrates
protein-protein interaction and domain-domain interactions into a joint
graph and maps interacting residues to exons.

DIGGER provides a great way to visualize PPI in structure context using a
dynamic graph visualization that can be toggled between a protein
isoform and a domain-centric view.  
By comparing the structure of different isoforms, DIGGER identifies the
surfaces specific to an isoform and missing in others. In this way,
isoform-specific interactions can be precisely identified.  
Furthermore, the tool maps the protein features encoded by a selected exon,
to judge the functional role of individual exons in the PPI.  
DIGGER is ideally suited to investigate the difference of interactions
between the isoforms, analyze the effect of isoform-switch, To or explore
how alternative splicing events such as exon skipping lead to altered
interactions of protein isoforms.   
DIGGER’s joint PPI and domain-domain interaction network can also be
used for subnetwork extraction, providing a basis for network analysis. From
a list of proteins or transcripts, the tool re-score the PPI edges based on
the structure evidence of the input proteins.  

# Workflow for alternative splicing analysis

# Isoform-level analysis

# Exon-level analysis

# Network-level analysis

# Running NEASE

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
