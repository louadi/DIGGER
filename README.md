# DIGGER
<p align="center">
  <img src="https://exbio.wzw.tum.de/digger/static/image/DIGGER.png" height="130">
</p>



Protein-protein interaction (PPI) networks are a key resource for systems biology. However, they do not consider the influence of alternative splicing, even though experimental evidence suggests that interaction partners are different for isoforms of the same protein. Domain Interaction Graph Guided ExploreR (DIGGER) integrates protein-protein interactions and domain-domain interactions into a joint graph and maps interacting residues to exons. DIGGER allows the users to query exons or isoforms individually or as a set to visually explore their interactions and it is available at: https://exbio.wzw.tum.de/digger




## Deploying DIGGER
To install DIGGER (running it for the first time) follow these steps:
```shell script
# First, follow this link, if you want to install docker-compose: 
# https://docs.docker.com/compose/install/

#  Clone this repository and change into the created directory
git clone https://github.com/louadi/DIGGER.git && cd DIGGER

# Download a copy of all the data files into domain/data
wget https://zenodo.org/record/3973368/files/data.zip    # extract the files

# Create a copy of the .env.sample file and edit the .env file
cp .env.sample .env   # now edit the .env file 

# Deploy and build the containers
docker-compose up -d --build

# Apply migrations to the database (make sure the containers are up an running, else you will get an error)
docker-compose exec web python manage.py migrate --noinput 

# Import all the datasets into the database
docker-compose exec web python manage.py import_datasets

# Collect all the static files
docker-compose exec web python manage.py collectstatic --no-input

# Enjoy your instance of DIGGER

```


## Cite

If you use DIGGER, please cite:


Zakaria Louadi, Kevin Yuan, Alexander Gress, Olga Tsoy, Olga Kalinina, Jan Baumbach, Tim Kacprowski*, Markus List*. DIGGER: exploring the functional role of alternative splicing in protein interactions, Nucleic Acids Research, https://doi.org/10.1093/nar/gkaa768  (*joint last author)



## Contact us
Zakaria Louadi: louadi@wzw.tum.de
