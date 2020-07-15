# DIGGER

Domain Interaction Graph Guided ExploreR (DIGGER) integrates protein-protein interaction and domain-domain interactions into a joint graph and maps interacting residues to exons. DIGGER allows the users to query exons or isoforms individually or as a set to visually explore their interactions. 

## Deploying DIGGER
To install DIGGER (running it for the first time) follow these steps:
```shell script
# First, follow this link, if you want to install docker-compose: https://docs.docker.com/compose/install/

#  Clone this repository and change into the created directory
git clone https://github.com/louadi/DIGGER.git && cd DIGGER

# Download a copy of all the data files into domain/data
wget https://zenodo.org/record/3885677/files/data.rar      # extract the files

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
