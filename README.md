# DIGGER

[Some project description]


## Deploying DIGGER
To install Digger (running it for the first time) follow these steps:
```shell script
# First clone this repository and change into the created directory
git clone https://github.com/louadi/DIGGER.git && cd DIGGER

# Download a copy of all the data files into domain/data
wget .... [maybe publish a zip archive with our files?]

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
