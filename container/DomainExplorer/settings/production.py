# Python imports
import os
from os.path import join
from sqlalchemy import create_engine

# project imports
from .common import *
from .i18n import *

# Dummy variable to import so PyCharm does not remove it while optimizing imports
__dummy = LANGUAGE_CODE

# SECURITY WARNING: don't run with debug turned on in production!
debug_env = os.environ.get("DEBUG", "")
DEBUG = True if debug_env == 'True' else False

# Restrict access to only the hostname(s) this service is deployed to.
# Set ALLOWED_HOSTS in .env as a comma-separated list; falls back to the public host.
ALLOWED_HOSTS = [h.strip() for h in os.environ.get('ALLOWED_HOSTS', '').split(',') if h.strip()] \
    or ['exbio.wzw.tum.de']

# Database
# https://docs.djangoproject.com/en/2.2/ref/settings/#databases

# ##### DATABASE CONFIGURATION ############################
db_name = os.environ.get('POSTGRES_DB', 'postgres')
db_user = os.environ.get('POSTGRES_USER', 'postgres')
db_password = os.environ.get('POSTGRES_PASSWORD', 'postgres')
db_host = 'db'
db_port = '5432'
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql',
        'NAME': db_name,
        'USER': db_user,
        'PASSWORD': db_password,
        'HOST': db_host,
        'PORT': db_port,
    },
}

# --- Create database connection aka 'SQLAlchemie engine'
# Database for analysis (SQLalchemie engine)
# sqlite://<no_hostname>/<path>
# where <path> is relative:
DATABASE_ENGINE = create_engine(f'postgresql://{db_user}:{db_password}@{db_host}:{db_port}/{db_name}')

# ##### APPLICATION CONFIGURATION #########################

INSTALLED_APPS = DEFAULT_APPS
