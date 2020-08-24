# Python imports
from os.path import join
from sqlalchemy import create_engine

# project imports
from .common import *
from .i18n import *

# Dummy variable to import so PyCharm does not remove it while optimizing imports
__dummy = LANGUAGE_CODE

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True    # ToDo change for production

# Restrict access to only the hostname that this service is deployed to
ALLOWED_HOSTS = ['*']   # ToDo change for production

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
DATABASE_ENGINE = create_engine(f'postgres://{db_user}:{db_password}@{db_host}:{db_port}/{db_name}')

# ##### APPLICATION CONFIGURATION #########################

INSTALLED_APPS = DEFAULT_APPS
