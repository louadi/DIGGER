# Python imports
from os.path import join
from sqlalchemy import create_engine

# project imports
from .common import *
from .i18n import *

# Dummy variable to import so PyCharm does not remove it while optimizing imports
__dummy = LANGUAGE_CODE

# ##### DEBUG CONFIGURATION ###############################
DEBUG = True

# allow all hosts during development
ALLOWED_HOSTS = ['*']

# ##### DATABASE CONFIGURATION ############################
sqlite_file_path = join(PROJECT_ROOT, 'run', 'dev.sqlite3')
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': sqlite_file_path,
    }
}

# --- Create database connection aka 'SQLAlchemie engine'
# Database for analysis (SQLalchemie engine)
# sqlite://<no_hostname>/<path>
# where <path> is relative:
DATABASE_ENGINE = create_engine(f'sqlite:///{sqlite_file_path}')

# ##### APPLICATION CONFIGURATION #########################

INSTALLED_APPS = DEFAULT_APPS
