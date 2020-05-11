# project imports
import os

from .common import *

# Enable debugging in development
DEBUG = True

# Database
# https://docs.djangoproject.com/en/2.2/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASE_DIR, '../../db.sqlite3'),
    }
}

# --- Create database connection aka 'SQLAlchemie engine'
# Database for analysis (SQLalchemie engine)
# sqlite://<no_hostname>/<path>
# where <path> is relative:
DATABASE_ENGINE = create_engine('sqlite:///datasets.db')