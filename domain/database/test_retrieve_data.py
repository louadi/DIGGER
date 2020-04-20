import pandas as pd

from os import path
from sqlalchemy import create_engine, text

# --- Create database connection aka 'SQLAlchemie engine'
# sqlite://<no_hostname>/<path>
# where <path> is relative:
engine = create_engine('sqlite:///datasets.db')
metadata_engine = engine

# --- Get tables from database
query = """
        SELECT DISTINCT "Transcript stable ID_x", "u_ac_1", "Transcript stable ID_y", "u_ac_2"
        FROM ppi_data 
        WHERE "Exon stable ID_x"=:exon_id
        """
p1_sql = pd.read_sql_query(sql=text(query), con=engine, params={'exon_id': 'ENSE00003756482'})
p1_sql

assert(False)
