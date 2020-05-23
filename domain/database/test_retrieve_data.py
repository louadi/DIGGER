import pandas as pd
import numpy as np

from os import path
from sqlalchemy import create_engine, text

# --- Create database connection aka 'SQLAlchemie engine'
# sqlite://<no_hostname>/<path>
# where <path> is relative:
engine = create_engine('sqlite:///../../run/dev.sqlite3')


# --- Get ppi_data tables from database
def get_ppi_data():
    query = """
            SELECT DISTINCT "Transcript stable ID_x", "u_ac_1", "Transcript stable ID_y", "u_ac_2"
            FROM ppi_data 
            WHERE "Exon stable ID_x"=:exon_id
            """
    p1_sql = pd.read_sql_query(sql=text(query), con=engine, params={'exon_id': 'ENSE00003756482'})


# --- Get exons_to_domains_data tables from database
def getexons_to_domains_data():
    # The input transcript ID
    # Outputs
    # Exons: all exons in the transcript
    # D: mapping coding exons to their respective protein domains
    # uniaue domains in the protein
    transcript_ID = 'ENST00000368192'
    # SQL
    query = """
            SELECT * 
            FROM exons_to_domains_data 
            WHERE "Transcript stable ID"=:transcript_id 
            ORDER BY "Exon rank in transcript"
            """
    tdata = pd.read_sql_query(sql=text(query), con=engine, params={'transcript_id': transcript_ID})

    # OLD
    data_old = pd.read_csv('../data/final.csv')
    df_filter = data_old['Transcript stable ID'].isin([transcript_ID])
    tdata_old = data_old[df_filter].sort_values(by=['Exon rank in transcript'])
    exons = tdata.drop(columns=["Pfam ID", "Pfam start", "Pfam end"]).drop_duplicates()

    # COMPARE
    assert (np.array_equal(tdata.values, tdata_old.values))

    # Keep
    D = tdata[tdata["Pfam ID"].notna()].drop_duplicates()
    p = D["Pfam ID"].unique()

    return exons, D, p


# --- Get name from transcript id (mapper)
def transcript_id_to_name():
    transcript_ID = 'ENST00000368192'
    # SQL
    query = """
            SELECT "Transcript name" 
            FROM gene_info 
            WHERE "Transcript stable ID"=:transcript_id
            LIMIT 1
            """
    tdata = \
    pd.read_sql_query(sql=text(query), con=engine, params={'transcript_id': transcript_ID}).iloc[0, 0].split('-')[0]
    print(tdata)


# --- Get all transcripts that contains a specific domain in the genome
def transcript_with_domain():
    pfam_id = 'PF00178'
    # SQL
    query = """
            SELECT * 
            FROM exons_to_domains_data 
            WHERE "Pfam ID"=:pfam_id 
            ORDER BY "Exon rank in transcript"
            """
    tdata = pd.read_sql_query(sql=text(query), con=engine, params={'pfam_id': pfam_id})

    # OLD
    data_old = pd.read_csv('../data/final.csv')
    tdata_old = data_old[data_old['Pfam ID'].isin([pfam_id])].sort_values(by=['Exon rank in transcript'])

    # COMPARE
    assert (np.array_equal(tdata['Exon stable ID'].sort_values().values, tdata_old['Exon stable ID'].sort_values().values))

if __name__ == '__main__':
    transcript_with_domain()
