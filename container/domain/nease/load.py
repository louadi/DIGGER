import pandas as pd
import pickle5 as pickle
import os

pd.options.mode.chained_assignment = None  # default='warn'


# Helper function

def load_obj(data_folder):
    data_folder = os.path.join(os.path.dirname(__file__), data_folder)
    with open(data_folder + '.pkl', 'rb') as f:
        return pickle.load(f)


# this project uses python 3.6 and only supports pickle protocol 4, whereas we need 5
def load_df(path):
    with open(path, 'rb') as f:
        return pickle.load(f)


# Save data as dictionary        
# Databases
database_mapping = {}
Pathways = {}
# Join graph
network = {}
# The PPI
PPI = {}

elm = {}
elm_interactions = {}

pdb = {}

# hierarchy for the pathways (currently only reactome)
pathway_hierarchy = {}

here = os.path.dirname(__file__)

# add more species here
for species in ["Homo sapiens[human]", "Mus musculus[mouse]"]:
    name = species.split("[")[1][:-1]
    network[name] = load_obj(os.path.join(here, f"data/{species}/graph"))
    PPI[name] = load_obj(os.path.join(here, f"data/{species}/PPI"))
    non_coding = load_obj(os.path.join(here, f'data/{species}/non_coding'))
    database_mapping[name] = load_df(os.path.join(here, f"data/{species}/{name.capitalize()}"))
    Pathways[name] = load_df(os.path.join(here, f"data/{species}/pathways"))
    elm[name] = load_df(os.path.join(here, f"data/{species}/elm"))
    elm_interactions[name] = load_df(os.path.join(here, f"data/{species}/ELM_interactions"))
    pdb[name] = load_df(os.path.join(here, f"data/{species}/pdb"))
    pathway_hierarchy[name] = {}
    pathway_hierarchy[name]['Reactome'] = load_obj(os.path.join(here, f"data/{species}/pathway_hierarchy_reactome"))
    pathway_hierarchy[name]['KEGG'] = load_obj(os.path.join(here, f"data/{species}/pathway_hierarchy_kegg"))

