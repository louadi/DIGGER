import os
import pandas as pd
import pickle


def load_obj(name):
    with open('domain/data/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


PPI_all = {}
g2d_all = {}
protein_all = {}
DDI_all = {}
DomainG_all = {}
gid2name_all = {}


for organism in os.listdir('domain/data'):
    if not os.path.isdir('domain/data/' + organism):
        continue
    trivial_name = organism.split("[")[1][:-1]
    DomainG_all[trivial_name] = load_obj(organism + '/DomainG')
    gid2name_all[trivial_name] = load_obj(organism + '/gid2name')
    PPI_all[trivial_name] = load_obj(organism + '/PPI')
    g2d_all[trivial_name] = load_obj(organism + '/g2d')
    protein_all[trivial_name] = pd.read_csv('domain/data/' + organism + '/all_Proteins.csv')
    if os.path.isfile('domain/data/' + organism + '/DDI.pkl'):
        DDI_all[trivial_name] = load_obj(organism + '/DDI')