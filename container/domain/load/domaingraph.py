import csv
import networkx as nx
import pickle
from operator import itemgetter



def load_obj(name):
    with open('domain/data/'+name + '.pkl', 'rb') as f:
        return pickle.load(f)


        