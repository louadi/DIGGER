import pickle

import pandas as pd
import os
from domain.nease import nease
from django.conf import settings
import uuid

images_path = os.path.join(settings.MEDIA_ROOT, 'images/')
if not os.path.exists(images_path):
    os.makedirs(images_path)

if not os.path.exists('nease_events/'):
    os.makedirs('nease_events/')

nease_events = 'nease_events/'


def same_input(run_id):
    pass


def run_nease(data, organism, params):
    run_id = str(uuid.uuid4())
    image_path = images_path + run_id

    events = nease.run(data, organism, params.get('db_type', []),
                       params.get('p_value', 0.05),
                       params.get('rm_not_in_frame'),
                       params.get('divisible_by_3'),
                       params.get('min_delta', 0.1),
                       params.get('majiq_confidence', 0.95),
                       params.get('only_ddis', False),
                       params.get('confidences', []))

    events.get_stats(file_path=image_path)

    # save events to pickle
    events.save(nease_events + run_id)
    return events, run_id


def get_nease_events(run_id):
    events = nease.load(nease_events + run_id + '.pkl')
    return events


def nease_domains(events):
    return events.get_domains()


def nease_classic_enrich(events, databases):
    try:
        classic_enrich_table = events.classic_enrich(databases)
        classic_enrich_table['Genes'] = classic_enrich_table['Genes'].apply(lambda x: x.replace(';', ', '))
    except ValueError:
        classic_enrich_table = pd.DataFrame(
            columns=["Gene_set", "Term", "Overlap", "P-value", "Adjusted P-value", "Old P-value",
                     "Old Adjusted P-value", "Odds Ratio", "Combined Score", "Genes"])
    return classic_enrich_table
