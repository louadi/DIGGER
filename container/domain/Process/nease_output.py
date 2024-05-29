import pandas as pd
import os
from domain.nease import nease
from django.conf import settings
import uuid

images_path = os.path.join(settings.MEDIA_ROOT, 'images/')
if not os.path.exists(images_path):
    os.makedirs(images_path)


def run_nease(data, organism, params):
    run_id = str(uuid.uuid4())
    image_path = images_path + run_id
    print(image_path)

    events = nease.run(data, organism, params.get('db_type', []),
                       params.get('p_value', 0.05),
                       params.get('rm_not_in_frame'),
                       params.get('divisible_by_3'),
                       params.get('min_delta', 0.1),
                       params.get('majiq_confidence', 0.95),
                       params.get('only_ddis', False),
                       params.get('confidences', []))
    print(events)
    try:
        classic_enrich_table = events.classic_enrich(params.get('enrich_dbs', []))
        classic_enrich_table['Genes'] = classic_enrich_table['Genes'].apply(lambda x: x.replace(';', ', '))
    except ValueError:
        classic_enrich_table = pd.DataFrame(
            columns=["Gene_set", "Term", "Overlap", "P-value", "Adjusted P-value", "Old P-value",
                     "Old Adjusted P-value", "Odds Ratio", "Combined Score", "Genes"])

    events.get_stats(file_path=image_path)

    return events, run_id, classic_enrich_table
