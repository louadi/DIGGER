import json
import pickle
import traceback
from io import StringIO

import numpy as np
import pandas as pd
import os

from matplotlib import pyplot as plt

from domain.models import NeaseSaveLocationMapping
from domain.nease import nease
from django.conf import settings
import uuid

images_path = os.path.join(settings.MEDIA_ROOT, 'images/')
data_path = os.path.join(settings.MEDIA_ROOT, 'nease_tables/')
nease_path = 'nease_events/'
# The subdirectories contain files saved for one week, one month, and six months. To be very lenient, we calculated every month with 31 days.
days_to_folder = {"7": nease_path+"seven_days/", "31": nease_path+"thirtyone_days/", "186": nease_path+"onehundredeightysix_days/"}
default_path = days_to_folder["7"]

for path in [images_path, data_path] + list(days_to_folder.values()):
    if not os.path.exists(path):
        os.makedirs(path)


def run_nease(data, organism, params, file_name='', custom_name=''):
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

    # TODO: add function to better format tables for web view while keeping the original tables in csv
    domains: pd.DataFrame = events.get_domains()
    # check if domains is not empty
    if not domains.empty:
        domains.to_csv(f"{data_path}{run_id}_domains.csv")

    edges = events.get_edges()
    if not edges.empty:
        edges.to_csv(f"{data_path}{run_id}_edges.csv")

    info_tables = {'domains': domains, 'edges': edges}

    if not params.get('only_ddis', False):
        elm = events.get_elm()
        if not elm.empty:
            elm.to_csv(f"{data_path}{run_id}_elm.csv")
        pdb = events.get_pdb()
        if not pdb.empty:
            pdb.to_csv(f"{data_path}{run_id}_pdb.csv")
        info_tables.update({'elm': elm, 'pdb': pdb})

    # remove the unnamed column
    for key, value in info_tables.items():
        if 'Unnamed: 0' in value.columns:
            value.drop(columns=['Unnamed: 0'], inplace=True)

    # save events to pickle
    events.save(default_path + run_id)
    NeaseSaveLocationMapping(run_id=run_id, saved_for_days=7, file_name=file_name, custom_name=custom_name).save()
    return events, info_tables, run_id


def file_needs_cleaning(file_obj):
    file_obj.seek(0)
    for line in file_obj:
        if line.endswith(b'\t\n') or line.endswith(b'\t \n'):
            file_obj.seek(0)
            print(f"File {file_obj.name} needs cleaning.")
            return True
    # Reset the file pointer to the beginning
    file_obj.seek(0)
    return False


def read_extra_spaces(file_obj):
    cleaned_lines = []
    for line in file_obj:
        # Decode the line if it's bytes, strip trailing tab and space, and encode it back to bytes
        cleaned_line = line.decode('utf-8').rstrip('\t').rstrip()
        cleaned_lines.append(cleaned_line)

    cleaned_data = '\n'.join(cleaned_lines)
    df = pd.read_table(StringIO(cleaned_data))
    return df


def get_nease_events(run_id):
    days = NeaseSaveLocationMapping.get_saved_for_days(run_id)
    if days not in days_to_folder:
        file_path = default_path
    else:
        file_path = days_to_folder[str(days)]
    print(f"Loading events from {file_path + run_id + '.pkl'}")
    events = nease.load(file_path + run_id + '.pkl')
    domains = pd.read_csv(f"{data_path}{run_id}_domains.csv")
    edges = pd.read_csv(f"{data_path}{run_id}_edges.csv")
    info_tables = {'domains': domains, 'edges': edges}
    if os.path.exists(f"{data_path}{run_id}_elm.csv"):
        elm = pd.read_csv(f"{data_path}{run_id}_elm.csv")
        pdb = pd.read_csv(f"{data_path}{run_id}_pdb.csv")
        info_tables.update({'elm': elm, 'pdb': pdb})
    return events, info_tables


def nease_domains(events):
    return events.get_domains()


# this function got cut from the final version
def nease_classic_enrich(events, databases, run_id):
    events, _ = events
    try:
        classic_enrich_table = events.classic_enrich(databases, cutoff=events.get_p_value())
        classic_enrich_table['Genes'] = classic_enrich_table['Genes'].apply(lambda x: x.replace(';', ', '))
        classic_enrich_table.to_csv(f"{data_path}{run_id}_clenr.csv")
    except ValueError:
        classic_enrich_table = pd.DataFrame(
            columns=["Gene_set", "Term", "Overlap", "P-value", "Adjusted P-value", "Old P-value",
                     "Old Adjusted P-value", "Odds Ratio", "Combined Score", "Genes"])
        return classic_enrich_table

    # make plot displaying the most enriched terms
    terms = classic_enrich_table['Term'][:8]
    p_values = classic_enrich_table['Adjusted P-value'][:8]
    p_values = [-np.log10(x) for x in p_values]
    cut_off = -np.log10(events.get_p_value())

    create_plot(terms, p_values, cut_off, f"{images_path}{run_id}_clenr.jpg")

    return classic_enrich_table


def nease_enrichment(events, databases, run_id):
    events, _ = events
    try:
        enrich_table = events.enrich(databases)
        enrich_table.to_csv(f"{data_path}{run_id}_neenr.csv")
    except ValueError:
        enrich_table = pd.DataFrame(
            columns=["Gene_set", "Term", "Overlap", "P-value", "Adjusted P-value", "Old P-value",
                     "Old Adjusted P-value", "Odds Ratio", "Combined Score", "Genes"])
        return enrich_table
    enrich_table = enrich_table.sort_values(by='Nease score', ascending=False)
    terms = enrich_table['Pathway name'][:8]
    pvalues = enrich_table['adj p_value'][:8]
    pvalues = [-np.log10(x) for x in pvalues]
    cut_off = -np.log10(events.get_p_value())

    enrich_table.iloc[:, 0] = enrich_table.iloc[:, 0].apply(lambda x: f"""<p class="tooltips" 
    tooltip="<p class='my-1'>Inspect {x}:</p> 
    <button class='btn btn-secondary' onclick=&quot;analysePathway('{x}')&quot;>Analyse</button>
    <button class='btn btn-secondary' onclick=&quot;visualisePath('{x}', 0.8)&quot;>Visualise</button>"
    tooltip-position="right">{x}</p>""".replace("\n", ""))

    create_plot(terms, pvalues, cut_off, f"{images_path}{run_id}_neenr.jpg")

    return enrich_table


def pathway_info(events, pathway, run_id):
    events, _ = events
    try:
        pathway_info_table = events.path_analysis(pathway)
        if not isinstance(pathway_info_table, pd.DataFrame) or len(pathway_info_table) == 0:
            return pathway_info_table
        print("Table is:", len(pathway_info_table))
        pathway_info_table.to_csv(f"{data_path}{run_id}_path_{pathway}.csv")
    except ValueError as e:
        traceback.print_exc()
        pathway_info_table = {'errorMsg': str(e)}
        # pathway_info_table = pd.DataFrame(
        #     columns=["Spliced genes", "NCBI gene ID", "Gene is known to be in the pathway",
        #              "Percentage of edges associated to the pathway", "p_value", "Affected binding (edges)",
        #              "Affected binding (NCBI)"])
    return pathway_info_table


def visualise_path(events, pathway, k):
    events, _ = events
    try:
        pathway_visualisation = events.Vis_path(pathway, k=k)
    # except ValueError as e:
    #     print(e)
    #     pathway_visualisation = {'error_msg': 'Something went wrong while visualising the pathway'}
    except Exception as e:
        traceback.print_exc()
        # add only the message of the error, not the error itself
        pathway_visualisation = {'errorMsg': str(e)}
    return pathway_visualisation


def create_plot(terms, pvalues, cut_off, filename):
    plt.style.use('ggplot')
    plt.barh(terms[::-1], pvalues[::-1])
    plt.axvline(x=cut_off, color='b', linestyle='--')
    plt.xlabel('-log10(adjusted p-value)')
    plt.ylabel('Terms')
    # explain the red line in the legend
    plt.legend(['Cut-off', 'Adjusted p-value'])
    plt.savefig(filename, bbox_inches='tight')
    # flush the plot
    plt.clf()
    plt.close()

def change_save_timing(run_id, days):
    mapping = NeaseSaveLocationMapping.objects.get(run_id=run_id)
    current_days_folder = mapping.get_number_of_saved_for_days()
    if days not in days_to_folder:
        new_file_path = default_path
    else:
        new_file_path = days_to_folder[str(days)]
    if current_days_folder not in days_to_folder:
        old_file_path = default_path
    else:
        old_file_path = days_to_folder[str(current_days_folder)]
    # move the file
    os.rename(old_file_path + run_id + '.pkl', new_file_path + run_id + '.pkl')
    # update the database
    mapping.saved_for_days = int(days)
    mapping.save()

    return json.dumps(
        {"logmessage": "Changing the save timing from " + str(current_days_folder) + " to " + str(days) + " was successful.",
         "days_left": mapping.days_left()}
    )