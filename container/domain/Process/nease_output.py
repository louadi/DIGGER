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
from domain.nease.process import webify_table
from django.conf import settings
import uuid

images_path = os.path.join(settings.MEDIA_ROOT, 'images/')
data_path = os.path.join(settings.MEDIA_ROOT, 'nease_tables/')
nease_path = 'nease_events/'
# The subdirectories contain files saved for one week, one month, and six months.
# To be very lenient, we calculated every month with 31 days.
days_to_folder = {"0": nease_path+"zero_days/", "7": nease_path+"seven_days/", "31": nease_path+"thirtyone_days/",
                  "186": nease_path+"onehundredeightysix_days/"}
default_path = days_to_folder["7"]

for path in [images_path, data_path] + list(days_to_folder.values()):
    if not os.path.exists(path):
        os.makedirs(path)

# web_tables_options is a dictionary that contains the options for webifying the tables
web_tables_options = {
    'domains': {'link_col': ['Gene name', 'Gene stable ID', 'Exon stable ID', 'Pfam ID'],
                'link_prefix': ['https://www.ncbi.nlm.nih.gov/gene/', 'https://www.ensembl.org/id/',
                                'https://www.ensembl.org/id/', 'https://www.ebi.ac.uk/interpro/entry/pfam/'],
                'link_id': ['NCBI gene ID', None, None, None],
                'drop_col': ['NCBI gene ID']},
    'edges': {'link_col': ['Gene name', 'Affected binding'],
              'link_prefix': ['https://www.ncbi.nlm.nih.gov/gene/', 'https://www.ncbi.nlm.nih.gov/gene/'],
              'link_id': ['NCBI gene ID', 'Affected binding (NCBI)'],
              'drop_col': ['NCBI gene ID', 'Affected binding (NCBI)']},
    'elm': {'link_col': ['Gene name', 'Gene stable ID', 'ELMIdentifier'],
            'link_prefix': ['https://www.ncbi.nlm.nih.gov/gene/', 'https://www.ensembl.org/id/',
                            'http://elm.eu.org/elms/'],
            'link_id': ['entrezgene', None, None],
            'drop_col': ['entrezgene', 'ELM link']},
    'pdb': {'link_col': ['Gene name', 'Gene stable ID', 'Co-resolved interactions symbol'],
            'link_prefix': ['https://www.ncbi.nlm.nih.gov/gene/', 'https://www.ensembl.org/id/',
                            'https://www.ncbi.nlm.nih.gov/gene/'],
            'link_id': ['NCBI gene ID', None, 'Co-resolved interactions'],
            'drop_col': ['NCBI gene ID', 'Co-resolved interactions']},
    'pathway': {'link_col': ['Spliced genes', 'Affected binding (edges)'],
                'link_prefix': ['https://www.ncbi.nlm.nih.gov/gene/', 'https://www.ncbi.nlm.nih.gov/gene/'],
                'link_id': ['NCBI gene ID', 'Affected binding (NCBI)'],
                'drop_col': ['NCBI gene ID', 'Affected binding (NCBI)']}
}


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

    domains: pd.DataFrame = events.get_domains()
    # check if domains is not empty
    if not domains.empty:
        domains.to_csv(f"{data_path}{run_id}_domains.csv", index=False)
        domains = webify_table(domains, web_tables_options['domains'])

    edges = events.get_edges()
    if not edges.empty:
        edges.to_csv(f"{data_path}{run_id}_edges.csv", index=False)
        edges = webify_table(edges, web_tables_options['edges'])

    info_tables = {'domains': domains, 'edges': edges}

    if not params.get('only_ddis', False):
        elm = events.get_elm()
        if not elm.empty:
            elm.to_csv(f"{data_path}{run_id}_elm.csv", index=False)
            elm = webify_table(elm, web_tables_options['elm'])

        pdb = events.get_pdb()
        if not pdb.empty:
            pdb.to_csv(f"{data_path}{run_id}_pdb.csv", index=False)
            pdb = webify_table(pdb, web_tables_options['pdb'])

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
    info_tables = {}
    try:
        domains = webify_table(pd.read_csv(f"{data_path}{run_id}_domains.csv"), web_tables_options['domains'])
        edges = webify_table(pd.read_csv(f"{data_path}{run_id}_edges.csv"), web_tables_options['edges'])

        if os.path.exists(f"{data_path}{run_id}_elm.csv"):
            elm = webify_table(pd.read_csv(f"{data_path}{run_id}_elm.csv"), web_tables_options['elm'])
        else:
            elm = pd.DataFrame(columns=["Gene name", "Gene stable ID", "ELMIdentifier", "dPSI"])

        if os.path.exists(f"{data_path}{run_id}_pdb.csv"):
            pdb = webify_table(pd.read_csv(f"{data_path}{run_id}_pdb.csv"), web_tables_options['pdb'])
        else:
            pdb = pd.DataFrame(columns=["Gene name", "Gene stable ID", "Co-resolved interactions symbol"])

        info_tables = {'domains': domains, 'edges': edges, 'elm': elm, 'pdb': pdb}
    except Exception as e:
        traceback.print_exc()
        print(e)
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

    create_plot(terms, p_values, cut_off, f"{images_path}{run_id}_clenr")

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

    create_plot(terms, pvalues, cut_off, f"{images_path}{run_id}_neenr")

    return enrich_table


def pathway_info(events, pathway, run_id):
    events, _ = events
    try:
        pathway_info_table = events.path_analysis(pathway)
        if not isinstance(pathway_info_table, pd.DataFrame) or len(pathway_info_table) == 0:
            return pathway_info_table
        print("Table is:", len(pathway_info_table))
        pathway_info_table.to_csv(f"{data_path}{run_id}_path_{pathway}.csv")
        pathway_info_table = webify_table(pathway_info_table, web_tables_options['pathway'])
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
    # add line break to terms that are longer than 50 characters
    terms = cut_long_terms(terms)

    colors = ['#344552', '#355871', '#366B91', '#377DB0', '#448FC5', '#5E9FCD', '#78B0D5', '#92C0DD']
    plt.figure(figsize=(12, 5))
    plt.grid(color='#D2D2D2', linestyle='-', zorder=0)

    plt.barh(terms[::-1], pvalues[::-1], color=colors[::-1], zorder=3)
    plt.axvline(x=cut_off, color='red', linestyle='--', linewidth=1.5, zorder=5)
    plt.ylabel("Terms", fontsize=12)
    plt.xlabel("-log10(adjusted p-value)", fontsize=12)

    plt.legend(['Cut-off', 'Adjusted p-value'], fontsize=10, frameon=True)

    # Edit the spines for a cleaner look
    ax = plt.gca()
    ax.tick_params(bottom=False)
    ax.tick_params(left=False)
    ax.spines[['left', 'bottom']].set_color('#D2D2D2')
    ax.spines[['left', 'bottom']].set_linewidth(2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    plt.savefig(f"{filename}_thumb.jpg", bbox_inches='tight')
    plt.savefig(f"{filename}.jpg", bbox_inches='tight', dpi=1200)
    # flush the plot
    plt.clf()
    plt.close()


def cut_long_terms(terms):
    cut_terms = []
    for term in terms:
        if len(term) > 50:
            # find the first space after the 50th character
            space_index = term.find(' ', 50)
            if space_index == -1:
                cut_terms.append(term)
            else:
                cut_terms.append(term[:space_index] + '\n' + term[space_index + 1:])
        else:
            cut_terms.append(term)
    return cut_terms


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


def match_name_with_format(filename):
    name_matches = {'deltapsi': 'MAIJQ',
                    'maijq': 'MAIJQ',
                    'diff': 'Whippet',
                    'whippet': 'Whippet',
                    'mats': 'rmats'}

    for name, format in name_matches.items():
        if name in filename.lower():
            return f"{format} input"
    return None
