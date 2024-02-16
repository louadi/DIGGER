import os
import pickle
import random
import re
import timeit

import pandas as pd

from django.db import connection
from django.conf import settings
from django.http import HttpResponse, JsonResponse
from django.shortcuts import render, redirect
from django.utils.html import escape
from io import StringIO

from pandas.errors import ParserError

from .Process import exonstodomain as exd
from .Process import exon as ex
from .Process import process_data as pr
from .Process import transcript as tr
from .Process import gene as g
from .Process import network_analysis as nt

# --- Create folder
# Global jobs path
jobs_path = os.path.join(settings.MEDIA_ROOT, 'jobs')
if not os.path.exists(jobs_path):
    os.makedirs(jobs_path)



# Display transcripts of a gene
def gene(request, gene_ID, organism):
    transcript_table, gene_name = g.input_gene(gene_ID, organism)

    if transcript_table == []:
        return HttpResponse(' wrong entry or protein without any known Pfam domains')

    context = {
        'tb': transcript_table,
        'name': gene_name,
    }
    return render(request, 'visualization/gene.html', context)


# Display information of an exon (Exon-Level Analysis)
def exon(request, organism, exon_ID):
    v = ex.input_exon(exon_ID, organism)

    if v == True:
        return HttpResponse(' wrong entry or exon in a gene without any known Pfam domains')
    else:
        _, domains, gene_name, Ensemble_geneID, entrezID, tb_transc, table_domains, number = v

    # only if the exon code for domains with known interactions

    nodes_domainV = {}
    edges_domainV = {}
    switcher = []
    switcher_js = []
    first = []
    maxx = 0

    # Interactionview
    Interactiveview_select = {}
    Interactiveview_switch = []
    first_victim = []

    if number > 0:

        # ProteinView
        nodes, edges, pd_interaction = ex.vis_exon(domains, entrezID, gene_name, exon_ID, organism)

        # DomainView
        first = domains[0]
        for pfams in domains:
            n, e, _, _ = exd.vis_node_(entrezID + "." + pfams, organism)
            if len(e['original']) > maxx:
                maxx = len(e)
                first = pfams
            if len([val for sublist in e.values() for val in sublist]) != 0:
                # initialise keys with empty list if they don't exist
                for key in n.keys():
                    if key not in nodes_domainV:
                        nodes_domainV[key] = []
                for key in e.keys():
                    if key not in edges_domainV:
                        edges_domainV[key] = []

                nodes_domainV = {k: nodes_domainV[k] + v for k, v in n.items()}
                edges_domainV = {k: edges_domainV[k] + v for k, v in e.items()}
                switcher.append('<option value="' + pfams + '"> ' + pfams + '</option>')
                switcher_js.append('case "' + pfams + '": return node.source === "' + pfams + '";')

    else:
        nodes, edges, pd_interaction = {}, {}, []

    # PPI res interfaces on the exon:
    # table: HTML table with all PPIs that have res interface mapped to the exon
    # number_of_PPI: number of interactions

    table, number_of_PPI = ex.PPI_inter(exon_ID, gene_name, organism)

    if number > 0 and len(pd_interaction) > 0:
        # added to combine evidence of DDI and Residue in one final table
        if number_of_PPI > 0:
            ppi_from_res = table['Partner Protein'].unique()
            f = pd_interaction['Partner Protein'].isin(ppi_from_res)
            pd_interaction.loc[f, 'Residue evidence'] = '<center>&#9989;</center>'

        # InteractionView
        pd_interaction['_'] = 'Interaction with &nbsp;&nbsp;' + pd_interaction[
            "Partner Protein"] + '&nbsp;&nbsp; ( Score ' + pd_interaction["Score"].round(2).astype(str) + ' )'

        pd_interaction['selector'] = '<option value="' + pd_interaction['NCBI gene ID'].astype(str) + '"> ' + \
                                     pd_interaction['_'] + '</option>'

        pd_interaction['switcher'] = 'case "' + pd_interaction['NCBI gene ID'].astype(
            str) + '": return (node.id === "' + pd_interaction['NCBI gene ID'].astype(
            str) + '")   || (node.id === "' + entrezID + '")     ||  (node.origin==="' + pd_interaction[
                                         'NCBI gene ID'].astype(str) + '") ||  (node.origin==="' + entrezID + '")  ;'

        Interactiveview_select['original'] = pd_interaction[pd_interaction['Origin'] == 'original']['selector'].tolist()
        Interactiveview_select['predicted'] = pd_interaction[pd_interaction['Origin'] == 'predicted']['selector'].tolist()
        Interactiveview_switch = pd_interaction['switcher'].tolist()
        # the first protein to show
        first_victim = pd_interaction['NCBI gene ID'].tolist()[0]

        pd_interaction = pd_interaction[
            ["Affected Protein", 'Partner Protein', 'NCBI gene ID', 'Retained DDIs', 'Lost DDIs',
             'Percentage of lost domain-domain interactions', 'Residue evidence', "Protein-protein interaction",
             'Score', 'Origin']]

        pd_interaction = pd_interaction.rename(columns={

            "Percentage of lost domain-domain interactions": "% of missing DDIs",
            "Retained DDIs": "Retained domain-domain interactions",
            "Lost DDIs": "Missing domain-domain interactions",
            "Protein-protein interaction": "Protein-protein interaction",
            'Residue evidence': 'Residue-level evidence*'
        })

        pd_interaction = pd_interaction.to_html(table_id='Interaction_table', **settings.TO_HTML_RESPONSIVE_PARAMETERS)

    table = table.to_html(**settings.TO_HTML_PARAMETERS)

    context = {

        'tb1': tb_transc,
        'tb2': table_domains,
        'tb3': pd_interaction,
        'tb4': table,
        'name': gene_name,
        'exon_ID': exon_ID,
        'entrezID': entrezID,
        'gID': Ensemble_geneID,
        'dis': number > 0,

        "dis2": number == -1,
        'dis3': number_of_PPI != 0,

        # only a self loop for the domains>> no interactionView
        'dis4': number > 0 and len(pd_interaction) == 0,
        'long_table': number_of_PPI > 25,
        'pv_nodes': nodes,
        'pv_edges': edges,

        'first_domain': first,
        'switch1': switcher,
        'switch2': switcher_js,
        'Domainview_edges': edges_domainV,
        'Domainview_nodes': nodes_domainV,

        'Interactiveview_selec': Interactiveview_select,
        'first_vict': first_victim,
        'Interactiveview_switch': Interactiveview_switch,

        'enable_Proteinview': len(edges_domainV) > 70,
    }
    return render(request, 'visualization/exon.html', context)


# Display information of a transcript or a protein (Isoform-Level Analysis)
def transcript(request, P_id, organism):
    print("Currently in transcript view")

    out = tr.Protein_view(P_id, organism)

    if out == 0: return HttpResponse(' Wrong entry or protein without any known Pfam domains')
    if out == 1: return HttpResponse(' The selected protein does not have any interaction in the current PPI database')

    nodes, edges, unique, text1, tran_name, Ensemble_geneID, entrezID, gene_description, droped1, droped2, trID, p, \
        missed, pd_interaction, isoforms, co_partners = out

    # Interactionview
    Interactiveview_select = {}
    Interactiveview_switch = []
    first_victim = []

    if len(pd_interaction) != 0:
        pd_interaction['Residue evidence'] = ''

        pd_interaction.loc[
            pd_interaction["NCBI gene ID"].isin(co_partners), 'Residue evidence'] = '<span>&#9733;</span>'

        pd_interaction = pd_interaction.sort_values('Protein name')
        pd_interaction['_'] = '&nbsp;Interaction &nbsp; with &nbsp;' + pd_interaction[
            "Protein name"] + '&nbsp;&nbsp; ( Score ' + pd_interaction["Score"].round(2).astype(str) + ')&nbsp;&nbsp;' + \
                              pd_interaction['Residue evidence'] + '&nbsp;&nbsp;'

        pd_interaction['selector'] = '<option value="' + pd_interaction['NCBI gene ID'].astype(str) + '"> ' + \
                                     pd_interaction['_'] + '</option>'

        pd_interaction['switcher'] = 'case "' + pd_interaction['NCBI gene ID'].astype(
            str) + '": return (node.id === "' + pd_interaction['NCBI gene ID'].astype(
            str) + '")   || (node.id === "' + entrezID + '")     ||  (node.origin==="' + pd_interaction[
                                         'NCBI gene ID'].astype(str) + '") ||  (node.origin==="' + entrezID + '")  ;'

        Interactiveview_select['original'] = pd_interaction[pd_interaction['Origin'] == 'original']['selector'].tolist()
        Interactiveview_select['predicted'] = pd_interaction[pd_interaction['Origin'] == 'predicted']['selector'].tolist()
        Interactiveview_switch = pd_interaction['switcher'].tolist()
        # the first protein to show
        first_victim = pd_interaction['NCBI gene ID'].tolist()[0]

        pd_interaction = pd_interaction.rename(columns={

            "Protein name": "Partner Protein",
            'Residue evidence': 'Residue-level evidence',
            "Percentage of lost domain-domain interactions": "% of missing DDIs",
            "Retained DDIs": "Retained domain-domain interactions",
            "Lost DDIs": "Missing domain-domain interactions",
            "Protein-protein interaction": "Protein-protein interaction"
        })

        pd_interaction = pd_interaction[
            ["Selected Protein variant", 'Partner Protein', 'NCBI gene ID', 'Retained domain-domain interactions',
             'Missing domain-domain interactions', '% of missing DDIs', 'Residue-level evidence',
             "Protein-protein interaction", 'Score', 'Origin']]
        pd_interaction = pd_interaction.to_html(table_id='Interaction_table', **settings.TO_HTML_RESPONSIVE_PARAMETERS)

    # Get ID of missing domains with interactions
    if len(missed) != 0:
        missing_domains = missed['Pfam ID'].unique()
        missed = missed.to_html(**settings.TO_HTML_PARAMETERS)

    nodes_domainV = {}
    edges_domainV = {}
    switcher = []
    switcher_js = []
    first = unique[0]
    maxx = 0

    # DomainView for retained domains
    for pfams in unique:
        n, e, _, _ = exd.vis_node_(entrezID + "." + pfams, organism)
        if len(e['original']) > maxx:
            maxx = len(e)
            first = pfams
        if len([val for sublist in e.values() for val in sublist]) != 0:
            # initialise keys with empty list if they don't exist
            for key in n.keys():
                if key not in nodes_domainV:
                    nodes_domainV[key] = []
            for key in e.keys():
                if key not in edges_domainV:
                    edges_domainV[key] = []

            nodes_domainV = {k: nodes_domainV[k] + v for k, v in n.items()}
            edges_domainV = {k: edges_domainV[k] + v for k, v in e.items()}

            switcher.append('<option value="' + pfams + '"> ' + pfams + '</option>')
            switcher_js.append('case "' + pfams + '": return node.source === "' + pfams + '";')

    # DomainView for missing domains

    switcher_m = []
    if len(missed) != 0:
        for pfams in missing_domains:
            n, e, _, _ = exd.vis_node_(entrezID + "." + pfams, organism)
            if len(e['original']) > maxx:
                maxx = len(e)
                first = pfams
            if len([val for sublist in e.values() for val in sublist]) != 0:
                # initialise keys with empty list if they don't exist
                for key in n.keys():
                    if key not in nodes_domainV:
                        nodes_domainV[key] = []
                for key in e.keys():
                    if key not in edges_domainV:
                        edges_domainV[key] = []

                nodes_domainV = {k: nodes_domainV[k] + v for k, v in n.items()}
                edges_domainV = {k: edges_domainV[k] + v for k, v in e.items()}
                switcher_m.append('<option value="' + pfams + '"> ' + pfams + ' (missing in the isoform) </option>')
                switcher_js.append('case "' + pfams + '": return node.source === "' + pfams + '";')

    # get total length of edges_domainV to determine disable_Proteinview
    len_edges_domainV = len([val for sublist in edges_domainV.values() for val in sublist])

    context = {
        'dt': droped1,
        'text1': text1,
        'tran_name': tran_name,
        'gene_description': gene_description,
        'trID': trID,
        'gID': Ensemble_geneID,
        'path': p,
        'pv_nodes': nodes,
        'pv_edges': edges,
        'entrezID': entrezID,
        'dt2': droped2,
        'dt3': missed,
        'dt4': pd_interaction,
        "dt5": isoforms,

        'dis1': missed != [],
        'dis2': pd_interaction != [],

        'dis3': isoforms != [],

        'Interactiveview_select': Interactiveview_select,
        'first_vict': first_victim,
        'Interactiveview_switch': Interactiveview_switch,

        'first_domain': first,
        'switch1': switcher,
        'switch1_missing': switcher_m,
        'switch2': switcher_js,
        'Domainview_edges': edges_domainV,
        'Domainview_nodes': nodes_domainV,

        # make it make sense
        'disable_Proteinview': (len_edges_domainV > 200) or (len_edges_domainV > 130 and len(unique) + len(missed) == 1),
    }

    return render(request, 'visualization/transcript.html', context)


# isoform_level search box
def isoform_level(request):
    if "search" in request.GET:  # If the form is submitted
        # Get and sanitize the search_query
        search_query = request.GET['search'].strip()

        # Get organism (human, mouse)
        organism = request.GET.get('organism', None)

        # Try and parse the search_query as gene name from the database
        with connection.cursor() as cursor:
            cursor.execute("""SELECT ensembl_id FROM domain_gene_""" + organism + """ WHERE gene_symbol ILIKE %s""",
                           [search_query])
            row = cursor.fetchone()

        if row:
            search_query = row[0]

        search_query = search_query.split("+")[0]
        search_query = search_query.split("%")[0]
        search_query = search_query.split(".")[0]
        # regex to check if search query starts with ENS and ends with T or P
        if re.match(r'ENS\w*[T,P]\d+$', search_query):
            print("User input is a protein")
            return redirect(transcript, P_id=search_query, organism=organism)

        elif re.match(r'ENS\w*G\d+$', search_query):
            print("User input is a gene")
            return redirect(gene, gene_ID=search_query, organism=organism)

    return render(request, 'setup/isoform_level.html', )


# exon_level search box
def exon_level(request):
    if "search" in request.GET:  # If the form is submitted
        # Input and Exon ID
        print('-----------------------------------------------------------')
        search_query = request.GET['search']
        organism = request.GET.get('organism', None)
        search_query = search_query.replace(" ", "")
        search_query = search_query.split("+")[0]
        search_query = search_query.split("%")[0]
        search_query = search_query.split(".")[0]
        print(search_query)

        if re.match(r'ENS\w*[E]\d+$', search_query):
            return redirect(exon, organism=organism, exon_ID=search_query)

    if "search 2" in request.GET:  # If the form is submitted
        # Input coordinate of the exon
        # Check if coordinate are correct
        #  Example   ' ENSG00000266028  206437964 206437042 '

        print('-----------------------------------------------------------')
        search_query = request.GET['search 2']
        organism = request.GET.get('organism', None)

        search_query = search_query.split(" ")
        search_query = [x for x in search_query if x != '']
        # search_query[0]=search_query[0].split(".")[0]
        if len(search_query) == 3 and re.match(r'ENS\w*[G]\d+$', search_query[0]) and \
                search_query[1].isdigit() and search_query[2].isdigit():

            gene_ID = search_query[0]
            s1 = int(search_query[1])
            e1 = int(search_query[2])

            # Correct for very big inputs
            if abs(s1 - e1) < 3000:

                # map coordinates to exon
                exonID = pr.coordinate_to_exonID(gene_ID, s1, e1, organism)

                if exonID != []:
                    return redirect(exon, organism=organism, exon_ID=exonID)
                    # return exon(request,exonID)
                else:
                    return HttpResponse("<h1>No match</h1>")
    if "search 3" in request.GET:  # If option 3 is selected

        # Get and sanitize the search_query
        search_query = request.GET['search 3'].strip()
        organism = request.GET.get('organism', None)

        with connection.cursor() as cursor:
            cursor.execute("""SELECT ensembl_id FROM domain_gene_""" + organism + """ WHERE gene_symbol ILIKE %s""",
                           [search_query])
            row = cursor.fetchone()

        if row:
            search_query = row[0]

        search_query = search_query.split("+")[0]
        search_query = search_query.split("%")[0]
        search_query = search_query.split(".")[0]

        # Input search is a protein:
        if re.match(r'ENS\w*[TP]\d+$', search_query):
            return redirect(transcript, organism=organism, P_id=search_query)

        # Input search is a gene:
        elif re.match(r'ENS\w*G\d+$', search_query):
            return redirect(gene, organism=organism, gene_ID=search_query)

    return render(request, 'setup/exon_level.html', )


# PPI network analysis
def network(request):
    error_message = ""
    jump_div = ""

    # Option 1: List of Ensembl IDs
    if "option1" in request.POST:
        organism = request.POST.get('organism', None)
        input_query = []
        for element in request.POST['input'].split('\n'):
            element = element.strip()
            if element:
                input_query.append(element)

        input_query = list(set(input_query))

        # max input IDs
        if 2000 > len(input_query) > 1:
            if re.match(r'^ENS\w*[GTP]', input_query[0]):
                job_num = str(random.randrange(500))
                with open(f'{jobs_path}/{job_num}.txt', "wb") as fp:
                    pickle.dump(input_query, fp)
                return redirect(Multi_proteins, organism=organism, job=job_num)
        else:
            error_message = f"The input list must contain between 2 and 2000 Ensembl IDs."
            jump_div = 'option1'

    # Option 2: Upload file
    if "option2" in request.POST and 'gene-count-file' in request.FILES:
        error_message_suffix = ""
        organism = request.POST.get('organism', None)
        try:
            # --- Check input file for correct format
            # Try to decode as UTF-8, sanitize and parse as table
            try:
                file_string = escape(request.FILES['gene-count-file'].read().decode('UTF-8'))
                file_buffer = StringIO(file_string)
                # Parse as pandas dataframe
                transcript_count_df = pd.read_table(file_buffer)
            except UnicodeDecodeError:
                error_message_suffix = "could not be parsed as an text file"
                raise RuntimeError

            except ParserError:
                error_message_suffix = f"could not be parsed as an table file (CSV or TSV)"
                raise RuntimeError

            # Check input shape
            if transcript_count_df.shape[0] < 2 or transcript_count_df.shape[1] < 2:
                error_message_suffix = f"could not be parsed as table or has less than two rows and columns"
                raise RuntimeError

            # Kevin: Zakaria please insert the magic down below:)
            # Zaka: And this is where the magic happens :p
            # Elias: lets hope this magic still works after I changed the format checking

            # Check if the first row corresponds to transcript Ensembl IDs
            if not (re.match(r'^ENS\w*T', transcript_count_df.iloc[0, 0]) or re.match(r'^ENS\w*T',
                                                                                      transcript_count_df.iloc[1, 0])):
                error_message_suffix = f"must have Ensembl transcript IDs in the first column starting with \"ENST\""
                raise RuntimeError

            # --- Try parsing counts for the different options (search for FPKM, tpm or counts)
            # max_isoforms: the max number of isoforms to consider:
            max_isoforms = int(request.POST['transcript-count-max'])

            column_names = transcript_count_df.columns

            # Cufflinks file (or a similar thing)
            if "FPKM" in column_names:
                transcript_count_df = transcript_count_df.sort_values(by=['FPKM'], ascending=False)
                cut_rows = transcript_count_df.iloc[:, 0].unique()[:max_isoforms]
                print('Input matches cufflinks output ')

            # Kallisto output counts in tpm
            elif "tpm" in column_names:
                transcript_count_df = transcript_count_df.sort_values(by=['tpm'], ascending=False)
                cut_rows = transcript_count_df.iloc[:, 0].unique()[:max_isoforms]
                print('Input matches kallisto output ')

            # Generic count matrix
            elif "counts" in column_names:
                transcript_count_df = transcript_count_df.sort_values(by=['counts'], ascending=False)
                cut_rows = transcript_count_df.iloc[:, 0].unique()[:max_isoforms]
                print('Input with counts column ')

            # Could not find the row
            else:
                error_message_suffix = "does not contain a column with the counts. The column must be named either \"FPKM\", \"tpm\" or \"counts\""
                raise RuntimeError

            # and let DIGGER do the magic ;)
            job_num = str(random.randrange(500))
            with open(f'{jobs_path}/{job_num}.txt', "wb") as fp:
                pickle.dump(cut_rows, fp)  # Pickling
                print(f"Starting network analysis with {len(cut_rows)} rows")
            return redirect(Multi_proteins, organism=organism, job=job_num)

        except RuntimeError:
            print("Could not parse uploaded file correctly")
            error_message = f"The uploaded file \"{request.FILES['gene-count-file']}\" {error_message_suffix}."
            jump_div = 'option2'

    return render(request, 'setup/network.html', context={
        'error_message': error_message,
        'jump_div': jump_div
    })


def Multi_proteins(request, organism, job='0'):
    with open(f'{jobs_path}/{job}.txt', "rb") as fp:  # Unpickling
        inputs = pickle.load(fp)

    if re.match(r'^ENS\w*G', inputs[0]):
        info = nt.analysis_input_genes(inputs, organism)

    elif re.match(r'^ENS\w*[TP]', inputs[0]):
        info = nt.analysis_input_isoforms(inputs, organism)
    else:
        return HttpResponse("<h1>wrong entry</h1>")

    if not info:
        return HttpResponse("<h1>Too many inputs (max=2000 genes)</h1>")

    else:
        genes, missing, num_isoforms = info

        Net = nt.Construct_network(genes, missing, job, organism)

        if Net == 0:
            return HttpResponse("<h1>There is no known interaction between these proteins</h1>")

        else:
            nodes, edges, tab, tb_html = Net

    context = {

        'pv_nodes': nodes,
        'pv_edges': edges,
        "tab": tb_html,
        'ID': job,
        "genes_number": len(missing),
        "isoforms_num": num_isoforms,
        'interacted_nodes': len(nodes),

    }

    return render(request, 'visualization/network.html', context)


def get_organisms(request):
    directory_path = os.path.join(settings.PROJECT_ROOT, 'domain/data')

    if os.path.exists(directory_path) and os.path.isdir(directory_path):
        organisms = os.listdir(directory_path)
        true_organisms = []
        trivial_names = []
        for orga in organisms:
            if not os.path.isdir(os.path.join(directory_path, orga)):
                continue
            true_organisms.append(orga.split("[")[0])
            trivial_names.append(orga.split("[")[1][:-1])
        return JsonResponse({'organisms': true_organisms, 'trivial_names': trivial_names})
    else:
        return JsonResponse({'organisms': [], 'trivial_names': []})
