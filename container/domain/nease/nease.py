import ast
import os.path

import pandas as pd

from .functions import *
import statsmodels.api as sm
import gseapy as gp
import pickle5 as pickle
import traceback


# Main Code 

class run(object):

    def __init__(self,
                 input_data,
                 organism='human',
                 input_type='Standard',
                 p_value_cutoff=0.05,
                 remove_non_in_frame=True,
                 only_divisible_by_3=False,
                 min_delta=0.05,
                 Majiq_confidence=0.95,
                 only_DDIs=False,
                 confidences=None,
                 ):

        """
            Create an instance of NEASE object.
            
            
            Parameters
            ---------- 
            data: dataframe 
                Splicing (or diff. splicing) exons  or events. 
               
                Standard input: Ensemble gene ID    Start of exon  End of exon     dPsi (optional)
                for MAJOQ output, Please change the  input_type to "MAJIQ"
                
            organism : str, optional
                NEASE 1.0 support only Human.
            
            input_type: str, optional
                
                Either "Standard","Spycone",'Whippet','rmats','DEXSeq'or "MAJIQ", If you need support of more types of outputs. Please contact: louadi@wzw.tum.de
                
            p_value_cutoff: float, optional
               The p value cutoff used to compute NEASE scores. (default is 0.05)
               
            remove_non_in_frame: boolean, default: True
                   If True remove all exons that are predicted to disturb the ORF or known to result in a non-coding gene. If you would like to change this option and include all exons, please change the parameter remove_non_in_frame to False.
           
           only_divisible_by_3: boolean, default: False
                    If True remove exons not divisible by 3
             
            min_delta: float, optional
                min delta to consider in case your input has dPsi column. The value also corresponds  to logfold change in case of Dexeq.   (default is 0.05)
            
            Majiq_confidence: float, optional
                In case of input_type='MAJIQ'. the parameter P(dPSI > 20%) is needed. Check MAJIQ paper for details about this  (default is 0.95)  
                
            only_DDIs: boolean, optional, default: False
               if True only use DDI and exclude PDB and linear motif interactions

            confidences: list, optional, default: None
                List of confidence levels to consider in the Join graph. If None, only original interactions are considered.
            
                
         """
        self.data = []
        self.organism = organism
        self.confidences = confidences

        self.summary = {}

        supported_organisms = ['human', 'mouse']
        supported_input_types = ['MAJIQ', 'Standard', 'Spycone', 'Whippet', 'rmats', 'DEXSeq']

        if organism not in supported_organisms:
            print('Error: Please choose one of the  supported  organism: "Human" and "Mouse".')

        if input_type not in supported_input_types:
            raise ValueError('Input type not supported')

        else:

            # Open the Join graph and databases of the selected organism:
            self.mapping = database_mapping[organism].copy()
            self.path = Pathways[organism].copy()
            self.ppi = PPI[organism].copy()
            self.only_DDIs = only_DDIs
            self.input_type = input_type

            self.elm = elm[organism].copy()
            self.elm_interactions = elm_interactions[organism].copy()
            self.pdb = pdb[organism].copy()

            self.non_coding = non_coding
            self.data = []
            self.spliced_genes = []

            self.cutoff = p_value_cutoff

            # I don't want to break anything downstream so I save it in a tmp variable
            tmp = self.mapping['NCBI gene ID'].astype(str)
            self.entrez_name_map = dict(zip(tmp, self.mapping['Gene name']))

            # preprocess data to make analysis easier
            try:
                self.path['entrez_gene_ids'] = self.path['entrez_gene_ids'].apply(eval)
            except TypeError as e:
                print(e)
            if not self.confidences:
                self.confidences = ['original']
            else:
                self.confidences.append('original')
            curr_net = network[organism].copy()
            curr_net.remove_edges_from([x for x in curr_net.edges(data=True)
                                        if x[2]['confidence'] not in self.confidences])
            self.Join = curr_net

            if input_type == 'MAJIQ':
                print("Using MAJIQ output")
                # Processing Majiq output
                self.data, self.spliced_genes, self.elm_affected, self.pdb_affected = process_MAJIQ(input_data,
                                                                                                    self.mapping,
                                                                                                    Majiq_confidence,
                                                                                                    min_delta,
                                                                                                    self.only_DDIs,
                                                                                                    self)
                if len(self.data) == 0:
                    self.summary['error'] = 'Found no overlap with protein domains.'

            elif input_type == 'Standard':
                print("using Standard output")
                try:
                    self.data, self.spliced_genes, self.elm_affected, self.pdb_affected, self.symetric_genes = process_standard(
                        input_data, self.mapping, min_delta, self.only_DDIs, self, remove_non_in_frame,
                        only_divisible_by_3)
                    if len(self.data) == 0:
                        self.summary['error'] = ('Found no overlap with protein domains. Make sure that the genomic '
                                                 'coordinates of the exons correspond to the human genome build hg38 (GRCh38).')
                except Exception as e:
                    raise ValueError('Could not recognize the standard format. Did you select the correct input '
                                     'format? <br> Please make sure your table matches '
                                     'the standard format with these columns: <code>Gene ensembl ID</code>, '
                                     '<code>EXON START</code>, <code>EXON END</code>, <code>dPSI (optional)</code>. '
                                     '<br> Also ensure that the genomic coordinates of the exons correspond '
                                     'to the correct genome build (i.e. hg38 (GRCh38) in human).')

            elif input_type == 'Whippet':
                try:

                    print("Using Whippet output")
                    input_data['Probability'] = input_data['Probability'].apply(float)
                    data = input_data[input_data['Probability'] >= 0.9]
                    data = data[abs(data['DeltaPsi']) >= min_delta]
                    # data = input_data
                    if len(data) == 0:
                        raise ValueError('No significant events found in the input data.')
                    # data = data.rename_axis('Gene ID').reset_index()
                    data = data.rename(columns={'Gene': 'Gene ID', 'DeltaPsi': 'dPSI'})
                    data['tmp'] = data['Coord'].apply(lambda x: x.split(':')[1])

                    data['start'] = data['tmp'].apply(lambda x: x.split('-')[0])
                    data['end'] = data['tmp'].apply(lambda x: x.split('-')[1])
                    data = data[['Gene ID', 'start', 'end', 'dPSI']]

                except:
                    raise ValueError('Invalid file format.')

                self.data, self.spliced_genes, self.elm_affected, self.pdb_affected, self.symetric_genes = process_standard(
                    data, self.mapping, min_delta, self.only_DDIs, self, remove_non_in_frame, only_divisible_by_3)

            elif input_type == 'rmats':

                try:
                    print("Using rMATS output")
                    data = input_data[input_data['FDR'] <= p_value_cutoff]
                    data = data[['GeneID', 'exonStart_0base', 'exonEnd', 'IncLevelDifference']]
                    data['GeneID'] = data['GeneID'].apply(lambda x: x.split('.')[0])

                except:
                    raise ValueError('Invalid file format.')

                self.data, self.spliced_genes, self.elm_affected, self.pdb_affected, self.symetric_genes = process_standard(
                    data, self.mapping, min_delta, self.only_DDIs, self, remove_non_in_frame, only_divisible_by_3)

            elif input_type == 'DEXSeq':

                try:
                    data = input_data[input_data['padj'] <= p_value_cutoff]
                    data = data[[data.columns[1], 'genomicData.start', 'genomicData.end', 'log2fold_control_case']]
                    print('proceding with log2fold threshold: ' + str(min_delta))
                except:
                    raise ValueError('Invalid file format.')

                self.data, self.spliced_genes, self.elm_affected, self.pdb_affected, self.symetric_genes = process_standard(
                    data, self.mapping, min_delta, self.only_DDIs, self, remove_non_in_frame, only_divisible_by_3)

            elif input_type == 'Spycone':

                try:

                    self.data, self.spliced_genes = process_spycone(input_data, self.mapping)
                    # spycone only uses DDI for now
                    self.only_DDIs = True

                except:
                    raise ValueError('Invalid file format.')

            if len(self.data) == 0:  #
                self.summary['error'] = 'Found no overlap with protein domains. Analysis cancelled...'
                raise Exception("No overlap with protein domains. Did you select the right organism?")

            else:
                self.data = self.data.drop_duplicates(['Gene name', 'NCBI gene ID', 'Gene stable ID', 'Pfam ID'],
                                                      keep='first')

                # check interaction of the domains
                # prepare elm to doamin interactions file
                self.elm_interactions['interactor 1'] = self.elm_interactions['Interator gene 1'].astype('str') + "/" + \
                                                        self.elm_interactions['Elm id of gene 1']
                self.elm_interactions['interactor 2'] = self.elm_interactions['Interator gene 2'].astype('str') + "/" + \
                                                        self.elm_interactions['Domain of gene 2']

                self.data = exons_to_edges(self.data, self.Join, self.elm_interactions, self.organism)

                # Identify binding of affected domains = Edges in the PPI
                # get DMI adn DDI for spliced domains

                self.interacting_domains = affected_edges(self, self.Join, self.only_DDIs)

                #get all edges of a gene from DDIs
                # self.interacting_domains,: DMI and DDI
                # self.pdb_affected: co-resolved interactions
                self.g2edges = gene_to_edges(self.interacting_domains, self.pdb_affected, self.only_DDIs)

                self.summary['domain_affected'] = str(len(self.data['Domain ID'].unique()))
                if not self.only_DDIs:
                    self.summary['lin_motifs'] = str(len(self.elm_affected['ELMIdentifier'].unique()))
                    self.summary['residues'] = str(len(self.pdb_affected))
                self.summary['known_interactions'] = str(
                    len(self.data[self.data['Interacting domain']]['Domain ID'].unique()))

                self.summary['interaction_affected'] = str(
                    len([item for sublist in self.g2edges.values() for item in sublist]))

                print('Running enrichment analysis...')

                self.supported_database = list(self.path['source'].unique())
                self.enrichment = pathway_enrichment(self.g2edges, self.path, self.mapping, self.Join, organism,
                                                     p_value_cutoff, self.only_DDIs).reset_index(drop=True)
                print('NEASE enrichment done.')

    def save(self, file_path):
        folders = "/".join(file_path.split("/")[:-1])
        if not os.path.isdir(folders):
            os.makedirs(folders)
        # filter out the static files that can be read in during loading
        self.mapping = None
        self.path = None
        self.ppi = None
        self.elm = None
        self.elm_interactions = None
        self.pdb = None
        with open(file_path + ".pkl", 'wb') as f:
            pickle.dump(self, f)

    def get_stats(self, file_path=''):

        """
            Display the stats of affected domains by splicing.
            
            Parameters
            ---------- 
            file_path : str
                Path for saving the statistics figure.
            
        """

        if len(self.data) == 0:
            print('Processing failed')

        else:

            # initial number of genes with known affected features

            affecting_number = 0
            self.spliced_genes = list(set(self.spliced_genes))

            # number of genes with affected domains/number of all events (genes)
            genes_with_domain = self.data['Gene stable ID'].unique()
            domain_number = len(genes_with_domain)

            affecting_number = domain_number
            number_of_features = domain_number

            elm_number = 0
            pdb_number = 0

            if not self.only_DDIs:
                # genes with affected elm
                genes_with_elm = list(self.elm_affected['Gene stable ID'].unique())

                # number of genes with elm affected
                elm_number = len(genes_with_elm)

                # number of genes with affected pdb
                genes_with_pdb = list(self.pdb_affected['Gene stable ID'].unique())
                pdb_number = len(genes_with_pdb)

                affecting_number = len(list(set(list(genes_with_elm) + list(genes_with_pdb) + list(genes_with_domain))))

                number_of_features = number_of_features + elm_number + pdb_number

            affecting_percentage = round(affecting_number / len(self.spliced_genes), 2)

            stats_domains(affecting_percentage, number_of_features, domain_number, elm_number, pdb_number, file_path)
        return

    def get_domains(self):

        """
            Display the list of AS events affectting domains in NEASE format.
            
        Returns
        -------
        pd.DataFrame Object
        
        """

        if len(self.data) == 0:
            # return empty dataframe
            print('Processing failed')
            return pd.DataFrame()

        # elif self.organism=="Mouse":
        #         #no visualization available for mouse in DIGGER
        #         return self.data.drop(columns=['Domain ID','Visualization link'])
        # else:

        # DIGGER visualization available for Human and Mouse
        return self.data.drop(columns=['Domain ID', 'DDI', 'elm']).reset_index(drop=True)

    def get_elm(self):

        """
            Display list of linear motifs affected by AS.
            
        Returns
        -------
        pd.DataFrame Object
        
        """

        if self.only_DDIs:
            print('You ran NEASE with the option:  only_DDIs=True. No ELM data available.')
            return pd.DataFrame()

        elif self.elm_affected.empty:
            print('No known linear motif are affected by AS')

        else:
            self.elm_affected['ELM link'] = self.elm_affected.apply(create_elm_link, axis=1)

        return self.elm_affected.drop(columns=['ID', 'Affected binding (NCBI)'], errors='ignore').reset_index(drop=True)

    def get_pdb(self):

        """
            Display list of linear motifs affected by AS.
            
        Returns
        -------
        pd.DataFrame Object
        
        """
        if self.only_DDIs:
            print('You ran NEASE with the option: only_DDIs=True. No pdb data available.')

        elif self.pdb_affected.empty:
            print('No residue from the PDB database motif are affected by AS')

        else:
            pdb_affected = self.pdb_affected.rename(
                columns={"symbol": "Gene name", 'entrezgene': 'Co-resolved interactions'}).copy()
            # Convert IDs to names
            c = lambda x: [Entrez_to_name(gene, mapping_dict=self.entrez_name_map) for gene in list(set(x))]
            pdb_affected['Co-resolved interactions symbol'] = pdb_affected['Co-resolved interactions'].apply(c)

            a = lambda x: ", ".join([str(val) for val in x])
            pdb_affected['Co-resolved interactions'] = pdb_affected['Co-resolved interactions'].apply(a)
            pdb_affected['Co-resolved interactions symbol'] = pdb_affected['Co-resolved interactions symbol'].apply(a)

            return pdb_affected[['Gene name', 'NCBI gene ID', 'Gene stable ID', 'Co-resolved interactions symbol',
                                 'Co-resolved interactions']].reset_index(drop=True)
        return self.pdb_affected.rename(
            columns={"symbol": "Gene name", 'entrezgene': 'Co-resolved interactions'}).copy()

    def get_edges(self):

        """   
            Display affected interactions from AS. 
            
        Returns
        -------
        pd.DataFrame Object
        
        """
        if len(self.data) == 0:
            print('Processing failed')
            return pd.DataFrame()
        elif len(self.interacting_domains) == 0:
            print('No affected edges identified.')
        else:

            edges = self.interacting_domains.copy()
            # Convert IDs to names
            c = lambda x: [Entrez_to_name(gene, mapping_dict=self.entrez_name_map) for gene in list(set(x))]
            edges['Affected binding'] = edges['Affected binding (NCBI)'].apply(c)

            # count number of affected PPI for every domain
            count = lambda x: len(x)
            edges['Number of affected interactions'] = edges['Affected binding'].apply(count)

            a = lambda x: ", ".join(x)
            edges['Affected binding'] = edges['Affected binding'].apply(a)
            edges['Affected binding (NCBI)'] = edges['Affected binding (NCBI)'].apply(a)
            edges = edges.drop_duplicates()
            edges = edges.sort_values('Number of affected interactions', ascending=False)

            return edges[['Gene name', 'NCBI gene ID', 'Identifier', 'dPSI', 'Number of affected interactions',
                          'Affected binding', 'Affected binding (NCBI)']].reset_index(drop=True)
        return self.interacting_domains

    def get_databases(self):
        return {'classic_dbs': gp.get_library_name(organism=self.organism), 'nease_dbs': self.supported_database}

    def classic_enrich(self, gseapy_databases, non_symmetrical_only=False, cutoff=0.05):

        """
        Classic gene level enrichement using the python library gseapy.

        Parameters
        ----------
        gseapy_databases: str, list, tuple of Enrichr Library name(s).
                  or custom defined gene_sets (dict, or gmt file). For more details, please check: https://pypi.org/project/gseapy/

        non_symmetrical_only: Bool
                    Run classic gene set enrichment only on non_symmetrical exons.

         Returns
        -------
        pd.DataFrame Object

        Example:
            gseapy_databases=['KEGG_2019_Human', 'Reactome_2016','WikiPathways_2019_Human']

        """

        if not self.spliced_genes:
            print('No genes found on your input.')
            return

        if not gseapy_databases:
            raise ValueError('Please provide a list of gene set databases to run the enrichment analysis.')

        # get gene sets supported in gseapy
        gseapy_library = gp.get_library_name(organism=self.organism)

        # compare with user input with gseapy_library
        gene_set_database = list(set(gseapy_library).intersection(gseapy_databases))
        # if nothing is found, use all the gene sets in gseapy library
        if len(gene_set_database) == 0:
            gene_set_database = list(set(gseapy_library))

        if not gene_set_database:
            print('none of the gene set databases provided is supported in gseapy library.')
            print('please check https://pypi.org/project/gseapy/')
            return

        # run gene set enrichment
        if non_symmetrical_only:
            if self.input_type == "Majiq":
                print(
                    'Non-symmetrical exons enrichment is not available for Majiq output. Please use standard input. ')

            else:
                # run only on non-symteric
                gene_list = [Ensemb_to_name(x, self.mapping) for x in self.symetric_genes]
        else:
            # run on all genes
            gene_list = [Ensemb_to_name(x, self.mapping) for x in self.spliced_genes]

        enr = gp.enrichr(gene_list=gene_list, organism=self.organism, gene_sets=gene_set_database,
                         outdir=None, cutoff=cutoff)

        # add a "significance" column
        enr.results['Significant'] = ['yes' if x <= cutoff else 'no' for x in
                                      enr.results['Adjusted P-value']]

        return enr.results.sort_values('Adjusted P-value')

    def enrich(self, database=None, cutoff=0.05):

        """ 
        Run enrichement analysis
        
        
        Parameters
        ---------- 
        database: List 
            gene set sources for enrichment.
        cutoff: float
            The p value cutoff.
        
        
            
        Returns
        -------
        pd.DataFrame Object
            
        Example: 
                events=nease.run(table, organism='Human')
                events.enrich(database=['Reactome'])
        
        """
        if not database:
            database = ['PharmGKB', 'HumanCyc', 'Wikipathways', 'Reactome', 'KEGG', 'SMPDB', 'Signalink', 'NetPath',
                        'EHMN', 'INOH', 'BioCarta', 'PID']

        if len(self.data) == 0:
            print('Processing failed')
            return None
        elif len(self.interacting_domains) == 0:
            print('No affected edges identified.')
            return None

        # Check if user input matches the available databases
        database = [x for x in database if x in self.supported_database]

        if len(database) == 0:
            print('Please select a supported pathway database as an argument. ')
            print('supported databases for ' + self.organism + ' are :', [x for x in self.supported_database], ".")
            print('\n')
            return None

        enrich_results = self.enrichment[self.enrichment['Source'].isin(database)]

        # Correct for multiple testing
        # fdr_bh : Benjamini/Hochberg (non-negative)
        enrich_results['adj p_value'] = \
            sm.stats.multipletests(list(enrich_results['p_value']), method='fdr_bh', alpha=0.05)[1]
        #enrich_results['adj p_value']=sm.stats.fdrcorrection(list(enrich_results['p_value']),alpha=0.01)[1]

        # shift column 'score to the last position
        scores = enrich_results.pop('Nease score')
        enrich_results.insert(len(enrich_results.columns), 'Nease score', scores)

        print('NEASE enrichment for the pathway databases:\n', [x for x in database])
        num = len(enrich_results[enrich_results['adj p_value'] <= cutoff])
        if num == 0:
            print('No enrichment found with the cutoff ' + str(cutoff) + '.')
        else:
            print("Found " + str(num) + " enriched pathways after multiple testing correction.\n")

        # add a "significance" column
        enrich_results['Significant'] = ['yes' if x <= cutoff else 'no' for x in
                                         enrich_results['adj p_value']]

        enrich_results = enrich_results.drop_duplicates()
        return enrich_results.sort_values(['Nease score', 'p_value'], ascending=[False, True]).reset_index(
            drop=True)

    def path_analysis(self, path_id):

        '''
        Run enrichment analysis on a specific pathway with details of impact of AS.
        
        Parameters
        ---------- 
        path_id: str 
            The ID of the pathway f interest.
            Run enrich() to find enriched pathways and their IDs.
            
            
        Returns
        -------
        pd.DataFrame Object
        
        exanple:
                # Run general enrichment
                events=nease.run(table, organism='Human')
                events.enrich(database=['Reactome']
                
                # Pathway specific analysis
                events.path_analysis('R-HSA-388396')

        
        '''

        if self.data.empty:
            raise ValueError('Processing failed')
        elif self.interacting_domains.empty:
            raise ValueError('No affected edges identified.')

        path_info = self.enrichment[self.enrichment['Pathway ID'] == path_id]

        if len(path_info) == 0:
            raise ValueError('No pathway with the given ID found.')

        else:
            # run enrichment
            enrich, _ = single_path_enrich(path_id, self.path, self.g2edges, self.mapping, self.organism,
                                           self.only_DDIs, self.entrez_name_map, self.Join)

            if len(enrich) == 0:
                raise ValueError('No enrichment or genes found for the selected pathway.')
            else:
                return enrich.sort_values(['p_value']).reset_index(drop=True)

    def Vis_path(self,
                 path_id,
                 k=0.8):

        '''
               Visualize the network module of a specific pathway.
               
               
            Parameters
            ---------- 
            path_id: str 
                The ID of the pathway f interest.
                Run enrich() to find enriched pathways and their IDs.
                
            file: str
                 A string representing a local file path for the html file.

            k: float 
                 -  Float (default=None))is a parameter to be tuned by the user:
                            Position nodes using Fruchterman-Reingold force-directed algorithm.
                            Optimal distance between nodes. If None the distance is set to 1/sqrt(n) where n is the number of nodes. 
                            Increase this value to move nodes farther apart.
                    Link: networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.spring_layout.html.
                
            auto_open: Boolean 
            
            
            Returns
            -------
            pd.DataFrame Object

            exanple:
                    # Run general enrichment
                    events=nease.run(table, organism='Human')
                    events.enrich(database=['Reactome']

                    # Pathway specific analysis
                    events.path_analysis('R-HSA-388396')

               
            path_id:    -  str: An unique pathway id. Run  enrich() to get all pathways and their ids.
            K:          -  Float (default=None))is a parameter to be tuned by the user:
                            Position nodes using Fruchterman-Reingold force-directed algorithm.
                            Optimal distance between nodes. If None the distance is set to 1/sqrt(n) where n is the number of nodes. 
                            Increase this value to move nodes farther apart.
                        Link: networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.spring_layout.html.
            file         - A string representing a local file path.
            '''
        if not k:
            k = 0.8
        k = float(k)

        if self.data.empty:
            raise ValueError('Processing failed')
        elif self.interacting_domains.empty:
            raise ValueError('No affected edges identified.')

        path_info = self.enrichment[self.enrichment['Pathway ID'] == path_id]

        if len(path_info) == 0:
            raise ValueError('No pathway with the given ID found.')

        path_name = list(path_info['Pathway name'])[0]
        print('Enrichment of the pathway: ' + path_name + '.\n')
        print('Overall p_value: ', list(path_info['p_value'])[0])
        print('\n')
        # run enrichment

        try:
            enrich, affected_graph = single_path_enrich(path_id, self.path, self.g2edges, self.mapping, self.organism,
                                                        self.only_DDIs, self.entrez_name_map, self.Join)
        except:
            traceback.print_exc()
            enrich = pd.DataFrame()

        if len(enrich) == 0:
            raise ValueError('No enrichment or genes found for the selected pathway.')

        # Get genes of the pathway (Entrez IDs)
        path_genes = list(self.path[self.path['external_id'] == path_id]['entrez_gene_ids'])[0]

        significant = list(enrich[enrich['p_value'] <= 0.05]['NCBI gene ID'].unique())

        if not isinstance(path_genes, list):
            path_genes = ast.literal_eval(path_genes)

        try:
            graph_data, missing_flag = extract_subnetwork(path_genes,
                                                          self.ppi,
                                                          list(enrich['NCBI gene ID'].unique()),
                                                          self.spliced_genes,
                                                          k,
                                                          self.mapping,
                                                          affected_graph,
                                                          significant,
                                                          self.entrez_name_map,
                                                          self.organism)
        except Exception as e:
            print(e)
            traceback.print_exc()
            raise ValueError('Something went wrong while creating the network.')

        print("extracted subnetwork")

        path_info = self.enrichment[self.enrichment['Pathway ID'] == path_id]
        path_name = list(path_info['Pathway name'])[0]

        fig_text = ""

        if missing_flag:
            # add a warning message to the figure text
            fig_text += "⚠️ Some genes could not be translated, network might be incomplete."

        fig = go.Figure(data=graph_data,
                        layout=go.Layout(
                            title=f"</b>{path_name}</b>",
                            titlefont_size=16,
                            showlegend=False,
                            hovermode='closest',
                            margin=dict(b=20, l=5, r=5, t=40),
                            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                            paper_bgcolor='white',
                            plot_bgcolor='white',
                            updatemenus=[dict(
                                buttons=list([
                                    dict(label="Show All Edges",
                                         method="update",
                                         args=[{"visible": [True, True, True, True]}]),
                                    dict(label="Hide Background Edges",
                                         method="update",
                                         args=[{"visible": [True, False, False, True]}])
                                ]),
                                direction="down",
                                showactive=True
                            )],
                            annotations=[dict(
                                text=fig_text,
                                showarrow=False,
                                xref="paper", yref="paper",
                                x=0.0,
                                y=0,
                                xanchor='left',
                                yanchor='bottom',
                                font=dict(size=14)
                            )]
                        ))

        # get the html as a string
        try:
            html = fig.to_html(full_html=True, include_plotlyjs='cdn',
                               config={'toImageButtonOptions':
                                           {'format': 'png', 'filename': f"vis_{path_name}", 'scale': 2}})
        except Exception as e:
            print(e)
            raise ValueError('Something went wrong while creating the network.')

        return html

    def get_p_value(self):
        return self.cutoff

    def vis_pathway_connection(self, enrichment, databases, k=1):
        supported_dbs = ['Reactome', 'KEGG']
        if not any([x in databases for x in supported_dbs]):
            return None

        enrichment_filtered = None
        pathway_graph = None
        network_pathway_names = []
        for db in supported_dbs:
            if db not in databases:
                continue
            if enrichment_filtered is None:
                enrichment_filtered = enrichment[enrichment['Source'] == db].copy()
                pathway_graph = pathway_hierarchy[self.organism][db]
                network_pathway_names.append(db)
            else:
                enrichment_filtered = pd.concat([enrichment_filtered, enrichment[enrichment['Source'] == db].copy()])
                pathway_graph = nx.compose(pathway_graph, pathway_hierarchy[self.organism][db])
                network_pathway_names.append(db)

        try:
            graph_data = all_pathway_network(enrichment_filtered, pathway_graph, k,
                                             db_name=','.join(network_pathway_names))
        except Exception as e:
            print(e)
            return None

        if graph_data is None:
            return None

        fig = go.Figure(data=graph_data,
                        layout=go.Layout(
                            title=f"</b>Connections of significant pathways for {', '.join(network_pathway_names)}</b>",
                            titlefont_size=16,
                            showlegend=False,
                            hovermode='closest',
                            margin=dict(b=20, l=5, r=5, t=40),
                            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                            paper_bgcolor='white',
                            plot_bgcolor='white',
                        ))
        html = fig.to_html(full_html=True, include_plotlyjs='cdn',
                           config={'toImageButtonOptions':
                                       {'format': 'png', 'filename': 'pathway_connections', 'scale': 2}})
        return html


def load(file_path):
    nease_object = pickle.load(open(file_path, 'rb'))
    # load the static files
    nease_object.mapping = database_mapping[nease_object.organism].copy()
    nease_object.path = Pathways[nease_object.organism].copy()
    nease_object.ppi = PPI[nease_object.organism].copy()

    nease_object.elm = elm[nease_object.organism].copy()
    nease_object.elm_interactions = elm_interactions[nease_object.organism].copy()
    nease_object.pdb = pdb[nease_object.organism].copy()

    return nease_object
