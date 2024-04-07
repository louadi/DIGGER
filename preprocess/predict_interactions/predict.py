# first clean protein interactions
# then use ppidm to predict and finally use graph attributes to annotate the graph
import timeit

import clean_protein_interactions as cpi
import create_uniprot_pfam_map as cupm
import parse_yaml
import ppidm
import preprocess.extend_digger.graph_attributes as ga
import preprocess.extend_digger.ddi_network as dn
import os

needed_files = ['3did_flat', 'INTERACTION.txt', 'mart_export.txt', 'pdb_chain_pfam.tsv',
                'uniprot_sprot.dat', 'uniprot_trembl.dat']

if __name__ == '__main__':
    # check if os is linux and stop if it is not
    if not os.name == 'posix':
        print("This script only works with linux systems")
        exit(1)
    for file in needed_files:
        if not os.path.isfile(f"../sourcedata/{file}"):
            print(f"Missing file: {file}, please read the README.md in sourcedata for more information.")
            exit(1)

    # parse yaml file with information about sources
    tasks, organism, functions, additional_flags = parse_yaml.parse("../sourcedata/database_sources.yml")

    # purely debug
    if 'none' in functions:
        print(tasks)
        print(organism)
        print(functions)
        print(additional_flags)
        exit(0)

    start = timeit.default_timer()

    # clean protein interactions from all sources
    if 'all' in functions or 'clean' in functions:
        cpi.main(tasks)
    # create uniprot pfam map needed for ppidm
    if 'all' in functions or 'map' in functions:
        cupm.main()
    # predict domain-domain interactions
    if 'all' in functions or 'predict' in functions:
        ppidm.main(additional_flags['ignore_threshold'])
    # create nodes necessary for graph attributes
    if 'all' in functions or 'cumulate' in functions:
        dn.main()
    # extend DIGGER graph with confidence levels
    if 'all' in functions or 'extend' in functions:
        ga.main(organism)

    stop = timeit.default_timer()
    print(f"Took : {stop - start} seconds for {functions}")