# after running this bash command on both uniprot_sprot.dat and uniprot_trembl.dat:
# awk '/^AC   /{ accession = $2 } /^DR   Pfam;/{ print accession, $0 }' uniprot_sprot.dat > extracted_data.txt
# execute this script to get the needed pfam-seq-sp and pfam-seq-tr files
import os


def create_map(data_file):
    uniprot_map = dict()
    with open(data_file, 'r') as file:
        for line in file:
            line_sp = line.rstrip().split(sep=";")
            if line_sp[0] in uniprot_map:
                uniprot_map[line_sp[0]].append(line_sp[2].strip())
            else:
                uniprot_map[line_sp[0]] = list()
                uniprot_map[line_sp[0]].append(line_sp[2].strip())
    return uniprot_map


def add_alternatives(map, data_file):
    # add the alternative sequences to the map. first run:
    # grep -P 'AC   ' uniprot_sprot.dat | cut -c 6- | grep -E '.*;.*;.*' | sed 's/;//g' | sed 's/ /\t/g' > uniprot_alternatives.tsv
    # for both uniprot_sprot.dat and uniprot_trembl.dat
    alternatives = 0
    with open(data_file, 'r') as file:
        for line in file:
            line_sp = line.rstrip().split(sep="\t")
            if line_sp[0] in map:
                domains = map[line_sp[0]]
            else:
                continue
            for alt in line_sp[1:]:
                map[alt] = domains
                alternatives += 1
    print("Added", alternatives, "alternatives")
    return map


def write_map(map, file):
    with open(file, 'w') as f:
        for key, value in map.items():
            for v in value:
                f.write(v + '\t' + key + '\n')


def main():
    # check if uniprot_sprot.dat is in sourcedata
    if not os.path.exists('../sourcedata/extracted_data.txt'):
        print("extracting data from uniprot_sprot.dat")
        os.system("awk '/^AC   /{ accession = $2 } /^DR   Pfam;/{ print accession, $0 }' "
                  "../sourcedata/uniprot_sprot.dat > ../sourcedata/extracted_data.txt")

    if not os.path.exists('../sourcedata/extracted_data_trembl.txt'):
        print("extracting data from uniprot_trembl.dat")
        os.system("awk '/^AC   /{ accession = $2 } /^DR   Pfam;/{ print accession, $0 }' "
                  "../sourcedata/uniprot_trembl.dat > ../sourcedata/extracted_data_trembl.txt")

    if not os.path.exists('../resultdata'):
        os.makedirs('../resultdata')
    map_sprot = create_map('../sourcedata/extracted_data.txt')
    map_trembl = create_map('../sourcedata/extracted_data_trembl.txt')
    if not os.path.exists('../sourcedata/uniprot_alternatives.tsv'):
        os.system("grep -P 'AC   ' ../sourcedata/uniprot_sprot.dat | cut -c 6- | grep -E '.*;.*;.*' | sed 's/;//g' | "
                  "sed 's/ /\t/g' > ../sourcedata/uniprot_alternatives.tsv")
    map_sprot = add_alternatives(map_sprot, '../sourcedata/uniprot_alternatives.tsv')

    if not os.path.exists('../sourcedata/uniprot_alternatives_tr.tsv'):
        os.system("grep -P 'AC   ' ../sourcedata/uniprot_trembl.dat | cut -c 6- | grep -E '.*;.*;.*' | sed 's/;//g' | "
                  "sed 's/ /\t/g' > ../sourcedata/uniprot_alternatives_tr.tsv")
    map_trembl = add_alternatives(map_trembl, '../sourcedata/uniprot_alternatives_tr.tsv')
    write_map(map_sprot, '../sourcedata/pfam-seq-sp')
    write_map(map_trembl, '../sourcedata/pfam-seq-tr')


if __name__ == '__main__':
    main()
