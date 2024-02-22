import pickle
import requests
import json

def get_gene_symbol(entrez_ids, species):
    url = "https://mygene.info/v3/gene"
    data = {"ids": ",".join([str(x) for x in entrez_ids]), "species": species}
    response = requests.post(url, params=data)
    result = json.loads(response.text)
    symbols = {}
    for item in result:
        try:
            symbols[item['_id']] = item.get('symbol', None)
        except KeyError:
            symbols[item['query']] = None
    return symbols


def unmapped_entrez_ids(g2n, ppis):
    c = 0
    unmapped = []
    for node in ppis.nodes:
        try:
            _ = g2n[node]
        except KeyError:
            c += 1
            unmapped.append(node)
    print(f"{c} nodes not found out of {len(ppis.nodes)}")
    return unmapped


def main():
    ppi_graph = pickle.load(open('../domain/data/Homo sapiens[human]/PPI.pkl', 'rb'))
    g2n = pickle.load(open('../domain/data/Homo sapiens[human]/gid2name.pkl', 'rb'))
    unmapped_ids = unmapped_entrez_ids(g2n, ppi_graph)
    # split unmapped ids in batches of max 500
    splits = [unmapped_ids[i:i + 500] for i in range(0, len(unmapped_ids), 1000)]
    all_symbols = {}
    missing = 0
    for ids in splits:
        all_symbols.update(get_gene_symbol(ids, ['human', 'mouse']))
    for key, value in all_symbols.items():
        if value is None:
            all_symbols[key] = key
            print(f"Replaced {key} since it is None")
            missing +=1
    print(f"{missing} ids are still not converted")
    g2n.update(all_symbols)
    pickle.dump(g2n, open('../domain/data/Homo sapiens[human]/gid2name_up.pkl', 'wb'))


if __name__ == '__main__':
    main()
