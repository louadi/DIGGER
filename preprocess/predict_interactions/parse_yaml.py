import yaml


# parse yaml file
def parse(file):
    with open(file) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)

    counter = 1
    # extract mapping to decide what functions to call
    sources_to_parse = data['sources']
    tasks = []
    additional_params = {}
    for source in sources_to_parse['mitab']:
        try:
            name = list(source.keys())[0]
        except AttributeError:
            print("YAML not properly formatted, specify a source")
            exit(1)

        try:
            path = list(source.values())[0]['path']
        except KeyError:
            print(f"YAML not properly formatted, specify a path for {list(source.keys())[0]}")
            exit(1)

        try:
            relevant_columns = list(source.values())[0]['interactor_columns']
        except KeyError:
            print(f"YAML not properly formatted, specify interactor columns for {list(source.keys())[0]}")
            exit(1)

        tasks.append(("clean_mitab", path, f"../sourcedata/source{counter}_{name}", relevant_columns))
        counter += 1

    for source in sources_to_parse['string']:
        try:
            name = list(source.keys())[0]
        except AttributeError:
            print("YAML not properly formatted, specify a source (i.e. physical, rest, all)")
            exit(1)

        try:
            path = list(source.values())[0]['path']
        except KeyError:
            print(f"YAML not properly formatted, specify a path for {list(source.keys())[0]}")
            exit(1)

        try:
            mapping = list(source.values())[0]['mapping']
        except KeyError:
            print(f"YAML not properly formatted, specify mapping for {list(source.keys())[0]}")
            exit(1)

        tasks.append(("clean_string", path, f"../sourcedata/source{counter}_string-{name}", mapping))
        counter += 1

    for source in ['mippie', 'homology']:
        if source in sources_to_parse:
            try:
                name = sources_to_parse[source]['name']
            except KeyError:
                name = source

            try:
                path = sources_to_parse[source]['path']
            except KeyError:
                print(f"YAML not properly formatted, specify a path for mippie")
                exit(1)

            try:
                mapping = sources_to_parse[source]['mapping']
            except KeyError:
                print(f"YAML not properly formatted, specify mapping for {name}")
                exit(1)
            tasks.append((f"clean_{source}", path, f"../sourcedata/source{counter}_{name}", mapping))
            counter += 1

    if 'mint' in sources_to_parse:
        try:
            name = sources_to_parse['mint']['name']
        except KeyError:
            name = 'mint'

        try:
            path = sources_to_parse['mint']['path']
        except KeyError:
            print(f"YAML not properly formatted, specify a path for mint")
            exit(1)

        tasks.append(("clean_mint", path, f"../sourcedata/source{counter}_{name}"))
        counter += 1

    if 'params' in data:
        for param in data['params']:
            additional_params[list(param.keys())[0]] = list(param.values())[0]

    return tasks, data['organism'], data['functions'], additional_params
