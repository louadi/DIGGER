import yaml


def parse_mitab(parse_sources, counter=1):
    tasks = []
    for source in parse_sources['mitab']:
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

        tasks.append(("clean_mitab", path, f"sourcedata/source{counter}_{name}", relevant_columns))
        counter += 1
    return tasks, counter


def parse_mippie(parse_sources, counter=1):
    return parse_general(parse_sources, counter, 'mippie')


def parse_homology(parse_sources, counter=1):
    return parse_general(parse_sources, counter, 'homology')


def parse_string(parse_sources, counter=1):
    tasks = []
    for source in parse_sources['string']:
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

        tasks.append(("clean_string", path, f"sourcedata/source{counter}_string-{name}", mapping))
        counter += 1
    return tasks, counter


def parse_mint(parse_sources, counter=1):
    tasks = []
    try:
        name = parse_sources['mint']['name']
    except KeyError:
        name = 'mint'

    try:
        path = parse_sources['mint']['path']
    except KeyError:
        print(f"YAML not properly formatted, specify a path for mint")
        exit(1)

    tasks.append(("clean_mint", path, f"sourcedata/source{counter}_{name}"))
    counter += 1
    return tasks, counter


def parse_general(parse_sources, counter=1, col_name=None):
    tasks = []
    if col_name is None:
        return tasks, counter

    source = col_name
    try:
        name = parse_sources[col_name]['name']
    except KeyError:
        name = source

    try:
        path = parse_sources[col_name]['path']
    except KeyError:
        print(f"YAML not properly formatted, specify a path for mippie")
        exit(1)

    try:
        mapping = parse_sources[col_name]['mapping']
    except KeyError:
        print(f"YAML not properly formatted, specify mapping for {name}")
        exit(1)
    tasks.append((f"clean_{source}", path, f"sourcedata/source{counter}_{name}", mapping))
    counter += 1
    return tasks, counter


# feel free to extend this list by adding a new function and updating the dictionary
supported_sources = {'mitab': parse_mitab,
                     'string': parse_string,
                     'mippie': parse_mippie,
                     'homology': parse_homology,
                     'mint': parse_mint}


# parse yaml file
def parse(file):
    with open(file) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)

    sources_to_parse = {}
    counter = 1
    # extract mapping to decide what functions to call
    if 'sources' in data:
        sources_to_parse = data['sources']
    tasks = []
    additional_params = {}

    for supported_source in supported_sources:
        if supported_source not in sources_to_parse:
            continue
        source_tasks, counter = supported_sources[supported_source](sources_to_parse, counter)
        tasks.extend(source_tasks)

    if 'params' in data:
        for param in data['params']:
            additional_params[list(param.keys())[0]] = list(param.values())[0]

    return tasks, data['organism'], data['functions'], additional_params
