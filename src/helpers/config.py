import os
import json
import pandas as pd
from datetime import datetime


from src.helpers.folder_magic import search_for_file, search_for_files#, save_table
from src.helpers.std_out import send_message
#######################################
# all functions in this document are 
# necessary before the start of Steps
# 1. config functions
# 2. GO_ID's & disgnet data
# 3. include_exclude_files
#######################################

"""
1. config functions
"""
def get_config(config_item=None, what=None):
    match what:
        case 1 | "Amigo" | "amigo" | "AmiGo" | "AmiGo2" | "amigo2" | "Amigo2":
            return get_config_amigo(config_item)
        case 2 | "disgnet":
            return get_config_disgnet(config_item)
        case 3 | "hgnc" | "HGNC" | "hugo":
            return get_config_hgnc(config_item)
        case 4 | "ensembl" | "ensemble":
            return get_config_ensembl(config_item)
        case 5 | "gnomad" | "Gnomad" | "gnomAD" | "GnomAD":
            return get_config_gnomad(config_item)
        case 6 | "vep":
            return get_config_vep(config_item)
        case 7 | "figure" | "plots":
            return get_config_results(config_item)
        case _:
            return get_config_from_path(config_item)

# load from path
def get_config_from_path(config_path = None):
    config = {}
    try:
        if config_path is not None and os.path.isfile(config_path):
            with open(config_path, 'r') as f:
                config = json.load(f)
        elif isinstance(config_path, dict):
            config = config_path
    except Exception as e:
        print("!!!NO CONFIG!!!")
        print("\n",str(e),"\n")
        print("...trying to build it from scratch!")
        config = {}

    if not config:
        config = {
            "flags": {
                "download_data": False,
                "threads": 10,
                "top_genes": 200,
                "extra": False,
            },
            "absolute_file_paths": {
                "config": "",
                "data": "",
                "results": "",
            },
            "populations": ["afr", "nfe"]
        }
        this_folder = os.getcwd()
        for item in ["config", "data", "results"]:
            abs_path = os.path.join(this_folder, item)
            os.makedirs(abs_path, exist_ok=True)
            config["absolute_file_paths"][item] = abs_path

        config_file = os.path.join(this_folder, "config", "config.json")
        with open(config_file, 'w') as file:
            json.dump(config, file)

    print("config build")
    return config


def get_config_amigo(config):
    storage = config.get("absolute_file_paths")
    data_path = os.path.join(storage.get("data"), "Amigo2", "genes")
    result_path = os.path.join(storage.get("results"), "Amigo2")
    download = config.get("flags", {}).get("download_data", False)
    extra = config.get("flags", {}).get("extra", False)

    os.makedirs(data_path, exist_ok=True)
    os.makedirs(result_path, exist_ok=True)

    config_path = storage.get("config")
    include_exclude_file = search_for_file(config_path, "AmiGo2_include_exclude", ".txt")
    group_in_ex = load_include_exclude_txt(include_exclude_file)

    # overview_file = search_for_file(config_path, "overview", ".tsv")
    # if not overview_file:
    overview_file = search_for_file(config_path, "overview", ".txt")
    overview_df = get_overview(config_path, overview_file)
    # print(overview_df)
    # send_message(overview_df)
    return (overview_df, group_in_ex, data_path, result_path, download, extra)


def get_config_disgnet(config):
    storage = config.get("absolute_file_paths")
    data_path = os.path.join(storage.get("data"), "disgnet", "genes")
    result_path = os.path.join(storage.get("results"), "disgnet")
    config_path = storage.get("config")
    download = config.get("flags", {}).get("download_data", False)
    extra = config.get("flags", {}).get("extra", False)

    os.makedirs(result_path, exist_ok=True)

    include_exclude_file = search_for_file(config_path, "disgnet_include_exclude", ".txt")
    group_in_ex = load_include_exclude_txt(include_exclude_file)

    return (group_in_ex, data_path, result_path, download, extra)


def get_config_hgnc(config):
    storage = config.get("absolute_file_paths")
    data_path = os.path.join(storage.get("data"), "HGNC", "genes")
    result_path = os.path.join(storage.get("results"), "HGNC")
    download = config.get("flags", {}).get("download_data", False)

    os.makedirs(data_path, exist_ok=True)
    os.makedirs(result_path, exist_ok=True)

    config_path = storage.get("config")
    symbol_checker = [
            search_for_file(config_path, "look_up_hugo", ".tsv"),
            search_for_file(config_path, "check_amigo", ".csv"), 
            search_for_file(config_path, "check_disgnet", ".csv"), 
        ]
    if not download:
        symbol_checker.extend(search_for_file(result_path, "look_up_hugo", ".tsv"),)

    disgnet_df = prepare_disgnet_hgnc(config)

    return (disgnet_df, symbol_checker, data_path, result_path, download)


def get_config_ensembl(config):
    top_genes = config.get("flags",{}).get("top_genes", 200)
    download = config.get("flags", {}).get("download_data", False)
    storage = config.get("absolute_file_paths")
    data_path = os.path.join(storage.get("data"), "ensembl")
    result_path = storage.get("results")

    os.makedirs(data_path, exist_ok=True)
    #os.makedirs(result_path, exist_ok=True)

    return (data_path, result_path, top_genes, download)


def get_config_gnomad(config):
    storage = config.get("absolute_file_paths")
    data_path = os.path.join(storage.get("data"), "gnomAD")
    os.makedirs(data_path, exist_ok=True)
    download = config.get("flags", {}).get("download_data", False)
    populations = config.get("populations")
    return (data_path, populations, download)


def get_config_vep(config):
    storage = config.get("absolute_file_paths")
    data_path = os.path.join(storage.get("data"), "ensembl")
    os.makedirs(data_path, exist_ok=True)
    download = config.get("flags", {}).get("download_data", False)
    return (data_path, download)

def get_config_results(result_path):
    ts = datetime.now().strftime("%Y%m%dT%H")
    plots = os.path.join(result_path, f"{ts}_plots")
    hgnc_complete = search_for_file(result_path, "hgnc_df_complete", ".tsv")
    look_up_hugo = search_for_file(result_path, "look_up_hugo", ".tsv")
    amigo_sub = search_for_file(result_path, "amigo_df_sub", ".tsv")
    amigo_complete = search_for_file(result_path, "amigo_df_complete", ".tsv")
    disgnet_complete = search_for_file(result_path, "disgnet_df_complete", ".tsv")
    disgnet_sub = search_for_file(result_path, "disgnet_df_sub", ".tsv")

    return (result_path, plots, hgnc_complete, look_up_hugo, amigo_sub, amigo_complete, disgnet_complete, disgnet_sub)
"""
2. GO ID's & disgnet
"""
def get_overview(config, path_to_amigo_overview):

    # try:
    #     if path_to_amigo_overview.ends_with("tsv"):
    #         df = pd.read_csv(path_to_amigo_overview, sep="\t", dtype=str)
    #         return df
    # except Exception as _:
    #     print("No overview table found, looking for txt")


    overviewtxt = { "Accession" : [], "Name" : [], "Ontology" : [], "Synonyms" : [], "Alternate IDs" : [], "Definition" : [], "not_found" : [] }
    with open(path_to_amigo_overview, 'r') as f:
        last_line = ""
        for line in f:
            if line.startswith('#'):
                continue
            elif last_line == "":
                last_line = line.strip()
                continue
            elif last_line in overviewtxt:
                overviewtxt[last_line].append(line.strip())
                last_line = ""
            else:
                print(last_line, line)
                overviewtxt["not_found"].append(line.strip())
                last_line = ""

    not_found = overviewtxt.pop("not_found")
    if not_found:
        print("These Strings are unknown for the Amigo overview:")
        print(not_found)

    norm_accession = [ go_id.replace("GO:", "").strip() for go_id in overviewtxt["Accession"] ]
    overviewtxt["Accession"] = norm_accession

    df = pd.DataFrame.from_dict(overviewtxt)
    #save_table(df, config, path_to_amigo_overview)
    return df

    # if isinstance(path_to_amigo_overview, list):
    #     solution = []
    #     for p in path_to_amigo_overview:
    #         solution.append(get_overview(config, p))
    #     df = pd.concat(solution, ignore_index=True)
    #     table_name = os.path.basename(os.path.commonprefix(path_to_amigo_overview)) + "_overview_concat"
    #     save_table(df, config, table_name)
    #     return df

    # if not os.path.isfile(path_to_amigo_overview):
    #     overviewtxt = { "Accession" : [path_to_amigo_overview.replace("GO:", "").strip()], 
    #                    "Name" : ["NO_NAME"], "Ontology" : ["NO Ontology"], "Synonyms" : ["NO Synonym"],  "Alternate IDs" : ["NO alt ID"], "Definition" : ["NO Definition"] }
    #     return pd.DataFrame.from_dict(overviewtxt)




def load_include_exclude_txt(path):
    in_out_txt = {
        "group" : {},
        "include" : [],
        "exclude" : [],
        "extra" : []
    }
    with open(path, 'r') as f:
        group = ""
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            elif line.startswith('<'): 
                line = line[1:-1].strip()
                group = line
                continue
            elif line.startswith('--'): 
                line = line[2:].strip()
                in_out_txt["exclude"].append(line)
            elif line.startswith('++'): 
                line = line[2:].strip()
                in_out_txt["extra"].append(line)
                in_out_txt["group"][line] = group
            else:
                in_out_txt["include"].append(line)
                in_out_txt["group"][line] = group

    return in_out_txt


def build_disgnet_tables(disgnet_config):
    df = pd.DataFrame()
    group_in_ex, data_path, result_path, download, extra = disgnet_config
    resulttsv = search_for_file(result_path, "disgnet_df", "tsv")
    if resulttsv and not download:
        df = pd.read_csv(resulttsv, sep="\t", dtype=str)
    else:
        resulttsv = search_for_files(data_path, "search_result", "tsv")

        df_list = []
        for tsv in resulttsv:
            try:
                df = pd.read_csv(tsv, sep="\t", dtype=str)
                if df.empty:
                    continue
                df_list.append(df)
            except Exception as e:
                print(str(e))
        if df_list:
            df = pd.concat(df_list, ignore_index=True)
    if df.empty:
        return df

    df["group"] = df["disease_name"].map(lambda x : group_in_ex["group"].get(x, "NO GROUP"))
    file_name = os.path.join(result_path, "disgnet_df_complete.tsv")
    df.to_csv(file_name, sep="\t", index=False)

    if extra:
        df = df[df["disease_name"].isin(group_in_ex["include"] + group_in_ex["extra"])]
    else:
        df = df[df["disease_name"].isin(group_in_ex["include"])]

    df = df.drop_duplicates(ignore_index=True)
    file_name = os.path.join(result_path, "disgnet_df_sub.tsv")
    df.to_csv(file_name, sep="\t", index=False)
    return df

def prepare_disgnet_hgnc(config):
    send_message("disgnet started", 0, "disgnet")
    disgnet_config = get_config_disgnet(config)
    disgnet_df = build_disgnet_tables(disgnet_config)
    send_message( disgnet_df.shape[0], 1 ,"disgnet")
    send_message(disgnet_df.shape[0], 2 ,"disgnet")
    #group_in_ex, data_path, result_path, download, extra = disgnet_config
    hgnc_df = disgnet_df[["geneEnsemblIDs", "gene_symbol"]]
    send_message("finished", 0, "disgnet")
    return hgnc_df.drop_duplicates(ignore_index=True)