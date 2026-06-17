import requests
import os
from time import sleep
from datetime import datetime

from io import StringIO

import pandas as pd

from src.helpers.folder_magic import search_for_file
from src.helpers.table_magic import load_include_exclude_txt
from src.AmiGo2.db import get_tables_from_path, build_full_url_from_go_id, get_col, get_single_df_from_path

SEED = 42

def get_config(config):
    storage = config.get("absolute_file_paths")
    data_path = os.path.join(storage.get("data"), "Amigo2")
    config_path = storage.get("config")
    os.makedirs(data_path, exist_ok=True)

    overview = config.get("relative_file_paths").get("AmiGo2_overview")
    overview_file = os.path.join(config_path, overview)

    include_exclude = config.get("relative_file_paths").get("Amigo2_inclue_exclude")
    include_exclude_file = os.path.join(config_path, include_exclude)

    download = config.get("flags").get("download_data")

    return (data_path, overview_file, include_exclude_file, download)


def get_overview(path_to_amigo_overview):
    if type(path_to_amigo_overview) == list:
        solution = []
        for p in path_to_amigo_overview:
            solution.append(get_overview(p))
        return pd.concat(solution)

    try: #in case i ever fix this :D
        df = pd.DataFrame.from_csv(path_to_amigo_overview)
        return df

    except Exception as _:

        overviewtxt = {
            "Accession" : [],
            "Name" : [],
            "Ontology" : [],
            "Synonyms" : [], 
            "Alternate IDs" : [],
            "Definition" : [],
            "not_found" : []
        }

        with open(path_to_amigo_overview, 'r') as f:
            last_line = ""
            for line in f:
                #print(line)
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

        print(overviewtxt.pop("not_found"))

        norm_accession = [ 
            go_id.replace("GO:", "").strip()
            for go_id in overviewtxt["Accession"]
            ]
        overviewtxt["Accession"] = norm_accession

        return pd.DataFrame.from_dict(overviewtxt)


def get_tables(path_to_Amigo_db):
    return get_tables_from_path(path_to_Amigo_db)

def download_from_amigo2(url, columns, dir_path):
    #name = url[-40:].split("&fq=")[-1]  #looks a bit ugly, but should print the GO_id and then some
    #print(f"downloading GO ID {name}")
    r = requests.get(url, timeout=120)
    if r.status_code == 200:
        text = r.text
        try:
            global SEED
            SEED ^= SEED << 13
            SEED ^= SEED >> 7
            SEED ^= SEED << 17
            sleep(6 + SEED % 6)
            return pd.read_csv(StringIO(text), 
                            sep="\t", 
                            dtype=str, 
                            header=None,
                            names=columns)
        except Exception as e:
            print(str(e))
            print("couldn't read response")
            file_path = os.path.join(dir_path, "text.txt")
            os.makedirs(dir_path, exist_ok=True)
            with open(file_path, 'w') as f:
                f.write(text)
            print("response saved under: ", file_path)

    print("failed to download from:", url)
    return None


def get_term_name(df, dir_path):
    dir_name, term_name = os.path.split(dir_path)
    if term_name != "empty":
        return term_name, dir_path

    alternativ = df["has_participant_closure_label"].unique().tolist()
    if len(alternativ) == 1 and alternativ[0] != "":
        term_name = alternativ[0].replace(' ', '_')
        return term_name, os.path.join(dir_name, term_name)

    elif len(alternatives) == 2: #here I hope one is empty
        term_name = (alternativ[0] + alternativ[1]).replace(' ', '_')
        return term_name, os.path.join(dir_name, term_name)

    else:
        alternatives = df["annotation_extension_class_closure_label"].tolist()
        alternatives.sort()
        count_dict = {}

        for s in alternatives:
            labels = s.split(',')
            for l in labels:
                l = l.strip()
                if l in count_dict:
                    count_dict[l] += 1
                else:
                    count_dict[l] = 1

        possible_terms = alternativ
        for k, v in sorted(count_dict.items(), key=lambda item: item[1]):
            if v == len(alternatives):
                possible_terms.append(k)

        for terf in possible_terms:
            for s in alternatives:
                labels = s.split(',')
                if terf == labels[-1] or terf == labels[-2]:
                    term_name = terf.replac(' ', '_')
                    return term_name, os.path.join(dir_name, term_name)

    return "NO TERM", dir_path


def get_single_table(dir_path, url, columns, download, in_ex_group):
    if not download:
        file_path = os.path.join(dir_path, "data.tsv")
        df = get_single_df_from_path(file_path, columns + ["term"] + ["group"])
        if df is not None:
            return df

    df = download_from_amigo2(url, columns, dir_path)
    if df is None:
        return None

    term_name, dir_path = get_term_name(df, dir_path)
    if term_name in in_ex_group["exclude"] or term_name in in_ex_group["plusplus"]:
        return None

    os.makedirs(dir_path, exist_ok=True)
    df["term"] = term_name
    df["group"] = in_ex_group["group"][term_name] if term_name in in_ex_group["group"].keys() else "NO GROUP"


    file_path = os.path.join(dir_path, "data.tsv")
    if term_name == "NO TERM":
        i=0
        while os.path.isfile(file_path):
            file_path = os.path.join(dir_path, "data{i}.tsv")
            i+=1
    try:
        df.to_csv(file_path, index=False, sep="\t")
        return df
    except Exception as _:
        return None

def check_count(df, count_dict):
    if df is None:
        return []
    name_symbol = []
    for name, symbol in zip(df["bioentity_label"], df["bioentity_name"]):
        i = count_dict.get(name, 0)
        j = count_dict.get(symbol, 0)
        if i == 0:
            count_dict[name] = 0
        if j == 0:
            count_dict[symbol] = 0

        count_dict[name] += 1
        count_dict[symbol] += 1
        if i > 0 or j > 0:
            name_symbol.append((name, symbol))

    return name_symbol

def download_data(amigo_send, amigo_config, go_ids=[]):
    data_path, overview_file, include_exclude_file, download = amigo_config
    overview_df = get_overview(overview_file)
    offline_data = get_tables_from_path(data_path)
    go_ids += overview_df["Accession"].tolist()
    dir_names = overview_df["Name"].tolist()
    columns = get_col()
    in_ex_group = load_include_exclude_txt(include_exclude_file)
    count_dict = {}
    count_dict_2 = {}
    x = 120
    for url in build_full_url_from_go_id(go_ids):
        dir_name = "NO TERM"
        if dir_names:
            dir_name = dir_names.pop()
            if dir_name in in_ex_group["exclude"] or dir_name in in_ex_group["plusplus"]:
                continue
        dir_path = os.path.join(data_path, "genes", dir_name)
        if offline_data is not None:
            sub_df = offline_data[offline_data["term"] == dir_name]
            for name, symbol in check_count(sub_df, count_dict_2):
                amigo_send.send((name, symbol))
                sleep(1)
        df = get_single_table(dir_path, url, columns, download, in_ex_group)

        if df is None:
            df = get_single_table(dir_path, url, columns, download, in_ex_group)

        print(f"{datetime.now().strftime('%H%M')} Amigo2 got {os.path.basename(dir_path)}")
        for name, symbol in check_count(df, count_dict):
            amigo_send.send((name, symbol))
            sleep(1)


    amigo_send.send(("finished", "finished"))
    amigo_send.close()
    print("Amigo thread done")