import requests
import os
from time import sleep
from datetime import datetime

import pandas as pd

from src.helpers.folder_magic import search_for_files
from src.helpers.table_magic import load_include_exclude_txt
from src.AmiGo2.get_Amigo2_API import build_full_url_from_go_id, get_col, download_from_amigo2

def get_single_df_from_path(table_path, columns):
    if not os.path.isfile(table_path):
        return None
    try:
        df = pd.read_csv(table_path, sep="\t", dtype=str)
        for col in columns:
            if col in df.columns:
                continue
            elif col == "term":
                term = os.path.basename(os.path.dirname(table_path))
                df[col] = term
            else:
                df[col] = ""
        new_df = df[columns]
        return new_df
    except Exception as e:
        print("\n\n coulnd't get df from path", table_path)
        print(str(e))
        return None


def get_all_genes_from_path(amigo_data):
    columns = get_col() + ["term", "group"]
    path_datatsv = search_for_files(amigo_data, "data", "tsv")
    df_list = []
    for data_tsv in path_datatsv:
        df = get_single_df_from_path(data_tsv, columns)
        if df is not None:
            df_list.append(df[columns])

    if df_list:
        return pd.concat(df_list)
    else:
        return None


def get_overview(path_to_amigo_overview):

    #no idea what i did here ??
    if type(path_to_amigo_overview) == list:
        solution = []
        for p in path_to_amigo_overview:
            solution.append(get_overview(p))
        return pd.concat(solution)

    try: #in case i ever fix this :D
        df = pd.read_csv(path_to_amigo_overview)
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



def get_term_name(df, dir_path):
    dir_name, term_name = os.path.split(dir_path)
    if term_name != "empty":
        return term_name, dir_path

    alternativ_term_name = df["has_participant_closure_label"].unique().tolist()
    if len(alternativ_term_name) == 1 and alternativ_term_name[0] != "":
        term_name = alternativ_term_name[0].replace(' ', '_')
        return term_name, os.path.join(dir_name, term_name)

    elif len(alternativ_term_name) == 2: #here I hope one is empty
        term_name = (alternativ_term_name[0] + alternativ_term_name[1]).replace(' ', '_')
        return term_name, os.path.join(dir_name, term_name)

    else:
        complicated_term_name = df["annotation_extension_class_closure_label"].tolist()
        complicated_term_name.sort()
        count_dict = {}

        for s in complicated_term_name:
            labels = s.split(',')
            for l in labels:
                l = l.strip()
                if l in count_dict:
                    count_dict[l] += 1
                else:
                    count_dict[l] = 1

        possible_terms = alternativ_term_name
        for k, v in sorted(count_dict.items(), key=lambda item: item[1]):
            if v == len(complicated_term_name):
                possible_terms.append(k)

        for terf in possible_terms:
            for s in complicated_term_name:
                labels = s.split(',')
                if terf == labels[-1] or terf == labels[-2]:
                    term_name = terf.replace(' ', '_')
                    return term_name, os.path.join(dir_name, term_name)

    return "NO_TERM", dir_path


def get_single_table(dir_path, url, columns, force_download, in_ex_group):

    if not force_download:
        file_path = os.path.join(dir_path, "data.tsv")
        df = get_single_df_from_path(file_path, columns + ["term", "group"])
        if df is not None:
            return df

    df = download_from_amigo2(url, columns, dir_path)
    if df is None:
        return None

    term_name, dir_path = get_term_name(df, dir_path)
    if (
        term_name in in_ex_group["exclude"] or 
        term_name in in_ex_group["plusplus"] #or 
        #term_name == "NO_TERM"
    ):
        return None

    df["term"] = term_name
    df["group"] = in_ex_group["group"][term_name] if term_name in in_ex_group["group"].keys() else "NO_GROUP"
    save_table(df, dir_path, term_name)
    return df


def save_table(df, dir_path, term):
    file_path = os.path.join(dir_path, "data.tsv")
    os.makedirs(dir_path, exist_ok=True)
    if term == "NO_TERM":
        i=0
        while os.path.isfile(file_path):
            file_path = os.path.join(dir_path, f"data_{i}.tsv")
            i+=1
    df.to_csv(file_path, index=False, sep="\t")


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

def download_data(amigo_send, amigo_config, go_ids=None):
    if go_ids is None:
        go_ids = []
    data_path, overview_file, include_exclude_file, force_download = amigo_config
    overview_df = get_overview(overview_file)
    offline_data = get_all_genes_from_path(data_path)
    dir_names = overview_df["Name"].tolist()
    go_ids += overview_df["Accession"].tolist()

    columns = get_col()
    in_ex_group = load_include_exclude_txt(include_exclude_file)
    count_dict = {}
    count_dict_2 = {}
    x = 120
    for url in build_full_url_from_go_id(go_ids):
        dir_name = "NO_TERM"
        if dir_names:
            dir_name = dir_names.pop()
            if dir_name in in_ex_group["exclude"] or dir_name in in_ex_group["plusplus"]:
                continue
        dir_path = os.path.join(data_path, "genes", dir_name)
        if offline_data is not None:
            sub_df = offline_data[offline_data["term"] == dir_name]
            for name, symbol in check_count(sub_df, count_dict_2):
                amigo_send.send((name, symbol))
        df = get_single_table(dir_path, url, columns, force_download, in_ex_group)

        if df is None:
            df = get_single_table(dir_path, url, columns, force_download, in_ex_group)

        print(f"{datetime.now().strftime('%H%M')} Amigo2 got {os.path.basename(dir_path)}")
        for name, symbol in check_count(df, count_dict):
            amigo_send.send((name, symbol))

    amigo_send.send(("finished", "finished"))
    amigo_send.close()
    print("Amigo thread done")