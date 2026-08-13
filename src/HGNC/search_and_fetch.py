import requests
import json 
import os
import pandas as pd
import numpy as np
from io import StringIO
from time import sleep
from datetime import datetime

from collections import Counter

from src.helpers.folder_magic import search_for_files, check_string
from src.helpers.std_out import send_message

# 10/s is the rate at which i can shoot requests 

# can use this with search and fetch
# i want hgnc_id and ens
searchableFields = [
    "alias_symbol",
    "symbol", # this is the symbol
    "uniprot_ids", # we have this in amigo
    "ccds_id",
    "locus_type",
    "mgd_id",
    "ena",
    "rna_central_id",
    "location",
    "prev_name",
    "ensembl_gene_id", #disgnet has this sometimes
    "hgnc_id",
    "status",
    "entrez_id",
    "alias_name",
    "vega_id",
    "mane_select",
    "ucsc_id",
    "refseq_accession",
    "omim_id",
    "rgd_id",
    "locus_group",
    "curator_notes",
    "name", #this we have as well
    "prev_symbol"
]

storedFields = [
    "status",
    "entrez_id",
    "ensembl_gene_id", #we want this
    "prev_name",
    "location",
    "mgd_id",
    "rna_central_id",
    "locus_type",
    "date_name_changed",
    "alias_symbol",
    "symbol",
    "gene_group",
    "prev_symbol",
    "gene_group_id",
    "gtrnadb",
    "name",
    "mamit-trnadb",
    "enzyme_id",
    "orphanet",
    "omim_id",
    "alias_name",
    "snornabase",
    "imgt",
    "date_modified",
    "hgnc_id",
    "merops",
    "mirbase",
    "horde_id",
    "uuid",
    "ena",
    "pseudogene.org",
    "ccds_id",
    "bioparadigms_slc",
    "prev_name",
    "uniprot_ids",
    "cosmic",
    "location",
    "date_symbol_changed",
    "cd",
    "_version_",
    "homeodb",
    "curator_notes",
    "locus_group",
    "pubmed_id",
    "rgd_id",
    "mane_select",
    "refseq_accession",
    "ucsc_id",
    "alias_name",
    "vega_id",
    "lncipedia",
    "date_approved_reserved",
    "lncrnadb",
    "lsdb",
    "agr"
]


def search_hugo(search_by, to_be_searched):
    base_url = "https://rest.genenames.org/search/"
    try:
        r = requests.get(base_url+search_by+'/'+to_be_searched, headers={"Accept":"application/json"}, timeout=120) 
        data = r.json().get("response", {}).get("docs", [])
        return pd.json_normalize(data)

    except Exception as e:
        print("\n\n search_failed \n\n")
        print(str(e))

    return None


def fetch_hugo(fetch_by, to_be_fetched, path_to_download, throw=True):
    if check_string(to_be_fetched):
        return None

    base_url = "https://rest.genenames.org/fetch/"
    try:
        r = requests.get(base_url+fetch_by+'/'+to_be_fetched, headers={"Accept":"application/json"}, timeout=120)
        data = r.json().get("response", {}).get("docs", [])
        df = pd.json_normalize(data)

        if df.empty:
            return None

        name = df["ensembl_gene_id"].iloc[0]
        hugo_gen_dir = os.path.join(path_to_download, name)
        os.makedirs(hugo_gen_dir, exist_ok=True)
        file_path = os.path.join(hugo_gen_dir, "hgnc_data.tsv")
        df.to_csv(file_path ,sep='\t', index=False)

        send_message(f"got {to_be_fetched}", 0, "hgnc")
        send_message(1,1,"hgnc")
        sleep(0.1)
        return df

    except Exception as e:
        if throw:
            send_message(f" - fetch failed {fetch_by}, {to_be_fetched}, {path_to_download}\n{str(e)}\n")

    return None


def load_table(path_to_table):
    df = pd.DataFrame()
    try:
        if isinstance(path_to_table, list):
            df_list = []
            for tsv in path_to_table:
                dff = load_table(tsv)
                if dff.empty:
                    continue
                df_list.append(dff)
            if df_list:
                df = pd.concat(df_list)
        elif path_to_table.endswith("csv"):
            df = pd.read_csv(path_to_table, sep=',', dtype=str)
        elif path_to_table.endswith("tsv"):
            df = pd.read_csv(path_to_table, sep='\t', dtype=str)
    except Exception as e:
        send_message(f" - couldn't load {path_to_table}\n{str(e)}\n")

    return df

def load_hgnc_symbol_check(path_to_hgnc_symbol):
    df = load_table(path_to_hgnc_symbol)

    if df.empty:
        return None, None

    df["Approved symbol"] = df["Approved symbol"].fillna("NO_SYMBOL")
    df["Approved name"] = df["Approved name"].fillna("NO_NAME")
    df["HGNC ID"] = df["HGNC ID"].fillna("NO_ID")
    df = df.drop_duplicates(ignore_index=True)
    vc = df["Input"].value_counts()

    unambiguous = df[ 
        df["Input"] == df["Approved symbol"] | #Input is identical to the approved symbol
        df["Input"].isin(vc.index[vc.eq(1)]) &  #there is only one approved symbol found, and it has at least a symbol, name or id
        (
            (df["Approved symbol"] != "NO_SYMBOL") |
            (df["Approved name"] != "NO_NAME") |
            (df["HGNC ID"] != "NO_ID")
        )
    ]

    ambiguous = df[ ~df["Input"].isin(unambiguous["Input"]) ]
    return unambiguous, ambiguous


def get_rename_dict(look_up_table):
    inputs = look_up_table["Input"].tolist()
    names = look_up_table["Approved symbol"].tolist()
    return {
        input_name : gene_name
        for input_name, gene_name in zip(inputs, names)
        }


def build_tables(path_to_HGNC, result_path=None):
    resulttsv = search_for_files(path_to_HGNC, "hgnc_data", "tsv")
    offline_data = load_table(resulttsv)
    if offline_data.empty:
        return None

    return offline_data


def fill_checked(already_checked, df):
    if df is None:
        return already_checked
    if df.empty:
        return already_checked
    sub_df = df[["ensembl_gene_id", "symbol", "name"]].dropna().drop_duplicates(ignore_index=True)
    for gene_id in sub_df["ensembl_gene_id"].unique():
        if not check_string(gene_id):
            already_checked["ensembl_gene_id"].add(gene_id)
    for symbol in sub_df["symbol"].unique():
        if not check_string(symbol):
            already_checked["symbol"].add(symbol)
    for name in sub_df["name"].unique():
        if not check_string(name):
            already_checked["name"].add(name)
    return already_checked

def send_gene_id(hgnc_send, df):
    if not check_string(df):
        hgnc_send.put(df)
    elif isinstance(df, pd.DataFrame):
        for ensembl_id in df["ensembl_gene_id"].unique():
            if check_string(ensembl_id):
                continue
            hgnc_send.put(ensembl_id)

def send_or_not(hgnc_send, df, colum, value):
    if check_string(value):
        return False
    if df is None:
        return False
    found_ = df[(df[colum] == value)]
    if found_.empty:
        return False
    else:
        send_gene_id(hgnc_send, df)
        return True

def fetch_and_check(hgnc_send, already_checked, data_path, to_look_for, value, throw=True):
    df = fetch_hugo(to_look_for, value, data_path, throw)
    if df is None:
        return True, already_checked
    else:
        send_gene_id(hgnc_send, df)
        return False, fill_checked(already_checked, df)

def fetch_disgnet(hgnc_send, disgnet_df, data_path, result_path, unambigious, download):
    already_checked = {
        "ensembl_gene_id": set(),
        "symbol" : set(),
        "name": set(), 
        }
    advanced_search = set()
    offline_data = None
    if not download:
        offline_data = build_tables(data_path)
        already_checked = fill_checked(already_checked, offline_data)

    if disgnet_df.empty:
        return advanced_search

    disgnet_unique = disgnet_df["geneEnsemblIDs"].unique()
    send_message(disgnet_unique.shape[0],2,"hgnc")

    for ensembl_gene_id in disgnet_unique:
        send_gene_id(hgnc_send, ensembl_gene_id)
        if ensembl_gene_id not in already_checked["ensembl_gene_id"]:
            missed, already_checked = fetch_and_check(hgnc_send, already_checked, data_path, "ensembl_gene_id", ensembl_gene_id)
            if missed:
                advanced_search.add(ensembl_gene_id)
                continue
        send_message(f"finished {ensembl_gene_id}", 0, "hgnc")
        send_message(1,1,"hgnc")

    disgnet_unique = disgnet_df["gene_symbol"].unique()
    send_message(disgnet_unique.shape[0],2,"hgnc")
    symbol_dict = get_rename_dict(unambigious)

    for symbol in disgnet_unique:
        if symbol in unambigious["Input"]:
            symbol = symbol_dict[symbol] #disgnet only has clean gene symbols
        if not send_or_not(hgnc_send, offline_data, "symbol", symbol):
            advanced_search.add(symbol)
        if symbol not in already_checked["symbol"]:
            missed, already_checked = fetch_and_check(hgnc_send, already_checked, data_path, "symbol", symbol)
            if missed:
                missed, already_checked = fetch_and_check(hgnc_send, already_checked, data_path, "name", symbol, False)
                if missed:
                    advanced_search.add(symbol)
        send_message(f"finished {symbol}", 0, "hgnc")
        send_message(1,1,"hgnc")
    return advanced_search, already_checked

def fetch_amigo(hgnc_send, hgnc_receive, already_checked, advanced_search, data_path, unambiguouse, download):
    count_dict = Counter()
    offline_data = None
    if not download:
        offline_data = build_tables(data_path)
        already_checked = fill_checked(already_checked, offline_data)
    hgnc_send.put("amigo")
    while (True):
        try:
            if not hgnc_receive.empty():
                amigo_df = hgnc_receive.get()
                if isinstance(amigo_df, str) and amigo_df == "finished" or count_dict["NO_ID"] > 60:
                    break
                if amigo_df.empty:
                    continue

                """
                Looking up BIOENTITY_LABEL aka Symbol
                """
                #df["bioentity_name", "bioentity_label"]
                #clean up name
                amigo_unique = amigo_df["bioentity_label"].unique()
                send_message(amigo_unique.shape[0],2,"hgnc")
                symbol_dict = get_rename_dict(unambiguouse)

                for symbol in amigo_unique:
                    if symbol in unambiguouse["Input"]:
                        symbol = symbol_dict[symbol]
                    if not send_or_not(hgnc_send, offline_data, "symbol", symbol):
                        advanced_search.add(symbol)
                    if symbol not in already_checked["symbol"]:
                        missed, already_checked = fetch_and_check(hgnc_send, already_checked, data_path, "symbol", symbol)
                        if missed:
                            missed, already_checked = fetch_and_check(hgnc_send, already_checked, data_path, "name", symbol, False)
                            if missed:
                                advanced_search.add(symbol)
                    send_message(f"finished {symbol}", 0, "hgnc")
                    send_message(1,1,"hgnc")

                amigo_unique = amigo_df["bioentity_name"].unique()
                send_message(amigo_unique.shape[0],2,"hgnc")
                for name in amigo_unique:
                    if not send_or_not(hgnc_send, offline_data, "name", name):
                        advanced_search.add(name)
                    if name not in already_checked["name"]:
                        missed, already_checked = fetch_and_check(hgnc_send, already_checked, data_path, "name", name)
                        if missed:
                            missed, already_checked = fetch_and_check(hgnc_send, already_checked, data_path, "symbol", name, False)
                            if missed:
                                advanced_search.add(name)
                    send_message(f"finished {name}", 0, "hgnc")
                    send_message(1,1,"hgnc")
            else:
                count_dict["NO_ID"] +=1
                sleep(180)

        except Exception as e:
            send_message(f" - something broke in hgnc loop\n{str(e)}\n")
            count_dict["NO_ID"] +=1
            sleep(180)

    return advanced_search, already_checked

def clean_advanced_search(advanced_search, already_checked, data_path, result_path):
    already_checked = fill_checked(already_checked, build_tables(data_path))
    final_search = []
    for item in advanced_search:
        if item in already_checked["ensembl_gene_id"] or item in already_checked["symbol"] or item in already_checked["name"] or check_string(item):
            continue
        final_search.append(item)
    if final_search:
        file_path = os.path.join(result_path, "hgnc_not_found.txt")
        with open(file_path, "a") as f:
            for item in final_search:
                f.write(item)
        return file_path
    return None

def download_hgnc_data(hgnc_receive, hgnc_send, hgnc_config):
    send_message("started", 0, "hgnc")
    disgnet_df, symbol_checker, data_path, result_path, download = hgnc_config
    unambiguouse, ambiguouse = load_hgnc_symbol_check(symbol_checker)

    advanced_search, already_checked = fetch_disgnet(hgnc_send, disgnet_df, data_path, result_path, unambiguouse, download)
    advanced_search, already_checked = fetch_amigo(hgnc_send, hgnc_receive, already_checked, advanced_search, data_path, unambiguouse, download)

    hgnc_send.put("finished")

    file_path = clean_advanced_search(advanced_search, already_checked, data_path, result_path)
    if file_path is not None:
        send_message(f"hgnc_coulnd't identify these {file_path}", 0, "hgnc")

    send_message("waiting for clean_up", 0, "hgnc")
    counter = 0
    final_list = set()
    while (True):
        try:
            if not hgnc_send.empty():
                ensembl_id = hgnc_receive.get()

                if check_string(ensembl_id):
                    continue

                if ensembl_id == "finished" or counter > 60:
                    hgnc_send.close()
                    break

                final_list.add(ensembl_id)

            else:
                counter += 1
                sleep(180)

        except Exception as e:
            counter +=1
            if counter > 60:
                send_message(f" - something broke in hgnc final loop\n{str(e)}\n")
                break
            sleep(300)


    df = build_tables(data_path)
    if df is not None and final_list:
        #all relevant look_up_data
        final_hugo = df[df["ensembl_gene_id"].isin(final_list)].drop_duplicates(ignore_index=True)
        best_hugo = os.path.join(result_path, "hgnc_df.tsv")
        final_hugo.to_csv(best_hugo, sep="\t", index=False)

        look_up_hugo = final_hugo[["ensembl_gene_id", "hgnc_id", "symbol", "name"]]
        look_up_hugo["Input"] = final_hugo["symbol"]
        look_up_hugo = look_up_hugo.rename(columns={"symbol": "Approved symbol", "name": "Approved name"})
        best_look_up_hugo = os.path.join(result_path, "look_up_hugo.tsv")
        look_up_hugo = look_up_hugo[["Input", "Approved symbol", "Approved name", "ensembl_gene_id", "hgnc_id"]].drop_duplicates(ignore_index=True)
        look_up_hugo.to_csv(best_look_up_hugo, sep="\t", index=False)
    send_message("finished", 0, "hgnc")