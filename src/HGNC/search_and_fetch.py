import requests
import json 
import os
import pandas as pd
import numpy as np
from io import StringIO
from time import sleep
from datetime import datetime

from collections import Counter

from src.helpers.folder_magic import search_for_files

# 10/s is the rate at which i can shoot requests 

# can use this with search and fetch
# i want hgnc_id and maybe ens
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

def get_config(config):
    storage = config.get("absolute_file_paths")
    data_path = os.path.join(storage.get("data"), "HGNC")
    hgnc = []
    hgnc_symbol_check = config.get("relative_file_paths").get("hgnc_symbol_check")
    for hgs in hgnc_symbol_check:
        hgnc.append(os.path.join(data_path, hgs))

    download = config.get("flags").get("download_data")
    top_genes = config.get("flags").get("top_genes")
    hgnc_data = config.get("relative_file_paths").get("hgnc_data")
    data_path = os.path.join(data_path, hgnc_data)
    os.makedirs(data_path, exist_ok=True)
    return (data_path, hgnc, download)


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


def fetch_hugo(fetch_by, to_be_fetched, path_to_download):
    base_url = "https://rest.genenames.org/fetch/"
    try:
        r = requests.get(base_url+fetch_by+'/'+to_be_fetched, headers={"Accept":"application/json"}, timeout=60)
        data = r.json().get("response", {}).get("docs", [])
        df = pd.json_normalize(data)

        if df.empty:
            return None

        #print(df.head())
        name = df["ensembl_gene_id"].iloc[0]
        hugo_gen_dir = os.path.join(path_to_download, name)
        os.makedirs(hugo_gen_dir, exist_ok=True)
        file_path = os.path.join(hugo_gen_dir, "hgnc_data.tsv")
        df.to_csv(file_path ,sep='\t', index=False)

        print(f"{datetime.now().strftime('%H%M')} hgnc got: {to_be_fetched}")
        # print("saved  to", file_path)
        sleep(0.2)
        return df

    except Exception as e:
        print("\n\n fetch failed \n\n")
        print(f"{fetch_by}, {to_be_fetched}, {path_to_download}")
        print(str(e))

    return None


def load_hgnc_symbol_check(paths_to_hgnc_symbol):
    df_list = []
    for symbol_check in paths_to_hgnc_symbol:
        try:
            df_list.append(pd.read_csv(symbol_check, sep=',', dtype=str))
        except Exception as e:
            print(str(e))
            print("couldn't load ", symbol_check)

    if not df_list:
        return None

    df = pd.concat(df_list).drop_duplicates()
    df["Approved symbol"] = df["Approved symbol"].fillna("NO SYMBOL")
    df["Approved name"] = df["Approved name"].fillna("NO NAME")
    df["HGNC ID"] = df["HGNC ID"].fillna("NO ID")
    vc = df["Input"].value_counts()

    #there is only one approved symbol found, and it has at least a symbol, name or id
    unambiguous_0 = df[ 
        df["Input"].isin(vc.index[vc.eq(1)]) & 
        (
            (df["Approved symbol"] != "NO SYMBOL") |
            (df["Approved name"] != "NO NAME") |
            (df["HGNC ID"] != "NO ID")
            )
        ]

    #Input is identical to the approved symbol
    unambiguous_1 = df[
        ~df["Input"].isin(unambiguous_0["Input"]) & 
        df["Input"] == df["Approved symbol"]
        ]

    ambiguous = df[
        ~df["Input"].isin(unambiguous_0["Input"]) & 
        ~df["Input"].isin(unambiguous_1["Input"])
        ]

    unambiguous = pd.concat([unambiguous_0, unambiguous_1])

    return unambiguous, ambiguous


def get_rename_dict(look_up_table):
    inputs = look_up_table["Input"].tolist()
    names = look_up_table["Approved symbol"].tolist()
    return {
        input_name : gene_name
        for input_name, gene_name in zip(inputs, names)
        }


def get_tables (path_to_HGNC):
    resulttsv = search_for_files(path_to_HGNC, "hgnc_data", "tsv")

    offline_data = []
    for tsv in resulttsv:
        try:
            df = pd.read_csv(tsv, sep='\t')
            if df.empty:
                continue
            offline_data.append(df)
        except Exception as _:
            continue
    if offline_data:
        df = pd.concat(offline_data).drop_duplicates()
        #print(df.head(10))
        return df
    else:
        return None


def fetch_disgnet(hgnc_send, disgnet_df, data_path, unambiguouse, download):
    offline_data = get_tables(data_path)
    already_checked = {
        "ensembl_gene_id": [],
        "symbol" : [],
        "name": [], 
        }
    advanced_search = []
    if not download and offline_data is not None:
        sub_df = offline_data[["ensembl_gene_id", "symbol", "name"]].dropna().drop_duplicates()
        already_checked["ensembl_gene_id"] += sub_df["ensembl_gene_id"].tolist()
        already_checked["symbol"] += sub_df["symbol"].tolist()
        already_checked["name"] += sub_df["name"].tolist()

    if disgnet_df.empty:
        return advanced_search

    for ensembl_gene_id in disgnet_df["geneEnsemblIDs"].unique():
        if (
            ensembl_gene_id == "NO ID" or 
            ensembl_gene_id == np.nan or 
            ensembl_gene_id == "nan" or 
            not isinstance(ensembl_gene_id, str) or
            not "ENS" in ensembl_gene_id
        ):
            continue
        hgnc_send.send(ensembl_gene_id)

        if download or ensembl_gene_id not in already_checked["ensembl_gene_id"]:
            df = fetch_hugo("ensembl_gene_id", ensembl_gene_id, data_path)
            if df is None:
                advanced_search.append(ensembl_gene_id)
                continue
            else:
                sub_df = df[["ensembl_gene_id", "symbol", "name"]].dropna().drop_duplicates()
                already_checked["ensembl_gene_id"] += sub_df["ensembl_gene_id"].tolist()
                already_checked["symbol"] += sub_df["symbol"].tolist()
                already_checked["name"] += sub_df["name"].tolist()

    symbol_dict = get_rename_dict(unambiguouse)
    #print(symbol_dict)
    for symbol in disgnet_df["gene_symbol"].unique():
        if symbol in unambiguouse["Input"]:
            symbol = symbol_dict[symbol] #disgnet only has clean gene symbols
            #print(genes)
        elif offline_data is not None:
            #try to find the ensembl id
            ensembl_id_list = offline_data[(offline_data["symbol"] == symbol)]
            if ensembl_id_list.empty:
                advanced_search.append(symbol)
                continue
            else:
                for ensembl_id in ensembl_id_list["ensembl_gene_id"].unique():
                    hgnc_send.send(id)
        if download or symbol not in already_checked["symbol"]:
            df = fetch_hugo("symbol", symbol, data_path)
            if df is None:
                df = fetch_hugo("name", symbol, data_path)
                if df is None:
                    advanced_search.append(symbol)
                    continue

            if df is not None:
                sub_df = df[["ensembl_gene_id", "symbol", "name"]].dropna().drop_duplicates()
                already_checked["ensembl_gene_id"] += sub_df["ensembl_gene_id"].tolist()
                already_checked["symbol"] += sub_df["symbol"].tolist()
                already_checked["name"] += sub_df["name"].tolist()
                for genes in sub_df["ensembl_gene_id"].unique():
                    hgnc_send.send(genes)

    return advanced_search

def fetch_amigo(hgnc_send, hgnc_receive, data_path, unambiguouse, advanced_search, download):
    name, symbol = "", ""
    offline_data = get_tables(data_path)
    df_list = []
    count_dict = Counter()
    already_checked = {
        "ensembl_gene_id" : [], 
        "symbol" : [], 
        "name" : []
    }
    if not download and offline_data is not None:
        sub_df = offline_data[["ensembl_gene_id", "symbol", "name"]].dropna().drop_duplicates()
        already_checked["ensembl_gene_id"] = sub_df["ensembl_gene_id"].tolist()
        already_checked["symbol"] = sub_df["symbol"].tolist()
        already_checked["name"] = sub_df["name"].tolist()

    symbol_dict = get_rename_dict(unambiguouse)
    while (True):
        try:
            if hgnc_receive.poll(timeout=180):
                name, symbol = hgnc_receive.recv()
            else:
                count_dict["NO ID"] +=1
            if name == "finished" and symbol == "finished" or count_dict["NO ID"] > 10:
                hgnc_receive.close()
                break
        except Exception as e:
            print(str(e))
            count_dict["NO ID"] +=1
            sleep(180)

        """
        Looking up BIOENTITY_LABEL aka Symbol
        """
        #clean up name
        if symbol in unambiguouse["Input"]:
            symbol = symbol_dict[symbol]
            #print(symbol)

            # #it is either already send or I find it in offline, data
            # #or both, in which case ensembl needs to filter out the double
            # ensembl_id = offline_data[(offline_data["symbol"] == symbol)]
            # if not ensembl_id.empty:
            #     for id in ensembl_id["ensembl_gene_id"].unique():
            #         if id in count_dict.keys():
            #             count_dict[id] += 1
            #         else:
            #             count_dict[id] = 1
            # continue
        if download or symbol not in already_checked["symbol"]:
            df = fetch_hugo("symbol", symbol, data_path)
            if df is None:
                df = fetch_hugo("name", symbol, data_path)
            if df is not None:
                sub_df = df[["ensembl_gene_id", "symbol", "name"]].dropna().drop_duplicates()
                already_checked["ensembl_gene_id"] += sub_df["ensembl_gene_id"].tolist()
                already_checked["symbol"] += sub_df["symbol"].tolist()
                already_checked["name"] += sub_df["name"].tolist()
                for ensembl_id in sub_df["ensembl_gene_id"].unique():
                    count_dict[ensembl_id] += 1
                continue
            else:
                advanced_search.append(symbol)

        elif offline_data is not None:
            ensembl_id_list = offline_data[(offline_data["symbol"] == symbol)]
            if ensembl_id_list.empty:
                advanced_search.append(symbol)
            else:
                for ensembl_id in ensembl_id_list["ensembl_gene_id"].unique():
                    count_dict[ensembl_id] += 1
                if not download:
                    continue
        """
        Looking up BIOENTITY NAME aka Name
        """
        # if name in already_checked["name"]:
        #     continue
            # ensembl_id = offline_data[(offline_data["name"] == name)]
            # if not ensembl_id.empty:
            #     for id in ensembl_id["ensembl_gene_id"].unique():
            #         hgnc_send.send(id)
            # continue
        if download or name not in already_checked["name"]:
            df = fetch_hugo("name", name, data_path)
            if df is None:
                df = fetch_hugo("symbol", name, data_path)
            if df is not None:
                sub_df = df[["ensembl_gene_id", "symbol", "name"]].dropna().drop_duplicates()
                already_checked["ensembl_gene_id"] += sub_df["ensembl_gene_id"].tolist()
                already_checked["symbol"] += sub_df["symbol"].tolist()
                already_checked["name"] += sub_df["name"].tolist()
                for ensembl_id in sub_df["ensembl_gene_id"].unique():
                    count_dict[ensembl_id] += 1
                continue
            else:
                advanced_search.append(name)

        elif unambiguouse[(unambiguouse["Approved name"] == name)].empty and offline_data is not None:
            ensembl_id_list = offline_data[(offline_data["name"] == name)]
            if ensembl_id_list.empty:
                advanced_search.append(name)
            else:
                for ensembl_id in ensembl_id_list["ensembl_gene_id"].unique():
                    count_dict[ensembl_id] += 1


    for ensembl_id, _ in count_dict.most_common(200):
        hgnc_send.send(ensembl_id)
    #amigo_disgnet_done = pd.concat(df_list + [already_visited]).drop_duplicates()
    return advanced_search


def download_hgnc_data(hgnc_receive, hgnc_send, disgnet_df, hgnc_config):
    data_path, hgnc, download = hgnc_config
    unambiguouse, ambiguouse = load_hgnc_symbol_check(hgnc)
    #print(unambiguouse)
    print("starting hugo")
    advanced_search = fetch_disgnet(hgnc_send, disgnet_df, data_path, unambiguouse, download)
    print("disgnet done")
    advanced_search = fetch_amigo(hgnc_send, hgnc_receive, data_path, unambiguouse, advanced_search, download)
    print("amigo done")
    hgnc_send.send("finished")
    hgnc_send.close()
    #TO-DO
    #when I'm here, Amigo has run through, so now i can check amigo and disgnet by hand for advanced search
    file_path = os.path.join(data_path, "advanced_search_open.txt")
    ts = datetime.now().strftime("%Y%m%dT%H")
    with open(file_path, "a+") as f:
        f.write(f"\n\n\n{str(ts)}")
        for item in advanced_search:
            f.write(str(item))
            f.write("\n")

    print("hgnc thread done")


# def clean_up(list_of_df, path_to_HGNC, save=True):
#     hgnc_symbol_check = load_simple_search_data(path_to_HGNC)
#     rest_to_be_cleaned, look_up_table = easy_clean_up(hgnc_symbol_check)
#     # do_complex_clean(list_of_df, rest_to_be_cleaned, look_up_table, path_to_HGNC)

#     if save:
#         save_look_up(look_up_table, path_to_HGNC)

#     replace_dict = get_rename_dict(look_up_table)
#     for df in list_of_df:
#         df.replace(to_replace=replace_dict, inplace=True)

#     return list_of_df, rest_to_be_cleaned



