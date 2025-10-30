import requests
import json 
import os
import pandas as pd
from io import StringIO
from time import sleep


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
    "ensembl_gene_id",
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

def get_gene_symbols(list_of_df):
    gene_symbols = ["bioentity_label", "gene", "gene_symbol"]
    return_set = set()
    for df in list_of_df:
        filtered_list = [col for col in gene_symbols if col in df.columns]
        return_set |= set(df[filtered_list[0]].unique())

    return return_set


def get_rename_dict(look_up_table):
    rename_dict = {} 
    inputs = look_up_table["Input"].unique().tolist()
    names = look_up_table["Approved symbol"].unique().tolist()
    for input, gene in zip(inputs, names):
        rename_dict[input] = gene
    return rename_dict


def load_simple_search_data(path_to_HGNC):
    path_data = os.path.join(path_to_HGNC, "data")
    df_list = []
    for entry in os.scandir(path_data):
        if not entry.is_file():
            continue
        elif not entry.name.endswith(".csv"):
            continue
        elif not entry.name.startswith("hgnc-symbol-check"):
            continue
        else:
            hgnc_csv = os.path.join(path_data, entry.name)
            try:
                df_list.append(pd.read_csv(hgnc_csv, sep=',', dtype=str))
            except Exception as e:
                print(str(e))
                print("couldn't load ", hgnc_csv)
                continue

    # for df in df_list:
    #     print(len(df["Input"].unique()))
    return pd.concat(df_list)

def search_hugo(search_by, to_be_searched):
    return_df = None

    try:
        base_url = "https://rest.genenames.org/search/"
        r = requests.get(base_url+search_by+'/'+to_be_searched, headers={"Accept":"application/json"}, timeout=60)
        decoded = r.json()
        data = decoded.get("response", {}).get("docs", [])
        return_df = pd.json_normalize(data)
    except Exception as e:
        print(str(e))

    return return_df

def fetch_hugo(fetch_by, to_be_fetched, path_to_download=False, download=False):
    return_df = None
    if download:
        try:
            base_url = "https://rest.genenames.org/fetch/"
            r = requests.get(base_url+fetch_by+'/'+to_be_fetched, headers={"Accept":"application/json"}, timeout=60)
            decoded = r.json()
            data = decoded.get("response", {}).get("docs", [])
            return_df = pd.json_normalize(data)
        except Exception as e:
            print(str(e))

    if path_to_download:
        try:
            hugo_gen_dir = os.path.join(path_to_download, "data", to_be_fetched)
            os.makedirs(hugo_gen_dir, exist_ok=True)
            file_path = os.path.join(hugo_gen_dir, "hgnc_data.tsv")
            if download:
                return_df.to_csv(file_path ,sep='\t')
            else:
                return_df = pd.read_csv(file_path, sep="\t", dtype=str)
        except Exception as e:
            print(str(e))

    return return_df

def save_look_up(look_up_table, path_to_HGNC):
    file_path = os.path.join(path_to_HGNC, "look_up_hugo.tsv")
    look_up_table.to_csv(file_path, sep='\t', index=False)

def test():
    path_to_HGNC = os.path.join(os.getcwd(), "HGNC")
    hgnc_symbol_check = load_simple_search_data(path_to_HGNC)
    app = hgnc_symbol_check["Approved symbol"].unique().tolist()
    inp = hgnc_symbol_check["Input"].unique().tolist()
    print(hgnc_symbol_check[~hgnc_symbol_check["Approved symbol"].isin(inp)]["Approved symbol"].unique().tolist())

def more_look_up_data(df):
    return df

def easy_clean_up(hgnc_symbol_check):
    #input is already approved_symbol
    simple = hgnc_symbol_check[hgnc_symbol_check["Input"] == hgnc_symbol_check["Approved symbol"]].drop_duplicates()
    #we need to look further
    rest = hgnc_symbol_check[~hgnc_symbol_check["Input"].isin(simple["Input"])].copy()

    rest["Approved symbol"].fillna("NO SYMBOL")
    rest["Approved name"].fillna("NO NAME")
    rest["HGNC ID"].fillna("NO ID")

    complex_search = []
    for _, row in rest.iterrows():
        input_r = row["Input"]
        if len(rest[rest["Input"] == input_r])>1:
            complex_search.append(input_r)
        elif row["Match type"] == "Unmatched":
            complex_search.append(input_r)
        elif (row["Approved symbol"] == "NO SYMBOL" and 
              row["HGNC ID"] == "NO ID" and 
              row["Approved name"] =="NO NAME"):
            complex_search.append(input_r)
        else:
            simple.loc[len(simple)] = row

    return more_look_up_data(rest[rest["Input"].isin(complex_search)]), simple


def load_look_up(path_to_HGNC):
    file_path = os.path.join(path_to_HGNC, "look_up_hugo.tsv")
    return pd.read_csv(file_path, sep='\t')

def load_HGNC(top_200_dataset, path_to_HGNC, load=True):
    look_up = load_look_up(path_to_HGNC)
    return load_data([top_200_dataset], look_up, path_to_HGNC, download=load)

def fetch_offline (gene_symbols, sub_look, load_by_symbol, path_to_HGNC):
    path_to_HGNC = os.path.join(path_to_HGNC, "data")

    if not os.path.isdir(path_to_HGNC):
        return [], []

    offline_data = []
    dir_to_visit = [path_to_HGNC]
    while dir_to_visit:
        current = dir_to_visit.pop(0)
        # print(current)
        try:
            for entry in os.scandir(current):

                if entry.name in ["bin", "include", "lib", "overview.txt", "data.tsv", "include_exclude.txt", "ensembl", "AmiGo2", "disgnet", "figures", "gnomAD"]:
                    continue

                if entry.is_dir():
                    dir_to_visit.append(os.path.join(current, entry.name))
                    continue

                if entry.name != "hgnc_data.tsv":
                    continue

                file_path = os.path.join(current, entry.name)
                offline_data.append( pd.read_csv(file_path, sep='\t') )

        except Exception as _:
            continue

    df = pd.concat(offline_data)
    symbol = []
    hgnc_id = []
    approved_name = []

    for symbol in load_by_symbol:
        if True:
            return

def load_data(list_of_df, look_up_table, path_to_HGNC, download=True):
    gene_symbols = get_gene_symbols(list_of_df)
    sub_look = look_up_table[look_up_table["Input"].isin(gene_symbols)]
    load_by_symbol = sub_look["Approved symbol"].unique().tolist()

    if not download:
        return fetch_offline(gene_symbols, sub_look, load_by_symbol, path_to_HGNC)

    HGNC_data = []
    next_to_download = []
    for symbol in load_by_symbol:
        df_hugo_symbol = fetch_hugo("symbol", symbol, path_to_HGNC, download=download)

        if df_hugo_symbol is None:
            next_to_download.append(symbol)
        elif df_hugo_symbol.empty:
            next_to_download.append(symbol)
        else:
            HGNC_data.append(df_hugo_symbol)
        sleep(0.1)
    load_by_id = sub_look[sub_look["Approved symbol"].isin(next_to_download)]
    load_by_id["HGNC ID"] = load_by_id["HGNC ID"].str.replace("HGNC:", '')
    # next_to_download = []

    for id in load_by_id["HGNC ID"].unique().tolist():
        df_hugo_symbol = fetch_hugo("hgnc_id", id, path_to_HGNC, download=download)
        if df_hugo_symbol.empty() or df_hugo_symbol is None:
            next_to_download.append(symbol)
        else:
            HGNC_data.append(df_hugo_symbol)
        sleep(0.1)

    load_by_name = sub_look[sub_look["HGNC ID"].isin(next_to_download)]
    for name in load_by_name["Approved name"].unique().tolist():
        df_hugo_symbol = fetch_hugo("name", name, path_to_HGNC, download=download)
        if df_hugo_symbol.empty() or df_hugo_symbol is None:
            next_to_download.append(symbol)
        else:
            HGNC_data.append(df_hugo_symbol)
        sleep(0.1)
    return pd.concat(HGNC_data), sub_look[sub_look["Approved name"].isin(next_to_download)].copy()


# def complex_search_offline(gene, list_of_df):
#     gene_symbols = ["bioentity_label", "gene", "gene_symbol"]
#     small_df = []
#     for df in list_of_df:
#         filtered_list = [col for col in gene_symbols if col in df.columns]
#         small_df.append(df[df[filtered_list] == gene])

#     return pd.concat(small_df, axis=1, ignore_index=True)


# def do_complex_clean(list_of_df, rest_to_be_cleaned, look_up_table, path_to_HGNC):
#     gene_symbols = rest_to_be_cleaned["Input"].unique().tolist()


#     for gene in gene_symbols:
#         set_of_identifiers = complex_search_offline(gene, list_of_df)
#         set_of_identifiers.fillna()
#     return



def clean_up(list_of_df, path_to_HGNC, save=True):
    hgnc_symbol_check = load_simple_search_data(path_to_HGNC)
    rest_to_be_cleaned, look_up_table = easy_clean_up(hgnc_symbol_check)
    # do_complex_clean(list_of_df, rest_to_be_cleaned, look_up_table, path_to_HGNC)

    if save:
        save_look_up(look_up_table, path_to_HGNC)

    replace_dict = get_rename_dict(look_up_table)
    for df in list_of_df:
        df.replace(to_replace=replace_dict, inplace=True)

    return list_of_df, rest_to_be_cleaned

# clean_up(None, os.path.join(os.getcwd(), "HGNC"))



#     HGNC_approved_symbol = ["Approved symbol", "symbol"]
# def build_look_up_table(list_of_df, path_to_HGNC):



