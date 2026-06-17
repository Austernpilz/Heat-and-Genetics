import requests
import pandas as pd
import json
import os
import numpy as np
from time import sleep
from datetime import datetime
#http://www.ensembl.org/biomart/martview/ad4dbf2f9ae74dbf5b9cda391d970be9?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id_version|hsapiens_gene_ensembl.default.feature_page.description|hsapiens_gene_ensembl.default.feature_page.start_position|hsapiens_gene_ensembl.default.feature_page.end_position|hsapiens_gene_ensembl.default.feature_page.chromosome_name|hsapiens_gene_ensembl.default.feature_page.hgnc_id|hsapiens_gene_ensembl.default.feature_page.entrezgene_id|hsapiens_gene_ensembl.default.feature_page.uniprot_gn_id&FILTERS=hsapiens_gene_ensembl.default.filters.hgnc_id."HGNC:5970"&VISIBLEPANEL=resultspanel

from src.ensembl.extending_variants import extend_data
def get_config(config):
    storage = config.get("absolute_file_paths")
    data_path = os.path.join(storage.get("data"), "ensembl")
    os.makedirs(data_path, exist_ok=True)
    download = config.get("flags").get("download_data")
    return (data_path, download)


def fetch_from_ensembl(ensembl_id, path_to_ensembl):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{ensembl_id}?expand=1"
    try:
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, timeout=180)
        if r.status_code != 200:
            print("failed to load ",ensembl_id)
        decoded = r.json()
        id_dir = os.path.join(path_to_ensembl, ensembl_id)
        os.makedirs(id_dir, exist_ok=True)
        p = os.path.join(id_dir, "ensembl_data.json")
        with open(p, 'w') as file:
            json.dump(decoded, file)
        print(f"{datetime.now().strftime('%H%M')} ensembl got: {ensembl_id}")
        return p

    except Exception as e:
        print(f"\n\n fetch failed")
        print(f"{path_to_ensembl}, {ensembl_id}")
        print(str(e))

    return None

#this function downloads the data if necessary
def get_from_path(ensembl_id, path_to_ensembl):
    path_to_id = os.path.join(path_to_ensembl, ensembl_id, "ensembl_data.json")
    if os.path.isfile(path_to_id):
        try:
            print(f"{datetime.now().strftime('%H%M')} ensembl got: {ensembl_id}")
            return pd.read_json(path_to_id)
        except Exception as e:
            print(str(e))

    path_to_data = fetch_from_ensembl(ensembl_id, path_to_ensembl)
    if path_to_data is not None:
        return pd.read_json(path_to_data)
    else:
        return None


def download_data(ensembl_receive, ensembl_send, ensembl_config):
    already_visited = set()
    data_path, download = ensembl_config
    ensembl_id = ""
    print("starting ensembl")
    while (True):
        try:
            ensembl_id = ensembl_receive.recv()
            if ensembl_id == "finished":
                ensembl_receive.close()
                break
            elif ensembl_id in already_visited or not isinstance(ensembl_id, str) or ensembl_id is np.nan:
                continue
            else:
                already_visited.add(ensembl_id)
                ensembl_send.send(ensembl_id)
        except Exception as _:
            sleep(90) #to build up the previous processes

        if download:
            _ = fetch_from_ensembl(ensembl_id, data_path)
            sleep(0.1)
        else:
            _ = get_from_path(ensembl_id, data_path)

    ensembl_send.send("finished")
    ensembl_send.close()
    print("ensembl thread done")


