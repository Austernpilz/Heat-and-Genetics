import requests
import pandas as pd
import json
import os
import numpy as np
from time import sleep
from datetime import datetime
#http://www.ensembl.org/biomart/martview/ad4dbf2f9ae74dbf5b9cda391d970be9?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id_version|hsapiens_gene_ensembl.default.feature_page.description|hsapiens_gene_ensembl.default.feature_page.start_position|hsapiens_gene_ensembl.default.feature_page.end_position|hsapiens_gene_ensembl.default.feature_page.chromosome_name|hsapiens_gene_ensembl.default.feature_page.hgnc_id|hsapiens_gene_ensembl.default.feature_page.entrezgene_id|hsapiens_gene_ensembl.default.feature_page.uniprot_gn_id&FILTERS=hsapiens_gene_ensembl.default.filters.hgnc_id."HGNC:5970"&VISIBLEPANEL=resultspanel

from src.ensembl.fetch_ensembl import fetch_from_ensembl
from src.ensembl.extending_variants import extend_data

def get_config(config):
    storage = config.get("absolute_file_paths")
    data_path = os.path.join(storage.get("data"), "ensembl")
    os.makedirs(data_path, exist_ok=True)
    download = config.get("flags").get("download_data")
    return (data_path, download)


#this function downloads the data if necessary
def get_from_path(ensembl_id, path_to_ensembl):
    path_to_id = os.path.join(path_to_ensembl, ensembl_id, "ensembl_data.json")
    if os.path.isfile(path_to_id):
        try:
            df = pd.read_json(path_to_id)
            print(f"{datetime.now().strftime('%H%M')} ensembl got: {ensembl_id}")
            return df
        except Exception as e:
            print(str(e))

    return None


def download_data(ensembl_receive, ensembl_send, ensembl_config):
    already_visited = set()
    data_path, download = ensembl_config
    ensembl_id = ""
    counter = 0
    print("starting ensembl")
    while (True):
        try:
            if ensembl_receive.poll(timeout=600):
                ensembl_id = ensembl_receive.recv()
            else:
                counter += 1
                ensembl_id = "NO ID"

            if ensembl_id in already_visited:
                already_visited.add(ensembl_id)
                continue

            elif ensembl_id == "finished" or (
                ensembl_id == "NO ID" and counter > 10
            ):
                ensembl_receive.close()
                break

            elif (
                not isinstance(ensembl_id, str) or 
                ensembl_id is np.nan or
                ensembl_id == "NO ID"
            ):
                continue

            else:
                already_visited.add(ensembl_id)
                ensembl_send.send(ensembl_id)

        except Exception as _:
            sleep(90) #to build up the previous processes

        if not download:
            df = get_from_path(ensembl_id, data_path)
            if df is not None:
                continue
        #no data found or force download is true
        _ = fetch_from_ensembl(ensembl_id, data_path)

    ensembl_send.send("finished")
    ensembl_send.close()
    print("ensembl thread done")