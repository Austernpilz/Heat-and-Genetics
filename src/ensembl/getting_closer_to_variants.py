import pandas as pd
import os
from time import sleep
from collections import Counter
import json
#http://www.ensembl.org/biomart/martview/ad4dbf2f9ae74dbf5b9cda391d970be9?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id_version|hsapiens_gene_ensembl.default.feature_page.description|hsapiens_gene_ensembl.default.feature_page.start_position|hsapiens_gene_ensembl.default.feature_page.end_position|hsapiens_gene_ensembl.default.feature_page.chromosome_name|hsapiens_gene_ensembl.default.feature_page.hgnc_id|hsapiens_gene_ensembl.default.feature_page.entrezgene_id|hsapiens_gene_ensembl.default.feature_page.uniprot_gn_id&FILTERS=hsapiens_gene_ensembl.default.filters.hgnc_id."HGNC:5970"&VISIBLEPANEL=resultspanel

from src.ensembl.fetch_ensembl import fetch_from_ensembl
from src.ensembl.extending_variants import extend_data
from src.helpers.std_out import send_message
from src.helpers.folder_magic import check_string

def check_if_exists(path_to_ensembl, ensembl_id):
    path_to_id = os.path.join(path_to_ensembl, ensembl_id, "ensembl_data.json")
    return os.path.isfile(path_to_id)

def get_from_path(ensembl_id, path_to_ensembl):
    if check_if_exists(path_to_ensembl, ensembl_id):
        try:
            path_to_id = os.path.join(path_to_ensembl, ensembl_id, "ensembl_data.json")
            data = {}
            with open(path_to_id) as f:
                data = json.load(f)
            return data
        except Exception as e:
            send_message(f" - couldn't read {ensembl_id} from {path_to_ensembl} \n{str(e)}\n")

    return None


def download_data(ensembl_receive, ensembl_send, ensembl_config):
    send_message("starting", 0, "ensembl")
    data_path, top_genes, download = ensembl_config
    already_visited = set()
    amigo_count = Counter()
    c = False
    while (True):
        try:
            if not ensembl_receive.empty():
                ensembl_id = ensembl_receive.get()

                if check_string(ensembl_id):
                    continue

                if ensembl_id == "finished" or amigo_count["time_out"] > 60:
                    break

                if ensembl_id == "amigo":
                    c = True
                    continue
                if c:
                    amigo_count[ensembl_id] += 1
                else:
                    ensembl_send.put(ensembl_id)
                    amigo_count[ensembl_id] += 1000
                    send_message(1, 2, "gnomad")

                if ensembl_id in already_visited:
                    continue

                if not download and check_if_exists(data_path, ensembl_id):
                    send_message(1,1,"ensembl")
                    send_message(f"got {ensembl_id}",0,"ensembl")
                    already_visited.add(ensembl_id)
                    continue

                _ = fetch_from_ensembl(ensembl_id, data_path)
                already_visited.add(ensembl_id)

            else:
                amigo_count["time_out"] += 1
                sleep(180)
        except Exception as _:
            amigo_count["time_out"] += 1
            sleep(180) #to build up the previous processes

    top_amigo = amigo_count.most_common()
    last_count, i = 0, 0
    for ensembl_id, count in top_amigo:
        if check_string(ensembl_id):
            continue
        i += 1
        if i >= top_genes and last_count > count:
            break

        ensembl_send.put(ensembl_id)
        ensembl_receive.put(ensembl_id)
        send_message(1, 2, "gnomad")
        last_count = count

    ensembl_send.put("finished")
    ensembl_receive.put("finished")

    send_message("finished", 0, "ensembl")