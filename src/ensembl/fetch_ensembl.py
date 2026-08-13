import json
import os
from time import sleep
from datetime import datetime

import requests
import pandas as pd
import numpy as np

from src.helpers.std_out import send_message
from src.helpers.folder_magic import check_string

def fetch_from_ensembl(ensembl_id, path_to_ensembl):
    if check_string(ensembl_id):
        return None

    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{ensembl_id}?expand=1"
    try:
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, timeout=600)

        if not r.ok:
            r.raise_for_status()

        decoded = r.json()
        id_dir = os.path.join(path_to_ensembl, ensembl_id)
        os.makedirs(id_dir, exist_ok=True)
        p = os.path.join(id_dir, "ensembl_data.json")
        with open(p, 'w') as file:
            json.dump(decoded, file)

        send_message(1,1,"ensembl")
        send_message(f"got {ensembl_id}",0,"ensembl")
        sleep(0.1)
        return p

    except Exception as e:
        send_message(f" - fetch failed {path_to_ensembl}, {ensembl_id}\n{str(e)}\n")

    return None


def fetch_hgvs_data(path_to_data, hgvs, download):
    if check_string(hgvs):
        return False

    os.makedirs(path_to_data, exist_ok=True)
    p = os.path.join(path_to_data, f"variant_by_hgvs_{hgvs}.json")
    if not os.path.isfile(p) or download:
        try:
            #print('loading_variation_data ', hgvs)
            server = "https://rest.ensembl.org"
            options = {
                "AlphaMissense" : 1, 
                "ClinPred" : 1, 
                "CADD" : 1, 
                "LoF" : 1, 
                "Paralogues": {"clnsig" : "ignore"}, 
                "REVEL" : 1,
                "SpliceAI" : 1, 
                }
            ext = f"/vep/homo_sapiens/hgvs/{hgvs}?"
            r = requests.get(server+ext, params=options, headers={ "Content-Type" : "application/json"}, timeout=600)

            if not r.ok:
                r.raise_for_status()

            decoded = r.json()
            with open(p, 'w') as file:
                json.dump(decoded, file)
            send_message(1,1,"vep")
            send_message(f"got {hgvs}",0,"vep")
            sleep(0.1)

        except Exception as e:
            send_message(f" - hgvs fetch failed\n{str(e)}\n")
            return False
    return True

def fetch_rsid_data(path_to_data, rsid, download):
    if isinstance(rsid, list):
        worked = False
        for r in rsid:
            worked |= fetch_rsid_data(path_to_data, r, download)
        return worked

    if check_string(rsid):
        return False

    p = os.path.join(path_to_data, f"variant_by_rsID_{rsid}.json")
    if not os.path.isfile(p) or download:
        try:
            server = "https://rest.ensembl.org"
            options = {
                "AlphaMissense" : 1, 
                "ClinPred" : 1, 
                "CADD" : 1, 
                "LoF" : 1, 
                "Paralogues": {"clnsig" : "ignore"}, 
                "REVEL" : 1,
                "SpliceAI" : 1, 
                }
            ext = f"/vep/homo_sapiens/id/{rsid}"
            r = requests.get(server+ext, params=options,  headers={ "Content-Type" : "application/json"}, timeout=600)

            if not r.ok:
                r.raise_for_status()

            decoded = r.json()
            with open(p, 'w') as file:
                json.dump(decoded, file)

            send_message(1,1,"vep")
            send_message(f"got {rsid}", 0, "vep")
            sleep(0.1)
        except Exception as e:
            send_message(f"- rsid fetch failed {rsid}\n{str(e)}\n")
            return False

    return True



def fetch_pop_data(path_to_data, rsid, download):
    if isinstance(rsid, list):
        results = []
        for r in rsid:
            results.extend(fetch_pop_data(path_to_data, r, download))
        return results

    if check_string(rsid):
        return []

    p = os.path.join(path_to_data, f"populations_{rsid}.json")
    if not os.path.isfile(p) or download:
        try:
            server = "https://rest.ensembl.org"
            ext = f"/variation/homo_sapiens/{rsid}?pops=1"
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, timeout=600)

            if not r.ok:
                r.raise_for_status()

            decoded = r.json()
            with open(p, 'w') as file:
                json.dump(decoded, file)

            send_message(f" - got POP data from {rsid}", 0, "vep")
            sleep(0.1)

        except Exception as e:
            send_message(f" - pop fetch failed {rsid}\n{str(e)}\n")
            return []

    return [p]


def translate_to_rsid(path_to_data, hgvs, download):
    if check_string(hgvs):
        return None

    p = os.path.join(path_to_data, f"aternative_names_for_{hgvs}.json")
    if not os.path.isfile(p) or download:
        try:
            server = "https://rest.ensembl.org"
            ext = f"/variant_recoder/human/{hgvs}"

            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, timeout=120)

            if not r.ok:
                r.raise_for_status()

            decoded = r.json()
            with open(p, 'w') as file:
                json.dump(decoded, file)

            send_message(f"got {hgvs} to rsid", 0, "vep")
            sleep(0.1)

        except Exception as e:
            send_message(f" - translation failed {hgvs}\n{str(e)}\n")
            return None

    return p

