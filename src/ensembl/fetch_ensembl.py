import json
import os
#from time import sleep
from datetime import datetime

import requests
import pandas as pd
import numpy as np


def fetch_from_ensembl(ensembl_id, path_to_ensembl):
    if ensembl_id is None:
        return
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{ensembl_id}?expand=1"
    try:
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, timeout=600)
        if r.status_code != 200:
            print("failed to load ",ensembl_id)
        if not r.ok:
            r.raise_for_status()

        decoded = r.json()
        id_dir = os.path.join(path_to_ensembl, ensembl_id)
        os.makedirs(id_dir, exist_ok=True)
        p = os.path.join(id_dir, "ensembl_data.json")
        with open(p, 'w') as file:
            json.dump(decoded, file)
        print(f"{datetime.now().strftime('%H%M')} ensembl got: {ensembl_id}")
        #sleep(0.1)
        return p

    except Exception as e:
        print(f"\n\n fetch failed")
        print(f"{path_to_ensembl}, {ensembl_id}")
        print(str(e))

    return None


def fetch_hgvs_data(args):
    path_to_data, hgvs, download = args
    if hgvs is None:
        return

    id_dir = os.path.join(path_to_data, hgvs)
    os.makedirs(id_dir, exist_ok=True)
    #ts = datetime.now().strftime("%Y%m%dT%H")
    p0 = os.path.join(path_to_data, f"variant_by_hgvs_{hgvs}.json")
    p1 = os.path.join(id_dir, f"variant_by_hgvs_{hgvs}.json")
    if not (os.path.isfile(p0) | os.path.isfile(p1)) | download:
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
            ext = f"/vep/homo_sapiens/hgvs/{hgvs}?{options}"
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, timeout=600)

            if not r.ok:
                r.raise_for_status()

            decoded = r.json()
            with open(p1, 'w') as file:
                json.dump(r.json(), file)

            print(f"{datetime.now().strftime('%H%M')} VEP got {hgvs}")
            #sleep(0.1)

        except Exception as e:
            print("\n\n hgvs fetch failed")
            print(str(e))


def fetch_rsid_data(args):
    path_to_data, rsid, download = args
    if rsid is None:
        return

    id_dir = os.path.join(path_to_data, rsid)
    os.makedirs(id_dir, exist_ok=True)

    p0 = os.path.join(path_to_data, f"variant_by_rsID_{rsid}.json")
    p1 = os.path.join(id_dir, f"variant_by_rsID_{rsid}.json")
    if not (os.path.isfile(p0) | os.path.isfile(p1)) | download:
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
            with open(p1, 'w') as file:
                json.dump(r.json(), file)

            print(f"{datetime.now().strftime('%H%M')} VEP got {rsid}")
            #sleep(0.1)

        except Exception as e:
            print(str(e))
            print("\n\n rsid fetch failed")


def fetch_pop_data(args):
    path_to_data, rsid, download = args
    if rsid is None:
        return
    id_dir = os.path.join(path_to_data, rsid)
    os.makedirs(id_dir, exist_ok=True)
    #ts = datetime.now().strftime("%Y%m%dT%H")
    p0 = os.path.join(path_to_data, f"populations_{rsid}.json")
    p1 = os.path.join(id_dir, f"populations_{rsid}.json")
    if not (os.path.isfile(p0) | os.path.isfile(p1)) | download:
        try:
            server = "https://rest.ensembl.org"
            options = {
                "phenotypes" : 1, 
                "pops" : 1, 
                }
            ext = f"/variation/homo_sapiens/{rsid}?{options}"
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, timeout=600)

            if not r.ok:
                r.raise_for_status()

            decoded = r.json()
            with open(p1, 'w') as file:
                json.dump(r.json(), file)

            print(f"{datetime.now().strftime('%H%M')} VEP got POP data from {rsid}")
            #sleep(0.1)

        except Exception as e:
            print("\n\n pop fetch failed")
            print(str(e))


def translate_to_rsid(args):
    path_to_data, hgvs, download = args
    if hgvs is None:
        return []

    #ts = datetime.now().strftime("%Y%m%dT%H")
    id_dir = os.path.join(path_to_data, hgvs)
    os.makedirs(id_dir, exist_ok=True)
    p0 = os.path.join(path_to_data, f"aternative_names_for_{hgvs}.json")
    p1 = os.path.join(id_dir, f"aternative_names_for_{hgvs}.json")
    if not (os.path.isfile(p0) | os.path.isfile(p1)) | download:

        try:
            server = "https://rest.ensembl.org"
            ext = f"/variant_recoder/human/{hgvs}"

            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, timeout=600)

            if not r.ok:
                r.raise_for_status()

            decoded = r.json()
            with open(p1, 'w') as file:
                json.dump(decoded, file)

            print(f"{datetime.now().strftime('%H%M')} VEP got {hgvs} to rsid")
            #sleep(0.1)
            rsids = []
            #print(decoded)
            for found in decoded:
                for base, possible_ids in found.items():
                    rsids += possible_ids.get("id", [])
            return rsids
        except Exception as e:
            print("\n\n translation failed")
            print(str(e))
    elif os.path.isfile(p):
        rsids = []
        with open(p, 'r') as f:
            translate_hgvs = json.load(f)
            for found in translate_hgvs:
                for base, possible_ids in found.items():
                    rsids += possible_ids.get("id", [])
        return rsids

    return []

