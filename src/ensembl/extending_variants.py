import requests
import pandas as pd
import json
import os
import numpy as np
from time import sleep
from datetime import datetime


def open_json(path):
    if not path.endswith("json"):
        return []
    try:
        v = []
        with open(path, 'r') as f:
            v = json.load(f).get("variants", [])
        return v
    except Exception as e:
        print(str(e))
        return []

def extend_data (VEP_receive, VEP_send, VEP_config):
    print("extending variants")
    files = []
    already_visited = set()
    data_path, download = VEP_config
    data_path = os.path.join(data_path, "variants")
    while (True):
        try:
            files = VEP_receive.recv()
            if files == "finished":
                VEP_receive.close()
                break
            elif files in already_visited or not files:
                continue
            else:
                already_visited.add(files)
                ensembl_id = population_check(files, data_path, download)
                VEP_send.send(ensembl_id)
                #download_VEP_data(variants)
        except Exception as _:
            sleep(90) #to build up the previous processes

    VEP_send.send("finished")
    VEP_send.close()

def potential_hgvs_notations(variant):
    notations = []
    chrom = variant.get("chrom", None)
    transcript_id = variant.get("transcript_id", None)
    gene_id = variant.get("gene_id", None)

    hgvsc = variant.get("hgvsc", None)
    hgvsp= variant.get("hgvsp", None)
    hgvs = variant.get("hgvs", None)

    if hgvsc is not None:
        notations.append(str(chrom) + hgvsc)
        notations.append(str(transcript_id) + hgvsc)
        notations.append(str(gene_id) + hgvsc)

    if hgvsp is not None:
        notations.append(str(chrom) + hgvsp)
        notations.append(str(transcript_id) + hgvsp)
        notations.append(str(gene_id) + hgvsp)

    if hgvs is not None:
        notations.append(str(chrom) + hgvs)
        notations.append(str(transcript_id) + hgvs)
        notations.append(str(gene_id) + hgvs)

    return potential_hgvs_notations


def population_check(files, data_path, download):
    population_check = {}
    gene_dir = os.path.commonpath(files)
    data_path = os.path.join(data_path, os.path.basename(gene_dir))
    os.makedirs(data_path, exist_ok=True)
    for path in files:
        variant_list = open_json(path)
        if not variant_list:
            continue
        for variant in variant_list:
            if not isinstance(variant, dict):
                continue
            variant_id = variant.get("variant_id", None)
            rsids = variant.get("rsids", [])
            if found(variant):
                for rsid in rsids:
                    fetch_rsid_data(data_path, rsid)
                if not rsids:
                    already_run = []
                    for hgvs in potential_hgvs_notations(variant):
                        if hgvs not in already_run:
                            fetch_hgvs_data(data_path, hgvs)
                            already_run.append(hgvs)
            else:
                for rsid in rsids:
                    fetch_pop_data(data_path, rsid)
                if not rsids:
                    already_run = []
                    for hgvs in potential_hgvs_notations(variant):
                        if hgvs not in already_run:
                            rsids = translate_to_rsid(data_path, hgvs)
                            # for rsid in rsids:
                            #     fetch_pop_data(data_path, rsid)
                            already_run.append(hgvs)
    return os.path.basename(gene_dir)

def translate_to_rsid(path_to_data, hgvs):
    p = os.path.join(path_to_data, f"variant_ids_by_{hgvs}.json")
    try:
        server = "https://rest.ensembl.org"
        ext = f"/variant_recoder/human/{hgvs}"
 
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, timeout=120)
        decoded = r.json()
        with open(p, 'w') as file:
            json.dump(r.json(), file)

        print(f"{datetime.now().strftime('%H%M')} ensemble got {hgvs}")
        sleep(0.2)
    except Exception as e:
        print(str(e))

    return []

"""
potentiel data quests for more information
"""
def fetch_hgvs_data(path_to_data, hgvs):
    p = os.path.join(path_to_data, f"variant_by_hgvs_{hgvs}.json")
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
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, timeout=120)

        decoded = r.json()
        with open(p, 'w') as file:
            json.dump(r.json(), file)

        print(f"{datetime.now().strftime('%H%M')} ensemble got {hgvs}")
        sleep(0.2)
        return pd.json_normalize(decoded)

    except Exception as e:
        print(str(e))
        return None

def fetch_rsid_data(path_to_data, rsid):
    p = os.path.join(path_to_data, f"variant_by_rsID_{rsid}.json")
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
        r = requests.get(server+ext, params=options,  headers={ "Content-Type" : "application/json"}, timeout=120)

        # if not r.ok:
        #     print(r.status_code, r.content)

        decoded = r.json()
        with open(p, 'w') as file:
            json.dump(r.json(), file)
    # with open (pth, )json.dump(decoded, pth)
        print(f"{datetime.now().strftime('%H%M')} ensemble got {rsid}")
        sleep(0.2)
        return pd.json_normalize(decoded)
    except Exception as e:
        print(str(e))
        return None


def fetch_pop_data(path_to_data, rsid):
    try:
        server = "https://rest.ensembl.org"
        options = {
            "phenotypes" : 1, 
            "pops" : 1, 
            }
        ext = f"/variation/homo_sapiens/{rsid}?{options}"
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, timeout=120)

        # if not r.ok:
        #     print(r.status_code, r.content)

        decoded = r.json()
        id_dir = os.path.join(path_to_data, rsid)
        os.makedirs(id_dir, exist_ok=True)
        p = os.path.join(id_dir, f"populations_{rsid}.json")

        with open(p, 'w') as file:
            json.dump(r.json(), file)

        print(f"{datetime.now().strftime('%H%M')} ensemble got {rsid}")
        sleep(0.2)
        return pd.json_normalize(decoded)

    except Exception as e:
        print(str(e))
        return None


def found(variant, ancestry=["afr", "nfe"], cutoff=0.05):
    for item in ["joint", "genome", "exome"]:
        populations_an_ac = variant.get(item, None)
        if populations_an_ac is None:
            continue

        ac = populations_an_ac.get("ac", None)
        if ac is None:
            continue

        populations = populations_an_ac.get("populations", [])
        if not populations:
            continue

        for population in populations:
            pop_id = population.get("id", None)
            if pop_id not in ancestry:
                continue

            ac = population.get("ac", None)
            if ac is None:
                continue

            an = population.get("an", None)
            if an is None:
                continue

            try:
                af = float(ac) / float(an)
                # the variant exists in our population with the specific cutoff
                if af >= cutoff:
                    return True
            except Exception as _:
                continue

    return False