import requests
import pandas as pd
import json
import os
import numpy as np
from time import sleep
from datetime import datetime
import threading

from src.ensembl.fetch_ensembl import fetch_hgvs_data, fetch_rsid_data, fetch_pop_data, translate_to_rsid

def lamma(funct, args, threads=[]):
    task = threading.Thread(
        target = funct, 
        args = args
        )
    task.start()
    threads.append(task)
    return threads

def open_json(path):
    if not path.endswith("json") or "gnomAD_variants" not in path:
        return None

    try:
        v = None
        with open(path, 'r') as f:
            v = json.load(f).get("variants", None)
        return v
    except Exception as e:
        print(str(e))
        return None

def extend_data (VEP_receive, VEP_config): #VEP_send, 
    print("extending variants")
    files = []
    already_visited = set()
    data_path, download = VEP_config
    counter = 0
    while (True):
        try:
            if VEP_receive.poll(timeout=600):
                files = VEP_receive.recv()
            else:
                counter += 1
                files = "NO ID"

            if files == "finished":
                VEP_receive.close()
                break
            elif not files or files == "NO ID":
                if counter > 10:
                    break
                continue
            gene_dir = os.path.commonpath(files)
            ensembl_id = os.path.basename(gene_dir)
            if ensembl_id in already_visited:
                continue
            else:
                variant_path = os.path.join(data_path, ensembl_id, "variants")
                os.makedirs(variant_path, exist_ok=True)
                get_variant_data(files, variant_path, download)
                #VEP_send.send(variant_path)
                already_visited.add(ensembl_id)

        except Exception as _:
            sleep(180) #to build up the previous processes

    #VEP_send.send("finished")
    #VEP_send.close()
    print("VEP thread done")

def potential_hgvs_notations(variant):
    hgvsc = variant.get("hgvsc", None)
    hgvsp= variant.get("hgvsp", None)
    hgvs = variant.get("hgvs", None)

    transcript_id = variant.get("transcript_id", None)
    if transcript_id is not None:
        if hgvsc is not None:
            return f"{str(transcript_id)}:{str(hgvsc)}"
        elif hgvsp is not None:
            return f"{str(transcript_id)}:{str(hgvsp)}"
        elif hgvs is not None:
            return f"{str(transcript_id)}:{str(hgvs)}"

    gene_id = variant.get("gene_id", None)
    if gene_id is not None:
        if hgvsc is not None:
            return f"{str(transcript_id)}:{str(hgvsc)}"
        elif hgvsp is not None:
            return f"{str(transcript_id)}:{str(hgvsp)}"
        elif hgvs is not None:
            return f"{str(transcript_id)}:{str(hgvs)}"

    chrom = variant.get("chrom", None)
    if chrom is not None:
        if hgvsc is not None:
            return f"{str(transcript_id)}:{str(hgvsc)}"
        elif hgvsp is not None:
            return f"{str(transcript_id)}:{str(hgvsp)}"
        elif hgvs is not None:
            return f"{str(transcript_id)}:{str(hgvs)}"

    return None


def get_variant_data(files, data_path, download):
    for pathpath in files:
        variant_json = open_json(pathpath)
        if variant_json is None:
            continue
        for variant in variant_json:
            if not isinstance(variant, dict):
                continue
            variant_id = variant.get("variant_id", "NO ID")
            variant_path = os.path.join(data_path, variant_id)
            os.makedirs(variant_path, exist_ok=True)

            hgvs = potential_hgvs_notations(variant)
            rsids = variant.get("rsids", []) + translate_to_rsid(data_path, hgvs)
            unique = list(set(rsids))

            threads = lamma(fetch_hgvs_data, (data_path, hgvs, download))

            for rsid in unique:
                threads = lamma(fetch_pop_data, (data_path, rsid, download))
                threads = lamma(fetch_rsid_data, (data_path, rsid, download))

            for t in threads:
                t.join()