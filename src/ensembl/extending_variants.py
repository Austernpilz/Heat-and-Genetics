import requests
import pandas as pd
import json
import queue
import os
import numpy as np
from time import sleep
from datetime import datetime
import threading

from src.ensembl.fetch_ensembl import fetch_hgvs_data, fetch_rsid_data, fetch_pop_data, translate_to_rsid, fetch_from_ensembl
from src.ensembl.getting_closer_to_variants import check_if_exists
from src.helpers.std_out import send_message
from src.helpers.folder_magic import check_string

def lamma(funct, args, threads=[], task_max=1):
    while (len(threads) >= task_max):
            finish_task = threads.pop(0)
            finish_task.join()

    task = threading.Thread(
        target = funct, 
        args = args
        )
    task.start()
    threads.append(task)
    return threads

def open_json(path):
    if isinstance(path, list):
        variant_list = []
        for p in path:
            v = open_json(p)
            if v is not None:
                variant_list.extend(v)
        return variant_list

    if not path.endswith("json") or "gnomAD_variants" not in path:
        return None

    try:
        v = None
        with open(path, 'r') as f:
            v = json.load(f).get("variants", None)
        return v
    except Exception as e:
        send_message(str(e))
        return None

def update_visited(already_visited, processed):
    for item in processed:
        already_visited.add(item)
    return already_visited

def ensemble_now_(data_path, ensembl_id, already_visited, download):
    if ensembl_id in already_visited:
        return already_visited

    if not download and check_if_exists(data_path, ensembl_id):
        send_message(1,1,"ensembl")
        send_message(f"got {ensembl_id}",0,"ensembl")
        already_visited.add(ensembl_id)
        return already_visited

    _ = fetch_from_ensembl(ensembl_id, data_path)
    already_visited.add(ensembl_id)
    return already_visited

def extend_data (VEP_receive, VEP_send, VEP_config): #VEP_send, 
    send_message("starting", 0, "vep")
    already_visited = set()
    data_path, download = VEP_config
    counter = 0
    while (True):
        try:
            if not VEP_receive.empty():
                item = VEP_receive.get()
                if item == "finished":
                    break

                files, populations, found, not_sure, ensembl_id = item
                variant_path = os.path.join(data_path, ensembl_id, "variants")
                os.makedirs(variant_path, exist_ok=True)

                get_variant_data(files, found, variant_path, download)
                found.update(look_up_variant_data(files, not_sure, variant_path, populations, download))
                ensemble_now_(data_path, ensembl_id, already_visited, download)

                VEP_send.put((found, files, variant_path, data_path, ensembl_id))
                counter = 0
            else:
                counter += 1
                sleep(1)

        except Exception as _:
            counter += 1
            sleep(1) #to build up the previous processes

        if counter > 3600 and VEP_receive.empty():
            break

    send_message("finished", 0, "VEP")
    send_message("finished", 0, "ensembl")

def potential_hgvs_notations(variant, first=True):
    hgvsc = variant.get("hgvsc", None)
    hgvsp= variant.get("hgvsp", None)
    hgvs = variant.get("hgvs", None)

    transcript_id = variant.get("transcript_id", None)
    transcript_version = variant.get("transcript_version", None)
    if transcript_id is not None and transcript_version is not None and first:
        if hgvsc is not None:
            return f"{str(transcript_id)}.{str(transcript_version)}:{str(hgvsc)}"
        elif hgvsp is not None:
            return f"{str(transcript_id)}.{str(transcript_version)}:{str(hgvsp)}"
        elif hgvs is not None:
            return f"{str(transcript_id)}.{str(transcript_version)}:{str(hgvs)}"

    if transcript_id is not None and first:
        if hgvsc is not None:
            return f"{str(transcript_id)}:{str(hgvsc)}"
        elif hgvsp is not None:
            return f"{str(transcript_id)}:{str(hgvsp)}"
        elif hgvs is not None:
            return f"{str(transcript_id)}:{str(hgvs)}"

    gene_id = variant.get("gene_id", None)
    if gene_id is not None and first:
        if hgvsc is not None:
            return f"{str(gene_id)}:{str(hgvsc)}"
        elif hgvsp is not None:
            return f"{str(gene_id)}:{str(hgvsp)}"
        elif hgvs is not None:
            return f"{str(gene_id)}:{str(hgvs)}"

    chrom = variant.get("chrom", None)
    if chrom is not None:
        if hgvsc is not None:
            return f"{str(chrom)}:{str(hgvsc)}"
        elif hgvsp is not None:
            return f"{str(chrom)}:{str(hgvsp)}"
        elif hgvs is not None:
            return f"{str(chrom)}:{str(hgvs)}"

    return None

def get_rsids(variant, variant_path, download):
    rsids = variant.get("rsids", [])
    if not rsids:
        hgvs = potential_hgvs_notations(variant)
        translate = translate_to_rsid(variant_path, hgvs, download)
        if not translate:
            hgvs = potential_hgvs_notations(variant, False)
            translate = translate_to_rsid(variant_path, hgvs, download)
        rsids += get_rsids_from_translate(translate)
    return rsids

def get_rsids_from_translate(translate):
    if check_string(translate):
        return []
    try:
        v = None
        with open(translate, 'r') as f:
            v = json.load(f)
        if not v:
            return []
        id_values = []
        for variant in v:
            if not isinstance(variant, dict):
                continue
            for _, name in variant.items():
                if not isinstance(name, dict):
                    continue

                ids = name.get("id", [])
                if isinstance(ids, str):
                    ids = [ids]
                for variant_id in ids:
                    if isinstance(variant_id, str) and variant_id.startswith("rs"):
                        id_values.append(variant_id)

        return id_values
    except Exception as e:
        send_message(f" - coulnd't translate rsid {str(e)}")

    return []

def check_pop(results, pop_ids):
    for r in results:
        v = None
        try:
            with open(r, 'r') as f:
                v = json.load(f)
            if not v:
                continue
            if not isinstance(v, dict):
                continue
            pop_af = {}
            populations = v.get("populations", [])
            #send_message(populations)
            for p in populations:
                if not isinstance(p, dict):
                    continue
                p_id = p.get("population", "NO_ID")
                if not any( [f":{pop_id}" in p_id for pop_id in pop_ids] ):
                    continue

                allele = p.get("allele", None)
                frequency = p.get("frequency", 0)
                if allele is None:
                    allele = "null"
                if allele not in pop_af.keys():
                    pop_af[allele] = [frequency]
                else:
                    pop_af[allele].append(frequency)

            for allele, frequency in pop_af.items():
                if max(frequency) >= 0.05:
                    return True

        except Exception as e:
            send_message(f" - coulnd't check pop {str(e)}")

    return False

def look_up_variant(variant, variant_path, populations, download):
    rsids = get_rsids(variant, variant_path, download)
    results = fetch_pop_data(variant_path, rsids, download)
    return check_pop(results, populations)

def look_up_variant_data(files, not_sure, variant_path, populations, download):
    if not not_sure:
        return
    variant_json = open_json(files)
    if variant_json is None:
        return
    processed = set()
    found = set()
    for variant in variant_json:
        if not isinstance(variant, dict):
            continue
        variant_id = variant.get("variant_id", "NO_ID")
        if variant_id not in not_sure or variant_id in processed:
            continue

        if look_up_variant(variant, variant_path, populations, download):
            found.add(variant_id)
            if get_variant(variant, variant_path, download):
                send_message(f"got {variant_id}", 0, "vep")
            else:
                send_message(f"{variant_id} not found", 0, "vep")
        else:
            send_message(f"{variant_id} not relevant", 0, "vep")
        send_message(1, 1, "vep")
        processed.add(variant_id)
    return found

def get_variant(variant, variant_path, download):
    hgvs = potential_hgvs_notations(variant)
    if fetch_hgvs_data(variant_path, hgvs, download):
        return True

    hgvs = potential_hgvs_notations(variant, False)
    if fetch_hgvs_data(variant_path, hgvs, download):
        return True

    rsids = get_rsids(variant, variant_path, download)
    if fetch_rsid_data(variant_path, rsids, download):
        return True

    return False

def get_variant_data(files, found, variant_path, download):
    if not found:
        return
    variant_json = open_json(files)
    if variant_json is None:
        return
    processed = set()
    for variant in variant_json:
        if not isinstance(variant, dict):
            continue

        variant_id = variant.get("variant_id", "NO_ID")
        if variant_id not in found or variant_id in processed:
            continue

        if get_variant(variant, variant_path, download):
            send_message(f"got {variant_id}", 0, "vep")
        else:
            send_message(f"{variant_id} not found", 0, "vep")
        send_message(1, 1, "vep")
        processed.add(variant_id)
