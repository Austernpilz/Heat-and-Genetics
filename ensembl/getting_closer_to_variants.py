import requests
import pandas as pd
import json
import os
from time import sleep
#http://www.ensembl.org/biomart/martview/ad4dbf2f9ae74dbf5b9cda391d970be9?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id_version|hsapiens_gene_ensembl.default.feature_page.description|hsapiens_gene_ensembl.default.feature_page.start_position|hsapiens_gene_ensembl.default.feature_page.end_position|hsapiens_gene_ensembl.default.feature_page.chromosome_name|hsapiens_gene_ensembl.default.feature_page.hgnc_id|hsapiens_gene_ensembl.default.feature_page.entrezgene_id|hsapiens_gene_ensembl.default.feature_page.uniprot_gn_id&FILTERS=hsapiens_gene_ensembl.default.filters.hgnc_id."HGNC:5970"&VISIBLEPANEL=resultspanel

# import requests, sys
 
def fetch_from_ensembl(id, path_to_ensembl):
    id_dir = os.path.join(path_to_ensembl, "data", id)
    os.makedirs(id_dir, exist_ok=True)
    p = os.path.join(id_dir, "ensemble_data.json")
    try:
        print('download ', id)
        server = "https://rest.ensembl.org"
        ext = f"/lookup/id/{id}?expand=1"
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

        # if not r.ok:
        #     print(r.status_code, r.content)

        decoded = r.json()
        with open(p, 'w') as file:
            json.dump(r.json(), file)
    # with open (pth, )json.dump(decoded, pth)
        print("download succesfull")
        sleep(0.1)
        return pd.json_normalize(decoded)
    except Exception as e:
        print(str(e))
        return None

def get_from_path(unique_ensemble_ids, path_to_ensembl):
    ensemble_data = []
    for id in unique_ensemble_ids:
        path_to_id = os.path.join(path_to_ensembl, "data", id, "ensemble_data.json")
        if os.path.isfile(path_to_id):
            try:
                ensemble_data.append(pd.read_json(path_to_id))
                continue
            except Exception as e:
                print(str(e))

        df = fetch_from_ensembl(id, path_to_ensembl)
        if df is None:
            continue
        else:
            ensemble_data.append(df)

    return pd.concat(ensemble_data)


def get_data(hgnc_df, path_to_ensembl, download=True):
    unique_ensemble_ids = hgnc_df["ensembl_gene_id"].unique().tolist()

    if not download:
        return get_from_path(unique_ensemble_ids, path_to_ensembl)

    ensemble_data = []
    for id in unique_ensemble_ids:
        df = fetch_from_ensembl(id, path_to_ensembl)
        if df is None:
            continue
        else:
            ensemble_data.append(df)

    return pd.concat(ensemble_data)
#df = pd.DataFrame.from_dict(decoded, orient='columns', dtype="str")

# print(d.head(10))
# # df = pd.read_json(d)
# # print(df.head(10))
# print(decoded)
# print(df.head(10))

def fetch_hgvs_data(path_to_data, hgvs):
    p = os.path.join(path_to_data, f"variant_by_hgvs_{hgvs}.json")
    try:
        print('loading_variation_data ', id)
        server = "https://rest.ensembl.org"
        options = "Conservation, CADD, AlphaMissense, "
        ext = f"/vep/homo_sapiens/hgvs/{hgvs}?{options}"
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

        # if not r.ok:
        #     print(r.status_code, r.content)

        decoded = r.json()
        with open(p, 'w') as file:
            json.dump(r.json(), file)
    # with open (pth, )json.dump(decoded, pth)
        print("download succesfull")
        sleep(0.1)
        return pd.json_normalize(decoded)

    except Exception as e:
        print(str(e))
        return pd.DataFrame()

def fetch_rsid_data(path_to_data, rsid):
    p = os.path.join(path_to_data, f"variant_by_rsID_{rsid}.json")
    try:
        print('loading_variation_data ', id)
        server = "https://rest.ensembl.org"
        options = ""
        ext = f"/vep/homo_sapiens/id/{rsid}?{options}"
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

        # if not r.ok:
        #     print(r.status_code, r.content)

        decoded = r.json()
        with open(p, 'w') as file:
            json.dump(r.json(), file)
    # with open (pth, )json.dump(decoded, pth)
        print("download succesfull")
        sleep(0.1)
        return pd.json_normalize(decoded)
    except Exception as e:
        print(str(e))
        return pd.DataFrame()
