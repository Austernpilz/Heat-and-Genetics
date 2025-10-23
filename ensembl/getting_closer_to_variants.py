import requests
import pandas as pd
import json
import os
from time import sleep
#http://www.ensembl.org/biomart/martview/ad4dbf2f9ae74dbf5b9cda391d970be9?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id_version|hsapiens_gene_ensembl.default.feature_page.description|hsapiens_gene_ensembl.default.feature_page.start_position|hsapiens_gene_ensembl.default.feature_page.end_position|hsapiens_gene_ensembl.default.feature_page.chromosome_name|hsapiens_gene_ensembl.default.feature_page.hgnc_id|hsapiens_gene_ensembl.default.feature_page.entrezgene_id|hsapiens_gene_ensembl.default.feature_page.uniprot_gn_id&FILTERS=hsapiens_gene_ensembl.default.filters.hgnc_id."HGNC:5970"&VISIBLEPANEL=resultspanel

# import requests, sys
 
def fetch_from_ensembl(id, path_to_ensembl, download=True):
    id_dir = os.path.join(path_to_ensembl, "data", id)
    os.makedirs(id_dir, exist_ok=True)
    p = os.path.join(id_dir, "ensemble_data.json")
    if download:
        try:
            print('download ', id)
            server = "https://rest.ensembl.org"
            ext = f"/lookup/id/{id}?expand=1"
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

            if not r.ok:
                print(r.status_code, r.content)

            decoded = r.json()
            with open(p, 'w') as file:
                json.dump(r.json(), file)
        # with open (pth, )json.dump(decoded, pth)
            print("download succesfull")
            return pd.json_normalize(decoded)
        except Exception as e:
            print(str(e))
    else:
        try:
            return pd.read_json(p)
        except Exception as e:
            print(str(e))
    return None

def get_data(hgnc_df, path_to_ensembl, download=True):
    unique_ensemble_ids = hgnc_df["ensembl_gene_id"].unique().tolist()

    ensemble_data = []
    for id in unique_ensemble_ids:
        sleep(0.1)
        df = fetch_from_ensembl(id, path_to_ensembl, download=download)
        if df is None:
            continue
        elif df.empty:
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