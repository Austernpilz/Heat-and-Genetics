import requests
import pandas as pd
import json
import os
#http://www.ensembl.org/biomart/martview/ad4dbf2f9ae74dbf5b9cda391d970be9?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id_version|hsapiens_gene_ensembl.default.feature_page.description|hsapiens_gene_ensembl.default.feature_page.start_position|hsapiens_gene_ensembl.default.feature_page.end_position|hsapiens_gene_ensembl.default.feature_page.chromosome_name|hsapiens_gene_ensembl.default.feature_page.hgnc_id|hsapiens_gene_ensembl.default.feature_page.entrezgene_id|hsapiens_gene_ensembl.default.feature_page.uniprot_gn_id&FILTERS=hsapiens_gene_ensembl.default.filters.hgnc_id."HGNC:5970"&VISIBLEPANEL=resultspanel

# import requests, sys
 

server = "https://rest.ensembl.org"
ext = "/lookup/id/ENST00000300651?expand=1"

r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
# if not r.ok:
#   r.raise_for_status()
#   sys.exit()
 
decoded = r.json()
#df = pd.DataFrame.from_dict(decoded, orient='columns', dtype="str")

pth = os.path.join(os.getcwd(), "tet.json")
with open(pth, 'w') as file:
    json.dump(r.json(), file)
# with open (pth, )json.dump(decoded, pth)
d = pd.json_normalize(decoded)
print(d.head(10))
# df = pd.read_json(d)
# print(df.head(10))
print(decoded)
# print(df.head(10))