import requests
import json 
import os
import pandas as pd
from io import StringIO

# 10/s is the rate at which i can shoot requests 

# can use this with search and fetch
# i want hgnc_id and maybe ens
searchableFields = [
    "alias_symbol",
    "symbol", # this is the symbol
    "uniprot_ids", # we have this in amigo
    "ccds_id",
    "locus_type",
    "mgd_id",
    "ena",
    "rna_central_id",
    "location",
    "prev_name",
    "ensembl_gene_id",
    "hgnc_id",
    "status",
    "entrez_id",
    "alias_name",
    "vega_id",
    "mane_select",
    "ucsc_id",
    "refseq_accession",
    "omim_id",
    "rgd_id",
    "locus_group",
    "curator_notes",
    "name", #this we have as well
    "prev_symbol"
]

storedFields = [
    "status",
    "entrez_id",
    "ensembl_gene_id", #we want this
    "prev_name",
    "location",
    "mgd_id",
    "rna_central_id",
    "locus_type",
    "date_name_changed",
    "alias_symbol",
    "symbol",
    "gene_group",
    "prev_symbol",
    "gene_group_id",
    "gtrnadb",
    "name",
    "mamit-trnadb",
    "enzyme_id",
    "orphanet",
    "omim_id",
    "alias_name",
    "snornabase",
    "imgt",
    "date_modified",
    "hgnc_id",
    "merops",
    "mirbase",
    "horde_id",
    "uuid",
    "ena",
    "pseudogene.org",
    "ccds_id",
    "bioparadigms_slc",
    "prev_name",
    "uniprot_ids",
    "cosmic",
    "location",
    "date_symbol_changed",
    "cd",
    "_version_",
    "homeodb",
    "curator_notes",
    "locus_group",
    "pubmed_id",
    "rgd_id",
    "mane_select",
    "refseq_accession",
    "ucsc_id",
    "alias_name",
    "vega_id",
    "lncipedia",
    "date_approved_reserved",
    "lncrnadb",
    "lsdb",
    "agr"
]

def load_simple_search_data(path_to_HGNC):
    path_data = os.path.join(path_to_HGNC, "data")
    df_list = []
    for entry in os.scandir(path_data):
        if not entry.is_file():
            continue
        elif not entry.name.endswith(".csv"):
            continue
        elif not entry.name.startswith("hgnc-symbol-check"):
            continue
        else:
            hgnc_csv = os.path.join(path_data, entry.name)
            try:
                df_list.append(pd.read_csv(hgnc_csv))
            except Exception as e:
                print(str(e))
                print("couldn't load ", hgnc_csv)
                continue

    return pd.concat(df_list)

#col_names_should have tuples of columns that can be used for a look up
# first one is always the Gen Symbol 
# def simple_match_data(list_of_df, hgnc_simple_search, col_names):
#     for df, col in zip(list_of_df, col_names):
#         df["HGNC_ID"] = df[col[0]].map(lambda x : look_up_hgnc"group"].get(x, "NO_GROUP_TERM"))




path_to_HGNC = os.path.join(os.getcwd(), "HGNC")
amigo_HGNC = os.path.join(path_to_HGNC, "data/hgnc-symbol-check_amigo.csv")
disgnet_HNGC = os.path.join(path_to_HGNC, "data/hgnc-symbol-check_disgnet.csv")

df = load_simple_search_data(path_to_HGNC)
sub = df[df["Match type"] == "Approved symbol"]
# print(sub["Approved symbol"].unique().tolist())

r = requests.get("https://rest.genenames.org/fetch/symbol/MED1")

# override encoding by real educated guess as provided by chardet
r.encoding = r.apparent_encoding
# access the data
df = StringIO(r.text)
print(r.text)

# df = pd.read_csv(StringIO(r.json), 
#                          sep="\t", 
#                          dtype=str)

# print(df.head(10))
# hgnc_check_amigo = pd.read_csv(amigo_HGNC)
# hgnc_check_disgnet = pd.read_csv(disgnet_HNGC)
# print(hgnc_check_disgnet[hgnc_check_disgnet["Match type"] != "Approved symbol"])