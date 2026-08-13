from collections import Counter
from io import StringIO
import os
from time import sleep

import requests
import pandas as pd
import random

from src.helpers.folder_magic import check_string
from src.helpers.std_out import send_message



def wait(): #variable ratelimit, because i get blocked sometimes :(
    sleep(random.uniform(5, 35))


# sadly, we need to do it by hand, 
# because we have pseudo headers and requests can't handle them
REQUEST_ARG = {
    "qt" : "standard",
    "indent" : "on",
    "wt" : "csv",
    "rows" : "1000",
    "start" : "0",
    #list of colummns we want #if you give that as list you get an array of the single columns,
    # so just transformed I guess
    "fl" : "", 
    "facet" : "true" , 
    "facet.mincount" : "1" , 
    "facet.sort" : "count" , 
    "json.nl" : "arrarr" ,
    "facet.limit" : "25",   #i think this is responsible for the max number of ch
    "hl" : "true" , 
    "hl.simple.pre" : "%3Cem%22class%22=%22hilite%22%3E" , # '<em"class"="hilite">' : probably some formatting
    "hl.snippets" : "1000" , 
    "csv.encapsulator" : "" , 
    "csv.separator" : "%09" , 
    "csv.header" : "false" , #default is false, if true, only loads the headers
    "csv.mv.separator" : "%2C",
    "fq" : [] , #this is our filter option, can be empty
    "facet.field" : [ #pretty sure these are other possible filters
        "aspect" ,  #something onthology
        "type" , #protein, gene_product, mRNA, ... 
        "evidence_subset_closure_label" ,  #no idea
        "regulates_closure_label" , 
        "isa_partof_closure_label" , 
        "annotation_class_label" , 
        "qualifier" , 
        "annotation_extension_class_closure_label" , 
        "assigned_by" , 
        "panther_family_label" ,
        ],
    "q" : "%2A%3A%2A", #no idea just *:*    <- supposed to look cute, 
    #maybe this also decides on the pseude-header signs
}


list_of_possible_columns = [
    #completly useless
    "id" , 

    #probably good identifiers
    "bioentity" , "bioentity_name" , "bioentity_label", "bioentity_internal_id",

    #taxons
    "taxon" , "taxon_label" , 
    "taxon_closure" , "taxon_closure_label", 
    "taxon_subset_closure" , "taxon_subset_closure_label" ,
    "secondary_taxon" , "secondary_taxon_label" ,
    "secondary_taxon_closure" , "secondary_taxon_closure_label" , 

    "qualifier" , #sometimes protein, sometimes just NOT ???

    #describtions for the data i haven't looked at yet
    "panther_family" , "type" , "reference" , "date",
    "isa_partof_closure_label" , "synonym" , "aspect" , "source" ,
    "panther_family_label" , 
    "has_participant_closure" , "regulates_closure_label" ,
    "has_participant_closure_label" , 
    "regulates_closure" , 
    "isa_partof_closure" , "assigned_by" , 

    "evidence" , "evidence_type_closure" , 
    "evidence_subset_closure_label" , "evidence_label" , "evidence_subset_closure" , 
    "evidence_closure_label" , "evidence_closure" , "evidence_type" , "evidence_with" ,

    "annotation_extension_class_label" , "annotation_extension_class_closure_label" ,
    "annotation_extension_class" , "annotation_class" , "annotation_extension_class_closure" ,
    "annotation_extension_json" , "annotation_class_label" ,
    #mostly empty so far
    "bioentity_isoform", "geospatial_x" , "geospatial_y" , "geospatial_z" , "is_redundant_for"
]

standard = [
    "bioentity" ,
    "bioentity_name" ,
    "qualifier" ,
    "annotation_class" ,
    "annotation_extension_json" ,
    "assigned_by" ,
    "taxon" ,
    "evidence_type" ,
    "evidence_with" ,
    "panther_family" ,
    "type" ,
    "bioentity_isoform" ,
    "reference" ,
    "date"
]

extension_for_this_purpose = [
     #labels are easier for me the rest is codes you need to look up :(
    "bioentity_label",
    "taxon_label",
    #"taxon_subset_closure_label",
    "isa_partof_closure_label" ,
    "regulates_closure_label" , 
    "annotation_class_label" ,
    "annotation_extension_class_label",
    "annotation_extension_class_closure_label",
    "has_participant_closure_label" , #for future reference
    #"panther_family_label"
    ]


def get_col():
    return standard + extension_for_this_purpose


def download_from_amigo2(url, columns, dir_path):
    #name = url[-40:].split("&fq=")[-1]  #looks a bit ugly, but should print the GO_id and then some
    #print(f"downloading GO ID {name}")

    try:
        r = requests.get(url, timeout=180)
        text = r.text
        if r.status_code == 200:
            df = pd.read_csv(StringIO(text), 
                            sep="\t", 
                            dtype=str, 
                            header=None,
                            names=columns)
            return df
        else:
            send_message(f"AmiGo2 returned HTTP {r.status_code}")
            file_path = os.path.join(dir_path, "text.txt")
            os.makedirs(dir_path, exist_ok=True)
            with open(file_path, 'w') as f:
                f.write(text)
            send_message(f" - response saved under: {file_path}\n")

    except Exception as e:
        send_message(f"- couldn't read response from AmiGo2\n{str(e)}\n")
        wait()
    return None


def make_columns_to_string(columns):
    string_list = []
    seperator_encoding = "%2C"
    for item in columns:
        string_list.append(item + seperator_encoding)

    return "".join(string_list)[:-3] #delete the last sign

def get_request_args():
    very_long_string = ""

    global REQUEST_ARG
    request_args = REQUEST_ARG.copy()
    columns = get_col()
    request_args["fl"] = make_columns_to_string(columns)

    for key, value in request_args.items():
        if key == "fq":
            continue
        elif isinstance(value, list):
            for item in value:
                very_long_string += f"&{key}={item}"
        else:
            very_long_string += f"&{key}={value}"

    return very_long_string

def get_filter(go_id):
    filter_fq = [
            "document_category:%22annotation%22", 
            "taxon_subset_closure_label:%22Homo%20sapiens%22",
            f"isa_partof_closure:%22GO%3A{go_id}%22",
            ]
    return "".join( f"&fq={item}" for item in filter_fq )

#this is an iterator function, that gives the urls
def build_full_url_from_go_id(go_id):
    # filter parameters for db request
    # filter_fq = { #these : are pseudo headers and our filters
    #     "genes" : "document_category:%22annotation%22", 
    #     "humans" : "taxon_subset_closure_label:%22Homo%20sapiens%22",
    #     "GO_NR" : lambda go_id : f"isa_partof_closure:%22GO%3A{go_id}%22",
    # } # because of these pseudo headers, normal request module breaks the url
    base_url = "https://golr-aux.geneontology.io/solr/select?defType=edismax&"
    request_args = get_request_args()
    filter_fq = get_filter(go_id)

    return base_url + request_args + filter_fq


def get_term_name(df, dir_path=None):
    if dir_path is not None:
        dir_name, term_name = os.path.split(dir_path)
        if check_string(term_name):
            return term_name

    alternativ_term_name = df["has_participant_closure_label"].unique().tolist()
    if len(alternativ_term_name) == 1 and not check_string(alternativ_term_name[0]):
        term_name = alternativ_term_name[0].strip().replace(' ', '_')
        return term_name

    elif len(alternativ_term_name) == 2 and (check_string(alternativ_term_name[0]) or check_string(alternativ_term_name[1])):
        #here I hope one is empty
        term_name = (alternativ_term_name[0] + alternativ_term_name[1]).strip().replace(' ', '_')
        return term_name

    else:
        complicated_term_name = df["annotation_extension_class_closure_label"].tolist()
        complicated_term_name.sort()
        count_dict = Counter()

        for s in complicated_term_name:
            labels = s.split(',')
            for l in labels:
                count_dict[l.strip()] += 1

        possible_terms = alternativ_term_name
        for k, v in sorted(count_dict.items(), key=lambda item: item[1]):
            if v == len(complicated_term_name):
                possible_terms.append(k)

        for term_name in possible_terms:
            for s in complicated_term_name:
                labels = s.split(',')
                if labels and term_name == labels[-1]:
                    return term_name.strip().replace(' ', '_')
                if len(labels) >= 1:
                    if term_name == labels[-1] or term_name == labels[-2]:
                        return term_name.strip().replace(' ', '_')

    return None


def save_table(df, dir_path):
    os.makedirs(dir_path, exist_ok=True)
    try:
        file_path = os.path.join(dir_path, "data.tsv")
        if os.path.isfile(file_path):
            old_df = pd.read_csv(file_path, sep="\t", dtype=str)
            try:
                df = pd.concat([df, old_df], ignore_index=True).drop_duplicates(ignore_index=True)
            except Exception as e:
                i = 0
                while os.path.isfile(file_path):
                    file_path = os.path.join(dir_path, f"data_{i}.tsv")
                    i+=1
                send_message(f" - table already exists {file_path}\n{str(e)}\n - table saved under {file_path}")

        df.to_csv(file_path, index=False, sep="\t")

        return file_path
    except Exception as e:
        send_message(f"couldn't save table to dir_path {dir_path}\n{str(e)}\n")
        return None

def download_table_from_go_id(data_path, go_id, term_name=None):
    if check_string(term_name):
        term_name = f"GO_ID_{go_id}"
    url = build_full_url_from_go_id(go_id)
    dir_path = os.path.join(data_path, term_name)
    df = download_from_amigo2(url, get_col(), dir_path)
    if df is None:
        return None

    file_path = save_table(df, dir_path)
    wait()
    return file_path
