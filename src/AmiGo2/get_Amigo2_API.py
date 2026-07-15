from io import StringIO
import os
from time import sleep

import requests
import pandas as pd

SEED = 42

def wait(): #variable ratelimit, because i get blocked sometimes :(
    global SEED
    SEED ^= SEED << 13
    SEED ^= SEED >> 7
    SEED ^= SEED << 17
    sleep(6 + SEED % 6)

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
    "secondary_taxon" , "secondary_taxon_label" 
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
    "annotation_extension_json" , "annotation_class_label" 
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
    global standard, extension_for_this_purpose
    return standard + extension_for_this_purpose



def download_from_amigo2(url, columns, dir_path):
    #name = url[-40:].split("&fq=")[-1]  #looks a bit ugly, but should print the GO_id and then some
    #print(f"downloading GO ID {name}")
    r = requests.get(url, timeout=120)
    if r.status_code == 200:
        text = r.text
        try:
            return pd.read_csv(StringIO(text), 
                            sep="\t", 
                            dtype=str, 
                            header=None,
                            names=columns)
        except Exception as e:
            print(str(e))
            print("couldn't read response")
            file_path = os.path.join(dir_path, "text.txt")
            os.makedirs(dir_path, exist_ok=True)
            with open(file_path, 'w') as f:
                f.write(text)
            print("response saved under: ", file_path)

    print("failed to download from:", url)
    return None

def make_columns_to_string(columns):
    string_list = []
    seperator_encoding = "%2C"
    for item in columns:
        string_list.append(item + seperator_encoding)

    return "".join(string_list)[:-3] #delete the last sign


#this is an iterator function, that gives the urls
def build_full_url_from_go_id(go_ids):
    base_url = "https://golr-aux.geneontology.io/solr/select?defType=edismax&"

    columns = get_col()
    # filter parameters for db request
    filter_fq = { #these : are pseudo headers and our filters
        "genes" : "document_category:%22annotation%22", 
        "humans" : "taxon_subset_closure_label:%22Homo%20sapiens%22",
        "GO_NR" : lambda GO : f"isa_partof_closure:%22GO%3A{GO}%22"
    } # because of these pseudo headers, normal request module breaks the url

    col = make_columns_to_string(columns)
    global REQUEST_ARG
    request_args = REQUEST_ARG.copy()
    request_args["fl"] = col

    very_long_string = ""
    # and_ = "&" # eq_ = "="

    for key, value in REQUEST_ARG.items():
        if key == "fq":
            continue
        elif isinstance(value, list):
            for item in value:
                very_long_string += f"&{key}={item}"
        else:
            very_long_string += f"&{key}={value}"

    while go_ids:
        go_id = go_ids.pop()
        filter_list = [filter_fq["genes"], filter_fq["humans"], filter_fq["GO_NR"](go_id)]
        filter_string = "".join( f"&fq={item}" for item in filter_list )
        yield base_url + very_long_string + filter_string
