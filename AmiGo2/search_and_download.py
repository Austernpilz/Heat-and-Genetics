import requests
import pandas as pd

# import json #actually not a json
from io import StringIO

# sadly, we need to do it by hand, 
# because we have pseudo headers and requests can't handle them
REQUEST_ARG = {
    "qt" : "standard",
    "indent" : "on",
    "wt" : "csv",
    "rows" : "100000",
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
    "csv.header" : "true" , #default is false, but headers seem to be an important thing
    "csv.mv.separator" : "%2C",
    "fq" : [] , #this is our filter option, can be empty
    "facet.field" : [ #pretty sure these are other possible filters
        "aspect" , 
        "type" ,
        "evidence_subset_closure_label" , 
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

# exclude auslassen, ist bei default false
# evidence auslassen
# Amigo2 runterziehen
# Disgnet runterziehen
# alles aufbereiten
# amigo2 zu disgnet matchen
# gnom ad runterladen komplett auf curie (gen namen)
# gnomad data genetic ancestry group mit laden 
# gnom ad . vcf nach genetic ancestry group vorfiltern 
# alls über 0.05 
# größten unterschiede finden
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
    "panther_family_label" , "annotation_extension_class_closure_label" ,
    "has_participant_closure" , "regulates_closure_label" ,
    "annotation_extension_class_closure" , "has_participant_closure_label" , 
    "evidence" , "annotation_class_label" , "regulates_closure" , 
    "isa_partof_closure" , "evidence_type_closure" , 
    "evidence_subset_closure_label" , "evidence_label" , 
    "evidence_closure_label" , "evidence_closure" , 
    "annotation_extension_class_label" , "evidence_subset_closure" , 
    "annotation_extension_class" , "annotation_class" , 
    "annotation_extension_json" , "assigned_by" , "evidence_type" , "evidence_with" ,
       
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
# panther 
extension_for_this_purpose = [
     #labels are easier for me the rest is codes you need to look up :(
    "bioentity_label",
    "taxon_label",
    "taxon_subset_closure_label",
    "isa_partof_closure_label" ,
    "regulates_closure_label" , 
    "annotation_class_label" ,
    "annotation_extension_class_closure_label" ,
    "has_participant_closure_label" , #for future reference
    "panther_family_label"
    ]


def make_columns_to_string(columns):
    string_list = []
    comma_encoding = "%09" #"%2C" original, but let's try it like this
    for item in columns:
        string_list.append(item + comma_encoding)

    return "".join(string_list)[:-3] #delete the last tab sign
    

def build_list_of_filters(list_of_GO_ID, filter_fq, genes=True, humans=True): 
    return_list = []
    base_args = []
    if genes:
        base_args.append(filter_fq["genes"])
    if humans:
        base_args.append(filter_fq["humans"])
    for go in list_of_GO_ID:
        fq_argument_list = base_args + [filter_fq["GO_NR"](go)]
        return_list.append(fq_argument_list)
    return return_list #this a list, with the filter arguments

def build_param_string(base_url, dictionary_to_string):
    string_list = []
    # and_ = "&"
    # eq_ = "="
    for key, value in dictionary_to_string.items():
        if isinstance(value, list):
            for item in value:
                string_list.append("&" + key + "=" + item)
        else:
            string_list.append("&" + key + "=" + value)
    return base_url.join(return_str)
            

def get_list_of_full_url(url, request_params, list_of_filter):
    return_list = []
    for filter_q in list_of_filter:
        request_params["fq"] = filter_q
        return_list.append(build_param_string(url, request_params))
    return return_list

def download_data_to(path_to_download, url_list, list_of_GO_IDs):
    os.makedirs(path_to_download, exist_ok=True)
    for url, go_id in zip(url_list, list_of_GO_IDs):
        new_dir = os.path.join(path_to_download, go_id)
        os.makedirs(new_dir, exist_ok=True)

        r = request.get(url)
        if r.status_code == 200:
            file_path = os.path.join(new_dir, "data.tsv") #should become a tsv because of the tab, but let's see
            with open(file_path, 'wb') as f:
                f.write(r.content)
        else:
            print(go_id, " couldn't be downloaded, the folder will be empty")



def down_load_all(BASE_URL_DB, REQUEST_ARG, filter_fq, list_of_GO_IDs, columns) 
    columns_to_get = make_columns_to_string(columns)
    REQUEST_ARG["fl"] = columns_to_get

    list_of_filter = build_list_of_filters(list_of_GO_IDs, filter_fq)

    list_of_full_url = get_list_of_full_url(BASE_URL_DB, REQUEST_ARG, list_of_filter)
    this_folder = os.getcwd()

    download_data_to(this_folder, list_of_full_url, list_of_GO_IDs)


filter_fq = { #these : are pseudo headers
    "genes" : "document_category:%22annotation%22", 
    "humans" : "taxon_subset_closure_label:%22Homo%20sapiens%22",
    "GO_NR" : lambda GO : "isa_partof_closure:%22GO%3A{GO}%22"
}
list_of_GO_IDs = ["0016048"]

# standard + extension_for_this_purpose

# https://golr-aux.geneontology.io/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=bioentity%2Cbioentity_name%2Cqualifier%2Cannotation_class%2Cannotation_extension_json%2Cassigned_by%2Ctaxon%2Cevidence_type%2Cevidence_with%2Cpanther_family%2Ctype%2Cbioentity_isoform%2Creference%2Cdate%2Cbioentity_label%2Ctaxon_label%2Ctaxon_subset_closure_label%2Cisa_partof_closure_label%2Cregulates_closure_label%2Cannotation_class_label%2Cannotation_extension_class_closure_label%2Chas_participant_closure_label%2Cpanther_family_label&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem+%22class%22+%3D+%22hilite%22+%3E+&hl.snippets=1000&csv.encapsulator=&csv.separator=+&csv.header=true&csv.mv.separator=%2C&fq=document_category%3A%22annotation%22&fq=taxon_subset_closure_label%3A%22Homo+sapiens%22&fq=isa_partof_closure%3A%22GO%3A%7BGO%7D%22&facet.field=aspect&facet.field=type&facet.field=evidence_subset_closure_label&facet.field=regulates_closure_label&facet.field=isa_partof_closure_label&facet.field=annotation_class_label&facet.field=qualifier&facet.field=annotation_extension_class_closure_label&facet.field=assigned_by&facet.field=panther_family_label&q=%2A%3A%2A