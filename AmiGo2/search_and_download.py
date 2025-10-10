import requests
import pandas as pd
import os
from time import sleep

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
    "annotation_extension_class_label"
    "annotation_extension_class_closure_label" ,
    "has_participant_closure_label" , #for future reference
    #"panther_family_label"
    ]


def make_columns_to_string(columns):
    string_list = []
    seperator_encoding = "%2C" #tab-symbol
    #"%2C" original, but let's try it like this
    for item in columns:
        string_list.append(item + seperator_encoding)

    return "".join(string_list)[:-3] #delete the last tab sign
    
#print(make_columns_to_string(standard))

#this is an iterator object, that gives the urls
def build_full_url(base_url, param_dict, columns, go_ids):
    # filter parameters for db request
    filter_fq = { #these : are pseudo headers and our filters
        "genes" : "document_category:%22annotation%22", 
        "humans" : "taxon_subset_closure_label:%22Homo%20sapiens%22",
        "GO_NR" : lambda GO : f"isa_partof_closure:%22GO%3A{GO}%22"
    }

    col = make_columns_to_string(columns)
    param_dict["fl"] = col
    
    very_long_string = ""  
    # and_ = "&" # eq_ = "="

    for key, value in param_dict.items():
        if key == "fq":
            continue
        elif isinstance(value, list):
            for item in value:
                very_long_string += f"&{key}={item}"
        else:
            very_long_string += f"&{key}={value}"
    
    for go_id in go_ids:
        filter_list = [filter_fq["genes"], filter_fq["humans"], filter_fq["GO_NR"](go_id)]
        filter_string = "".join( f"&fq={item}" for item in filter_list )
        yield base_url + very_long_string + filter_string


def look_for_overview_GO(path_to_dir):
    if os.path.isfile(path_to_dir):
        if path_to_dir == "overview.txt":
            return path_to_dir
        return False

    if not os.path.isdir(path_to_dir):
        return False

    dir_to_visit = [path_to_dir]
    while dir_to_visit:
        current = dir_to_visit.pop(0)
        try:
            for entry in os.scandir(current):
                if entry.name in {"bin", "data", "include", "lib"}:
                    continue
                elif entry.is_file() and entry.name == "overview.txt":
                    return os.path.join(current, entry.name)
                elif entry.is_dir():
                    dir_to_visit.append(os.path.join(current, entry.name))
        except Exception as _:
            continue
    return False

#print(look_for_overview_GO(os.getcwd())) 

def get_overview(path_to_overview):
    path_to_overview = os.getcwd() if path_to_overview is None else path_to_overview
    path_to_overview = look_for_overview_GO(path_to_overview)
    
    if not path_to_overview:
        raise FileNotFoundError("overview.txt not found")
    overviewtxt = {
            "Accession" : [],
            "Name" : [],
            "Ontology" : [],
            "Synonyms" : [], 
            "Alternate IDs" : [],
            "Definition" : [],
            "not_found" : []
        }
    
    with open(path_to_overview, 'r') as f:
        last_line = ""
        for line in f:
            #print(line)
            if line and line.startswith('#'):
                continue

            elif last_line == "":
                last_line = line.strip()
                continue

            elif last_line in overviewtxt:
                overviewtxt[last_line].append(line.strip())
                last_line = ""

            else:
                print(last_line, line)
                overviewtxt["not_found"].append(line.strip())
                last_line = ""

    # if overviewtxt["not_found"]:
    #     print("couldn't match", overviewtxt["not_found"])

    print(overviewtxt.pop("not_found"))
    #print(overviewtxt)

    norm_accession = [ go_id.replace("GO:", "").strip() for go_id in overviewtxt["Accession"] ]       
    overviewtxt["Accession"] = norm_accession

    return pd.DataFrame.from_dict(overviewtxt)


# print(get_overview(None)) 

def get_go_ids(df_overview):
    return df_overview["Accession"].tolist()

# print(len(get_go_ids(get_overview(None))))


def get_single_table(path_to_download, url, col):
    name = url[-40:].split("&fq=")[-1]  #looks a bit ugly, but should print the GO_id and then some
    print(f"downloading GO ID {name}")
    r = requests.get(url, timeout=60)
    if r.status_code != 200:
        print("failed to download:", r)
        return None
    
    text = r.text
    # print(text)
    # print(url)
    try:
        df = pd.read_csv(StringIO(text), 
                         sep="\t", 
                         dtype=str, 
                         header=None,
                         names=col)
        df["term"] = os.path.basename(path_to_download)
    except Exception as e:
        print(str(e))
        print("couldn't read response")

    os.makedirs(path_to_download, exist_ok=True)
    file_path = os.path.join(path_to_download,"data.tsv")

    try:
        df.to_csv(file_path, index=False, sep="\t")
    except Exception as _:
        with open(file_path, 'w') as f:
            f.write(text)
        return None

    return df


def get_data(download_path = None, path_to_overview = None, is_downloaded = False):
    base_url = "https://golr-aux.geneontology.io/solr/select?defType=edismax&"

    global REQUEST_ARG, standard, extension_for_this_purpose, filter_fq
    columns = standard + extension_for_this_purpose
    overview_df = get_overview(path_to_overview)
    download_path = os.getcwd() if download_path is None else download_path

    if is_downloaded:
        return get_data_from_path(download_path, columns), overview_df

    go_ids = get_go_ids(overview_df)
    list_of_df = []
    for dir_name, url in zip(overview_df["Name"], build_full_url(base_url, REQUEST_ARG, columns, go_ids) ):
        dir_path = os.path.join(download_path, dir_name)
        df = get_single_table(dir_path, url, columns)
        if df is not None:
            print(f"got {dir_name}")
            list_of_df.append(df)

        #to not get blocked, if I use 2 or 3 i will be blocked after 30 tables or so
        sleep(5)

    if list_of_df:
        return pd.concat(list_of_df), overview_df
    else:
        print("apperently it either wasn't downloaded or pandas has some feelings about the data")
        print("let's try and get the data if it was downloaded")
        return get_data_from_path(download_path, columns), overview_df

#print(get_data(os.path.join(os.getcwd(), "try_out_data")))

def get_paths_to_data(path_datatsv = None):
    path_datatsv = os.path.join(os.getcwd()) if path_datatsv is None else path_datatsv

    if not os.path.isdir(path_datatsv):
        return []

    datatsv = []
    dir_to_visit = [path_datatsv]
    while dir_to_visit:
        current = dir_to_visit.pop(0)
        try:
            for entry in os.scandir(current):
                if entry.name in {"bin", "include", "lib", "disgnet", "gnomAD"}:
                    continue
                elif entry.is_dir():
                    dir_to_visit.append(os.path.join(current, entry.name))
                    continue
                elif not entry.is_file():
                    continue

                if entry.name in ["data.tsv", "data.csv"]:
                    datatsv.append(os.path.join(current, entry.name))

        except Exception as _:
            continue
    return datatsv

def get_data_from_path(download_path, columns):
    path_datatsv = get_paths_to_data(download_path)
    df_list = []
    for p in path_datatsv:
        try:
            df_list.append(
                pd.read_csv(p, 
                            sep="\t", 
                            dtype=str))
        except Exception as e:
            print(str(e))

    return pd.concat(df_list)

