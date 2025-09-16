import requests


#I think these symbol mean ...

#next flag
AND_ = "&"
column_specifier = "%2C"

# there is an sql database behind, so we can filter, which makes these urls shorter
get_from_database = "https://golr-aux.geneontology.io/solr/select?defType=edismax"

#seems to be always the first
first_args = [
    "qt=standard",
    "indent=on",
    "wt=csv",
    "rows=100000",
    "start=0",
]

#this gives us arguments for the Columnswe want to have in our csv
#at the end we specify how many char per columns, so not all columns are usefull (25)
second_string = "fl="

#probably something like convert too (2) column (C)


list_of_columns = [
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
    "panther_family_label" , 
]

# think all of them need to be sended, not sure. 
# # but they usually come in this order
# probably qualifiers for the document
third_string = [
    "facet=true" , "facet.mincount=1" , "facet.sort=count" , 
    "json.nl=arrarr" , "facet.limit=25", #this seems to be the limit of shown columns
    "hl=true" , "hl.simple.pre=%3Cem%20class%3D%22hilite%22%3E" , 
    "hl.snippets=1000" , "csv.encapsulator=" , 
    "csv.separator=%09" , "csv.header=true" , #default is false, but headers seem to be an important thing
    "csv.mv.separator=%7C"
]

# last string here we can sort
#fq = filter for
#facet.field = don't filter for??!?
fourth_string = lambda GO : [
            "fq=document_category:%22annotation%22" , #this is for genes, we could in theory look for other specifiers
            f"fq=isa_partof_closure:%22GO%3A{GO}%22" ,
            "facet.field=aspect" ,
            f"fq=taxon_subset_closure_label:%22Homo%20sapiens%22" , #here we can filter for HomoSapiens
            "facet.field=type" ,
            "facet.field=evidence_subset_closure_label" ,
            "facet.field=regulates_closure_label" ,
            "facet.field=isa_partof_closure_label" ,
            "facet.field=annotation_class_label" ,
            "facet.field=qualifier" ,
            "facet.field=annotation_extension_class_closure_label" ,
            "facet.field=assigned_by" ,
            "facet.field=panther_family_label" ,
            "q=*%3A*"
        ]


def return_list_as_string(l):
    return_str = ""
    for item in l:
        return_str+= AND_ + item
    return return_str

def build_all_the_strings_for_get(list_of_GO_ID): 
    return_str = []
    for item in list_of_GO_ID:    
        return_str.append(return_list_as_string(fourth_string(item)))
    return return_str



first_ = return_list_as_string(first_args)
print(first_)
base_url = get_from_database + first_

def get_columns_as_string(xs_col):
    return_col = "&fl="
    for item in xs_col:
        return_col += item + column_specifier
        #the first one seem to not need and, and the last one needs to be deleted.
    return return_col[:-3]


columns = get_columns_as_string(standard + extension_for_this_purpose)
first_half = base_url + columns
# print(first_half)

# GO_ID_test = ["1990845"]

# second_half = build_all_the_strings_for_get(GO_ID_test)

def give_full_url(url, list_of_GO_ID):
    return_list = []
    for seconds in build_all_the_strings_for_get(list_of_GO_ID):
        return_list.append(url + seconds)
    return return_list

# print(give_full_url(first_half, GO_ID_test))


for url in give_full_url(first_half, GO_ID_test):
    #print(url)
    r = requests.get(url)
    #print(r.text)

# https://golr-aux.geneontology.io/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=bioentity%2Cbioentity_name%2Cqualifier%2Cannotation_class%2Cannotation_extension_json%2Cassigned_by%2Ctaxon%2Cevidence_type%2Cevidence_with%2Cpanther_family%2Ctype%2Cbioentity_isoform%2Creference%2Cdate&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class%3D%22hilite%22%3E&hl.snippets=1000&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=%7C&fq=document_category:%22annotation%22&fq=isa_partof_closure:%22GO%3A0009408%22&fq=taxon_subset_closure_label:%22Homo%20sapiens%22&facet.field=aspect&facet.field=taxon_subset_closure_label&facet.field=type&facet.field=evidence_subset_closure_label&facet.field=regulates_closure_label&facet.field=isa_partof_closure_label&facet.field=annotation_class_label&facet.field=qualifier&facet.field=annotation_extension_class_closure_label&facet.field=assigned_by
# https://golr-aux.geneontology.io/solr/select?defType=edismax&qt=standardindent=onwt=csvrows=100000start=0&fl=bioentity%2Cbioentity_name%2Cqualifier%2Cannotation_class%2Cannotation_extension_json%2Cassigned_by%2Ctaxon%2Cevidence_type%2Cevidence_with%2Cpanther_family%2Ctype%2Cbioentity_isoform%2Creference%2Cdate%2Cbioentity_label%2Ctaxon_label%2Ctaxon_subset_closure_label%2Cisa_partof_closure_label%2Cregulates_closure_label%2Cannotation_class_label%2Cannotation_extension_class_closure_label%2Chas_participant_closure_label%2Cpanther_family_label&fq=document_category:%22annotation%22&fq=isa_partof_closure:%22GO%3A1990845%22&facet.field=aspect&fq=taxon_subset_closure_label:%22Homo%20sapiens%22&facet.field=type&facet.field=evidence_subset_closure_label&facet.field=regulates_closure_label&facet.field=isa_partof_closure_label&facet.field=annotation_class_label&facet.field=qualifier&facet.field=annotation_extension_class_closure_label&facet.field=assigned_by