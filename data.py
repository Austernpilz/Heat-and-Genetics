from AmiGo2 import search_and_download as sad
from disgnet import get_tables as dis

import os
import pandas as pd


this_folder = os.getcwd()
path_to_overview = os.path.join(this_folder, "AmiGo2")
download_path_amigo = os.path.join(path_to_overview, "data")

df_amigo, df_overview = sad.get_data(download_path_amigo, path_to_overview, True)
# print(df_overview.head(10))
# print(df_overview["Name"].unique().tolist())
#print(df_amigo.head(10))

path_to_disgnet = os.path.join(this_folder, "disgnet")
df_disgnet = dis.build_tables(path_to_disgnet)
# print(df_disgnet.head(10))
print(df_disgnet["disease_name"].unique().tolist())

exlude_amigo = [

]

plusplus_amigo = [
    "cytoplasm_protein_quality_control",
    "brown_fat_cell_differentiation",
    "renal_system_process_involved_in_regulation_of_blood_volume",

    ]

exclude_disgnet = [

]


plusplus_disgnet = [

]



# df_general = make_new_table(df_amigo, df_disgnet)



df_amigo
# print(len(df_general)) #4890 datapoints
# print(get_col_as_unique_and_count(df_general, "genes"))
# print(get_col_as_unique_and_count(df_general, "terme_general"))
# print(get_col_as_unique_and_count(df_general, "terme_specific"))
# print(df_general.count()) #1280 gene

# sinnvolle übergriffe suchen oder finden (neuer general term)
# alter general term wird specific term
# 
# 4 daten tables
# amigo
# amigo ++
# disgnet
# disgnet ++
# disgnet + amigo 
# disgnet + amigo 


# heatmap -> kurz erzeugen
# sankeyplot -> html, #bipartit netzwerk -> html
# table genname + ensamble ID + biomart
# gefilterte daten GO nochmal einzeln
# + 

# unique_count_genes = get_col_as_unique_and_count(df_general, "genes")
# unique_count_terme_general = get_col_as_unique_and_count(df_general, "terme_general")
# unique_count_terme_specific = get_col_as_unique_and_count(df_general, "terme_specific")
# print(unique_count_genes[unique_count_genes>9])


    # genes_amigo = get_col_as_list(df_amigo, "bioentity_label")
    # genes_disgnet = get_col_as_list(df_disgnet, "gene_symbol")

    # term_amigo = get_col_as_list(df_amigo, "annotation_class_label")
    # term_amigo_gen = get_col_as_list(df_amigo, "term")
    # term_disgnet = get_col_as_list(df_disgnet, "disease_name")

    # if (len(genes_amigo) != len(term_amigo) or 
    #     len(genes_disgnet) != len(term_disgnet)):
    #     print(len(genes_amigo), len(term_amigo), len(genes_disgnet), len(term_disgnet))
    #     print("something is off")
    #     return None
    # else:
    #     return pd.DataFrame(
    #         {
    #             "genes" : genes_amigo + genes_disgnet,
    #             "terme_general" : term_amigo_gen + term_disgnet,
    #             "terme_specific" : term_amigo + term_disgnet
    #         }
    #     )



# print(unique_count_terme_general[unique_count_terme_general>9])
# print(unique_count_terme_specific[unique_count_terme_specific>9])
# amigo2 zu disgnet matchen

# fragen für markus sammeln
# gnom ad runterladen komplett auf curie (gen namen)
# gnomad data genetic ancestry group mit laden 
# gnom ad . vcf nach genetic ancestry group vorfiltern 
# alls über 0.05 
# größten unterschiede finden


# interessante gene -> HGNC ID approved symbol -> gprofiler -> biomart ensemble 
# README schritte schreiben

#build_data_table(download_path)
