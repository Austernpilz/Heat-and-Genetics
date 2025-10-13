from AmiGo2 import search_and_download as sad
from disgnet import get_tables as dis
from figures import figures as fig

import data_utility as dut

import os
import pandas as pd


this_folder = os.getcwd()
path_to_amigo = os.path.join(this_folder, "AmiGo2")
download_path_amigo = os.path.join(path_to_amigo, "data")
amigo_in_out = os.path.join(path_to_amigo, "include_exclude.txt")

path_to_disgnet = os.path.join(this_folder, "disgnet")
disgnet_in_out = os.path.join(path_to_disgnet, "include_exclude.txt")

#True => Data is downloaded, False => Data needs to be downloaded
df_amigo, df_overview = sad.get_data(download_path_amigo, path_to_amigo, True)
df_disgnet = dis.build_tables(path_to_disgnet)

df_amigo_reduced, df_amigo_plusplus = dut.apply_include_exclude_txt(amigo_in_out, df_amigo, "term")
df_disgnet_reduced, df_disgnet_plusplus = dut.apply_include_exclude_txt(disgnet_in_out, df_disgnet, "disease_name")

rename_dict = {"gene_symbol": "gene", "bioentity_label": "gene", 
                   "group_term" : "term_general", 
                   "term" : "term_specific", "disease_name": "term_specific"}
df_combined_reduced = dut.make_new_table([df_amigo_reduced, df_disgnet_reduced], list(rename_dict.keys()), rename_dict)

fig.sankey_genes_groups(df_amigo_reduced, "bioentity_label", "group_term", "term", gene_cutoff=10, top_genes=100, top_general=20, top_specific=50, name="amigo_sankey_reduced_top_100")
fig.sankey_genes_groups(df_amigo_plusplus, "bioentity_label", "group_term", "term", gene_cutoff=20, top_genes=100, top_general=20, top_specific=50, name="amigo_sankey_plusplus_top_100")
fig.sankey_genes_groups(df_disgnet_reduced, "gene_symbol", "group_term", "disease_name", gene_cutoff=0, top_genes=100, top_general=20, top_specific=30, name="disgnet_sankey_reduced_top_100")
fig.sankey_genes_groups(df_combined_reduced, "gene", "term_general", "term_specific", gene_cutoff=10, top_genes=100, top_general=20, top_specific=50, name="combined_sankey_reduced_top_100")

fig.plot_bipartite_network(df_amigo_reduced, "bioentity_label", "group_term", gene_cutoff=10, max_genes=100, max_groups=30, name="amigo_sankey_reduced_top_100")
fig.plot_bipartite_network(df_amigo_plusplus, "bioentity_label", "group_term", gene_cutoff=20, max_genes=100, max_groups=30, name="amigo_sankey_plusplus_top_100")
fig.plot_bipartite_network(df_disgnet_reduced, "gene_symbol", "group_term",  gene_cutoff=0, max_genes=100, max_groups=30, name="disgnet_sankey_reduced_top_100")
fig.plot_bipartite_network(df_combined_reduced, "gene", "term_general", gene_cutoff=10, max_genes=100, max_groups=30, name="combined_sankey_reduced_top_100")


#path_to_HGNC = os.path.join(this_folder, "HGNC")

# i want to rund every data_set in 3 figures heatmap, sankey_plot, network (maybe 3 times)

# for dataset_to_plot in [df_disgnet_reduced, df_disgnet_plusplus]:
#     term_general = "group_term"
#     term_specific = "disease_name"
#     gen_term = "gene_symbol"

#     fig.sankey_genes_groups(dataset_to_plot, gen_term, term_specific, term_general, gene_cutoff=0)

#     #fig.plot_incidence_heatmap(dataset_to_plot, gen_term, term_general)
#     #fig.plot_incidence_heatmap(dataset_to_plot, gen_term, term_specific)

#     fig.plot_bipartite_network(dataset_to_plot, gen_term, term_general, gene_cutoff=0)
#     fig.plot_bipartite_network(dataset_to_plot, gen_term, term_general, gene_cutoff=0)



#     term_general = "group_term"
#     term_specific = "disease_name"
# print(df_overview.head(10))
# print(df_overview["Name"].unique().tolist())
#print(df_amigo.head(10))

# print(df_disgnet.head(10))
# print(df_disgnet["disease_name"].unique().tolist())



# df_general = make_new_table(df_amigo, df_disgnet)


# sinnvolle übergriffe suchen oder finden (neuer general term)
# alter general term wird specific term
# 
# 4 daten tables
# amigo
# amigo ++
# disgnet
# disgnet ++
# disgnet + amigo 
# disgnet + amigo ++


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


# interessante gene -> HGNC ID approved symbol -> gprofiler -> biomart ensemble  -> gnomAD
# README schritte schreiben

#build_data_table(download_path)



# network zusammen, sankey gen middle, zusammen
# hgnc -> id 
# ensemble -> http://www.ensembl.org/biomart/martview/ad4dbf2f9ae74dbf5b9cda391d970be9?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id_version|hsapiens_gene_ensembl.default.feature_page.description|hsapiens_gene_ensembl.default.feature_page.start_position|hsapiens_gene_ensembl.default.feature_page.end_position|hsapiens_gene_ensembl.default.feature_page.chromosome_name|hsapiens_gene_ensembl.default.feature_page.hgnc_id|hsapiens_gene_ensembl.default.feature_page.entrezgene_id|hsapiens_gene_ensembl.default.feature_page.uniprot_gn_id&FILTERS=hsapiens_gene_ensembl.default.filters.hgnc_id."HGNC:5970"&VISIBLEPANEL=resultspanel
# gnomad chr:start-end -> vcf -> filtern nach allel frequency african//european >= 0.05 
# ad4dbf2f9ae74dbf5b9cda391d970be9