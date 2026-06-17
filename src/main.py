#internal modules
# import argparse
# import os
# import json
from multiprocessing import Pipe
import threading

#global modules
import pandas as pd

#local modules
from src.helpers.folder_magic import get_config


from src.AmiGo2 import search_and_download as amigo
from src.disgnet import get_tables as dn
#from src.helpers import figures as fig
from src.HGNC import search_and_fetch as hugo
from src.ensembl import getting_closer_to_variants as ense
from src.gnomAD import variants as var
import src.helpers.table_magic as dut

def new_task(function, function_arguments, threads, task_max = 1):

    while (len(threads) >= task_max):
        finish_task = threads.pop(0)
        finish_task.join()

    task = threading.Thread(
        target = function, 
        args = function_arguments
        )
    task.start()
    threads.append(task)

    return threads

"""
load parameters
"""
config_file = get_config()#"/Users/m/Desktop/neue_Ablage/bsc_praktikum/tud/Heat_and_Genetics/config/config.json")
data_storage = config_file.get("absolute_file_paths").get("data")
config_storage = config_file.get("absolute_file_paths").get("config")
t = config_file.get("flags").get("threads")

"""
pipeline : DOWNLOAD
"""
amigo_send, hgnc_receive = Pipe()

#Step 1 Download Amigo Data 
amigo_config = amigo.get_config(config_file)
threads = new_task(amigo.download_data, (amigo_send, amigo_config), [], t)

#Step 2 Load Disgnet Data
disgnet_config = dn.get_config(config_file)
disgnet_df = dn.build_tables(disgnet_config)

#Step 3 hgnc
hgnc_send, ensembl_receive = Pipe()
hgnc_config = hugo.get_config(config_file)
threads = new_task(hugo.download_hgnc_data, (hgnc_receive, hgnc_send, disgnet_df, hgnc_config), threads, t)

#Step 4 ensembl
ensembl_send, gnomAD_receive = Pipe()
ensembl_config = ense.get_config(config_file)
threads = new_task(ense.download_data, (ensembl_receive, ensembl_send, ensembl_config), threads, t)

#Step 5 gnomAD
gnomAD_send, VEP_receive = Pipe()
gnomAD_config_download = var.get_config(config_file)
threads = new_task(var.download_data, (gnomAD_receive, gnomAD_send, gnomAD_config_download), threads, t)


#Step 6 collect gnomAD
VEP_send, gnomad_filter_receive = Pipe()
VEP_config = ense.get_config(config_file)
threads = new_task(ense.extend_data, (VEP_receive, VEP_send, VEP_config), threads, t)

#Step 7 filter gnomAD data
gnomAD_config_filter = var.get_config(config_file)
threads = new_task(var.simplify_df, (gnomad_filter_receive, gnomAD_config_filter), threads, t)

"""
pipline: CLEAN DATA
"""
#Step 6 finished work
for t in threads:
    t.join()

dut.save_results(config_file)


# ensemble, data enrichment variant effect predictor
# ensemble, what kind of variant?
# missens, nonsens, .....
# 

# amigo2_data = os.path.join(data_storage, "AmiGo2")
# amigo_df = sad.get_data(amigo2_data, config)
# # """
# # Paths to the different folders
# # """

# # 
# amigo_in_out = os.path.join(path_to_amigo, "include_exclude.txt")

# path_to_disgnet = os.path.join(data_storage, "disgnet")
# disgnet_in_out = os.path.join(path_to_disgnet, "include_exclude.txt")

# # path_to_HGNC = os.path.join(this_folder, "HGNC")
# # path_to_ensemble = os.path.join(this_folder, "ensembl")
# # path_to_gnomAD = os.path.join(this_folder, "gnomAD")


# """
# load basic datasets (amigo, disgnet)
# first filtering: reduced and extended dataset
# second filtering: compare to hgnc approved symbol (and take that, with id)
# build union and intersection
# """

# #False => Data is downloaded, True => Data needs to be downloaded
# df_disgnet = dis.build_tables(path_to_disgnet)

# df_amigo_reduced, df_amigo_plusplus = dut.apply_include_exclude_txt(amigo_in_out, df_amigo, "term")
# df_disgnet_reduced, df_disgnet_plusplus = dut.apply_include_exclude_txt(disgnet_in_out, df_disgnet, "disease_name")

# [df_amigo_reduced, df_disgnet_reduced], rest = hugo.clean_up([df_amigo_plusplus, df_disgnet_plusplus], path_to_HGNC)
# rename_dict = {
#     "gene_symbol": "gene", 
#     "bioentity_label": "gene", 
#     "group_term" : "term_general", 
#     "term" : "term_specific", 
#     "disease_name": "term_specific"}
# # print(rest["Input"].unique().tolist())
# # df_union = dut.make_new_table([df_amigo_plusplus, df_disgnet_plusplus], list(rename_dict.keys())+["HGNC ID"], rename_dict)
# # #print(len(df_union["gene"].unique()))

# unique_genes_amigo = df_amigo_reduced["bioentity_label"].unique()
# unique_genes_disgnet = df_disgnet_reduced["gene_symbol"].unique()
# df_union = dut.make_new_table([df_amigo_reduced, df_disgnet_reduced], list(rename_dict.keys())+["HGNC ID"], rename_dict)
# df_intersection = df_union[ 
#     df_union["gene"].isin(unique_genes_amigo) & 
#     df_union["gene"].isin(unique_genes_disgnet)
#     ]

# # print(len(df_union["gene"].unique()))
# # #all of disgnet and the rest from amigo
# # #205 becuase 5 are weird
# top_amigo = df_amigo_reduced["bioentity_label"].value_counts().index[:(210 - len(unique_genes_disgnet))].tolist()
# top_200_dataset = df_union[
#     df_union["gene"].isin(unique_genes_disgnet) | 
#     df_union["gene"].isin(unique_genes_amigo)
#     ]

# # print(top_200_dataset[top_200_dataset["gene"].isin(rest["Input"].unique())]["gene"].unique().tolist())

# """
# plot data
# """
# # special = df_intersection["gene"].unique().tolist()
# # fig.sankey_genes_groups(df_amigo_reduced, "bioentity_label", "group_term", "term", gene_cutoff=5, top_genes=50, top_general=20, top_specific=30, name="amigo_sankey_reduced_top_50")
# # fig.sankey_genes_groups(df_amigo_plusplus, "bioentity_label", "group_term", "term", gene_cutoff=10, top_genes=50, top_general=20, top_specific=30, name="amigo_sankey_plusplus_top_50")
# # fig.sankey_genes_groups(df_disgnet_reduced, "gene_symbol", "group_term", "disease_name", gene_cutoff=0, top_genes=50, top_general=20, top_specific=30, name="disgnet_sankey_reduced_top_50")
# # fig.sankey_genes_groups(top_200_dataset, "gene", "term_general", "term_specific", gene_cutoff=0, top_genes=50, top_general=20, top_specific=30, name="combined_sankey_reduced_top_50")

# # fig.plot_bipartite_network(df_amigo_reduced, "bioentity_label", "group_term", special_genes=special, gene_cutoff=5, max_genes=50, max_groups=20, name="amigo_network_reduced_top_50")
# # fig.plot_bipartite_network(df_amigo_plusplus, "bioentity_label", "group_term", special_genes=special, gene_cutoff=10, max_genes=50, max_groups=20, name="amigo_network_plusplus_top_50")
# # fig.plot_bipartite_network(df_disgnet_reduced, "gene_symbol", "group_term",  special_genes=special, gene_cutoff=0, max_genes=50, max_groups=20, name="disgnet_network_reduced_top_50")
# # fig.plot_bipartite_network(top_200_dataset, "gene", "term_general", special_genes=special, gene_cutoff=0, max_genes=100, max_groups=20, name="combined_network_reduced_top_50")
# # fig.plot_bipartite_network(top_200_dataset, "gene", "term_specific", special_genes=special, gene_cutoff=0, max_genes=50, max_groups=30, name="combined_network_reduced_top_50")


# """
# load hgnc data
# load ensemble data
# """

# df_HGNC, rest = hugo.load_HGNC(top_200_dataset, path_to_HGNC, False)
# # if not rest.empty:
# #     print("couldn't be loaded:", '\n', rest)

# top_200_dataset["gene_id"] = top_200_dataset.apply(
#     func= lambda row:  df_HGNC[df_HGNC["symbol"] == row["gene"]]["ensembl_gene_id"],
#     axis=1,
#     )
# save_200 = os.path.join(this_folder, "top_2000_genes.tsv")
# top_200_dataset.drop_duplicates(ignore_index= True, inplace=True)
# top_200_dataset.to_csv(save_200, sep='\t', index=False)
# # # print(df_HGNC)
# df_ensemble = ense.get_data(df_HGNC, path_to_ensemble, False)
# egid = df_HGNC["ensembl_gene_id"].unique().tolist()
# path_to_data = var.get_data(egid, path_to_gnomAD, ["afr", "nfe"], 0.05, 4,  True)

# vcf = os.path.join(this_folder, "testtest.vcf")
# dut.build_vcf(path_to_data, vcf)
# save_2000 = os.path.join(this_folder, "top_2000_variants.tsv")
# variants_df = dut.save_variants_to_csv(path_to_data, save_2000)

# merged_df = pd.merge(top_200_dataset, variants_df, on="gene_id", how="outer", validate="many_to_many")
# save_3000 = os.path.join(this_folder, "top2000_genes_w_filtered_variants.tsv")
# merged_df.to_csv(save_3000, sep='\t')
# merged_df = pd.read_csv(save_3000, sep='\t')
# fig.analyze_data(merged_df, False)
# fig.show_corr_top_scores(df_single, "hgvs", "joint.af_afr", "in_silico_predictors.cadd", "in_silico_predictors.cadd", 50)
# fig.show_corr_top_scores(df_single, "hgvs", "joint.af_nfe", "in_silico_predictors.cadd", "in_silico_predictors.cadd", 50)
# fig.show_corr_top_scores(df_single, "hgvs", "joint.af_afr", "in_silico_predictors.cadd", "joint.af_afr", 50)
# fig.show_corr_top_scores(df_single, "hgvs", "joint.af_nfe", "in_silico_predictors.cadd", "joint.af_nfe", 50)


# fig.scatter_allel_score(df_single, "hgvs", "joint.af_afr", "in_silico_predictors.cadd", "in_silico_predictors.cadd", 50)
# fig.scatter_allel_score(df_single, "hgvs", "joint.af_nfe", "in_silico_predictors.cadd", "in_silico_predictors.cadd", 50)
# fig.scatter_allel_score(df_single, "hgvs", "joint.af_afr", "in_silico_predictors.cadd", "joint.af_afr", 50)
# fig.scatter_allel_score(df_single, "hgvs", "joint.af_nfe", "in_silico_predictors.cadd", "joint.af_nfe", 50)














# import numpy as np
# def entropy(row):
#     p = row.values
#     p = p / p.sum()
#     return -(p * np.log2(p + 1e-9)).sum()

# df["pop_entropy"] = df[pop_af_cols].apply(entropy, axis=1)


# Interpretation:

# * **0**: variant almost exclusive to 1 population

# * **2 bits**: variant evenly distributed



# ##  **5. Standardized Z-difference to show deviation from average**

# from scipy.stats import zscore

# df_z = df[pop_af_cols].apply(zscore)
# df[[col + "_z" for col in pop_af_cols]] = df_z


# Variants with high |Z| in one population are outliers.


# # B) Automatically find **variants with strong population differences**

# #  **1. Variants enriched in one population (dominance + threshold)**


# Variant is strongly population-specific if:

# * its highest AF is >=5%
# * it is at least 4x higher than the next highest AF
# * its AF-range is large (> 0.10)


# df["max_af"] = df[pop_af_cols].max(axis=1)
# df["second_af"] = df[pop_af_cols].apply(
#     lambda r: sorted(r)[-2], axis=1
# )

# df["dominance_ratio"] = df["max_af"] / (df["second_af"] + 1e-9)
# strong_pop_specific = df[
#     (df["max_af"] >= 0.05) &
#     (df["dominance_ratio"] >= 4) &
#     (df["af_range"] >= 0.10)
# ]

# This returns variants that are "almost exclusively found in one population".


# #  **2. Variants with low entropy (population-specific)**

# pop_specific_variants = df[df["pop_entropy"] < 0.5]
# Meaning: AF concentrated in one population.


# #  **3. Variants with high entropy (uniform across populations)**

# pop_uniform_variants = df[df["pop_entropy"] > 1.5]

# Meaning: AF consistent across populations.


# # **4. Clustering to detect population-drift patterns**
# ### PCA / UMAP / K-means
# This groups variants with similar population patterns.


# from sklearn.cluster import KMeans

# X = df[pop_af_cols]
# kmeans = KMeans(n_clusters=4, random_state=42).fit(X)

# df["pop_cluster"] = kmeans.labels_


# * cluster vs term_general
# * cluster vs scores
# * cluster vs AF variables

# ##  **Population-dominance barplot**

# sns.countplot(data=df, x="pop_dominant")
# plt.title("Dominant Population per Variant")
# plt.show()

# ##  **Histogram of AF-range**

# sns.histplot(df["af_range"], kde=True)
# plt.title("Distribution of Population AF Differences")
# plt.show()

# ##  **Scatterplot (pop_entropy vs score1)**

# sns.scatterplot(data=df, x="pop_entropy", y="score1", hue="pop_dominant")
# plt.show()


# This often reveals strong patterns:
# * population-specific variants having different score distributions

# ##  **Heatmap of top population-specific variants**


# subset = strong_pop_specific.head(50)

# sns.heatmap(subset.set_index("variant_id")[pop_af_cols], cmap="magma", annot=False)
# plt.show()


# You now have:

# ##  (A) Population-focused metrics:

# 1. dominant population
# 2. AF ratios
# 3. AF range
# 4. entropy
# 5. Z-normalized AF

# ##  (B) Automatic detection code:

# * population-specific variants
# * evenly distributed variants
# * clustering into AF-based groups
# * selection filters with thresholds