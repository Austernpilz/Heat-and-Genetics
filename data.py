from AmiGo2 import search_and_download as sad
from disgnet import get_tables as dis
from figures import figures as fig
from HGNC import search_and_fetch as hugo
from ensembl import getting_closer_to_variants as ense
from gnomAD import variants as var
import data_utility as dut

import os
import pandas as pd

"""
Paths to the different folders
"""

this_folder = os.getcwd() #this excpects you to start from the Heat_and_Genetics_Folder
path_to_amigo = os.path.join(this_folder, "AmiGo2")
amigo_in_out = os.path.join(path_to_amigo, "include_exclude.txt")

path_to_disgnet = os.path.join(this_folder, "disgnet")
disgnet_in_out = os.path.join(path_to_disgnet, "include_exclude.txt")

path_to_HGNC = os.path.join(this_folder, "HGNC")
path_to_ensemble = os.path.join(this_folder, "ensembl")
path_to_gnomAD = os.path.join(this_folder, "gnomAD")


"""
load basic datasets (amigo, disgnet)
first filtering: reduced and extended dataset
second filtering: compare to hgnc approved symbol (and take that, with id)
build union and intersection
"""

#True => Data is downloaded, False => Data needs to be downloaded
df_amigo, df_overview = sad.get_data(path_to_amigo, False)
df_disgnet = dis.build_tables(path_to_disgnet)

df_amigo_reduced, df_amigo_plusplus = dut.apply_include_exclude_txt(amigo_in_out, df_amigo, "term")
df_disgnet_reduced, df_disgnet_plusplus = dut.apply_include_exclude_txt(disgnet_in_out, df_disgnet, "disease_name")

[df_amigo_reduced, df_disgnet_reduced], rest = hugo.clean_up([df_amigo_reduced, df_disgnet_reduced], path_to_HGNC)
rename_dict = {
    "gene_symbol": "gene", 
    "bioentity_label": "gene", 
    "group_term" : "term_general", 
    "term" : "term_specific", 
    "disease_name": "term_specific"}
print(rest["Input"].unique().tolist())
df_union = dut.make_new_table([df_amigo_plusplus, df_disgnet_plusplus], list(rename_dict.keys())+["HGNC ID"], rename_dict)
#print(len(df_union["gene"].unique()))

unique_genes_amigo = df_amigo_reduced["bioentity_label"].unique()
unique_genes_disgnet = df_disgnet_reduced["gene_symbol"].unique()
df_union = dut.make_new_table([df_amigo_reduced, df_disgnet_reduced], list(rename_dict.keys())+["HGNC ID"], rename_dict)
df_intersection = df_union[ 
    df_union["gene"].isin(unique_genes_amigo) & 
    df_union["gene"].isin(unique_genes_disgnet)
    ]

print(len(df_union["gene"].unique()))
#all of disgnet and the rest from amigo
#205 becuase 5 are weird
top_amigo = df_amigo_reduced["bioentity_label"].value_counts().index[:(210 - len(unique_genes_disgnet))].tolist()
top_200_dataset = df_union[
    df_union["gene"].isin(unique_genes_disgnet) | 
    df_union["gene"].isin(top_amigo)
    ]

print(top_200_dataset[top_200_dataset["gene"].isin(rest["Input"].unique())]["gene"].unique().tolist())

"""
plot data
"""
special = df_intersection["gene"].unique().tolist()
fig.sankey_genes_groups(df_amigo_reduced, "bioentity_label", "group_term", "term", gene_cutoff=5, top_genes=50, top_general=20, top_specific=30, name="amigo_sankey_reduced_top_50")
fig.sankey_genes_groups(df_amigo_plusplus, "bioentity_label", "group_term", "term", gene_cutoff=10, top_genes=50, top_general=20, top_specific=30, name="amigo_sankey_plusplus_top_50")
fig.sankey_genes_groups(df_disgnet_reduced, "gene_symbol", "group_term", "disease_name", gene_cutoff=0, top_genes=50, top_general=20, top_specific=30, name="disgnet_sankey_reduced_top_50")
fig.sankey_genes_groups(top_200_dataset, "gene", "term_general", "term_specific", gene_cutoff=0, top_genes=50, top_general=20, top_specific=30, name="combined_sankey_reduced_top_50")

fig.plot_bipartite_network(df_amigo_reduced, "bioentity_label", "group_term", special_genes=special, gene_cutoff=5, max_genes=50, max_groups=20, name="amigo_network_reduced_top_50")
fig.plot_bipartite_network(df_amigo_plusplus, "bioentity_label", "group_term", special_genes=special, gene_cutoff=10, max_genes=50, max_groups=20, name="amigo_network_plusplus_top_50")
fig.plot_bipartite_network(df_disgnet_reduced, "gene_symbol", "group_term",  special_genes=special, gene_cutoff=0, max_genes=50, max_groups=20, name="disgnet_network_reduced_top_50")
fig.plot_bipartite_network(top_200_dataset, "gene", "term_general", special_genes=special, gene_cutoff=0, max_genes=100, max_groups=20, name="combined_network_reduced_top_50")
fig.plot_bipartite_network(top_200_dataset, "gene", "term_specific", special_genes=special, gene_cutoff=0, max_genes=50, max_groups=30, name="combined_network_reduced_top_50")


"""
load hgnc data
load ensemble data
"""

df_HGNC, rest = hugo.load_HGNC(top_200_dataset, path_to_HGNC, True)
if not rest.empty:
    print("couldn't be loaded:", '\n', rest)

top_200_dataset["gene_id"] = top_200_dataset.apply(
    func= lambda row:  df_HGNC[df_HGNC["symbol"] == row["gene"]]["ensembl_gene_id"],
    axis=1,
    )
save_200 = os.path.join(this_folder, "top_200_genes.tsv")
top_200_dataset.drop_duplicates(ignore_index= True, inplace=True)
top_200_dataset.to_csv(save_200, sep='\t', index=False)
# print(df_HGNC)
df_ensemble = ense.get_data(df_HGNC, path_to_ensemble, True)
egid = df_HGNC["ensembl_gene_id"].unique().tolist()
# gnomad_dict = var.download_data(, path_to_gnomAD)
path_to_data = var.get_data(egid, path_to_gnomAD, ["afr", "nfe"], 0.05, 4,  True)
# path_to_data = os.path.join(path_to_gnomAD, "data")
vcf = os.path.join(this_folder, "test.vcf")
dut.get_vcf(path_to_data, vcf)
save_2000 = os.path.join(this_folder, "top_200_variants.tsv")
cute_df = dut.save_cute_dfs(path_to_data, save_2000)

merged_df = pd.merge(top_200_dataset, cute_df, on="gene_id", how="outer", validate="many_to_many")
save_3000 = os.path.join(this_folder, "merged_6000.tsv")
merged_df.to_csv(save_3000, sep='\t')

# fig.sankey_genes_groups(top_200_dataset, "joint.af_afr", "hgvsc", "gene", gene_cutoff=0, top_genes=50, top_general=20, top_specific=30, name="combined_sankey_reduced_top_50")
# fig.sankey_genes_groups(top_200_dataset, "joint.af_nfe", "hgvsc", "gene", gene_cutoff=0, top_genes=50, top_general=50, top_specific=30, name="combined_sankey_reduced_top_50")
# fig.sankey_genes_groups(top_200_dataset, "hgvsc", "in_silico_predictors.cadd", "joint.af_nfe", gene_cutoff=0, top_genes=50, top_general=50, top_specific=30, name="combined_sankey_reduced_top_50")
# fig.sankey_genes_groups(top_200_dataset, "gene", "in_silico_predictors.cadd", "joint.af_afr", gene_cutoff=0, top_genes=50, top_general=50, top_specific=30, name="combined_sankey_reduced_top_50")
# fig.sankey_genes_groups(top_200_dataset, "gene", "in_silico_predictors.cadd", "joint.af_nfe", gene_cutoff=0, top_genes=50, top_general=50, top_specific=30, name="combined_sankey_reduced_top_50")

# big_dict = var.big_loop(gnomad_dict)
# smaller_dict = var.clean(big_dict)
# smaller_dict.to_csv(os.path.join(path_to_gnomAD, "clean.tsv"), sep='\t')
#print(df_ensemble.head(10))
# print(
# len(gnomad_dict.keys()) == len(df_HGNC),
# len(df_HGNC) == len(top_200_dataset),
# len(gnomad_dict.keys()) == len(top_200_dataset)
# )


