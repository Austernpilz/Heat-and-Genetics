
import os
from datetime import datetime
# import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# import plotly.express as px
import plotly.graph_objects as go
# from plotly.subplots import make_subplots as sp
import plotly.io as pio
import networkx as nx
# import itertools as it

from src.helpers.folder_magic import search_for_file, check_string
from src.helpers.config import get_config_results
from src.helpers.std_out import send_message

def plot_first_results(result_path):
    result_path, plots, hgnc_complete, look_up_hugo, amigo_sub, amigo_complete, disgnet_complete, disgnet_sub = get_config_results(result_path)
    load_and_plot_simple(look_up_hugo, amigo_complete, disgnet_complete, hgnc_complete, plots, "complete_dataset")
    load_and_plot_simple(look_up_hugo, amigo_sub, disgnet_sub, hgnc_complete, plots, "reduced_dataset")


def load_and_plot_simple(luh, amigo, disgnet, hgnc, plots, name):
    try:
        amigo_df, disgnet_df = clean_symbols(luh, amigo, disgnet)
        df = build_combined_table(hgnc, amigo_df, disgnet_df)
        file_name = os.path.join(plots, name)
        plot_plot(df, name=file_name)
    except Exception as e:
        send_message(f" - couldn't plot {name}\n{str(e)}\n")


def plot_plot(df, name=""):
    sankey_genes_groups(df, max_genes=10, max_group=10, max_term=10, name=name)
    sankey_genes_groups(df, max_genes=20, max_group=10, max_term=20, name=name)
    sankey_genes_groups(df, max_genes=30, max_group=10, max_term=20, name=name)
    sankey_genes_groups(df, max_genes=40, max_group=20, max_term=30, name=name)
    sankey_genes_groups(df, max_genes=50, max_group=20, max_term=30, name=name)
    sankey_genes_groups(df, max_genes=50, max_group=10, max_term=10, name=name)


    plot_bipartite_network(df, max_genes=20, max_group=10, name=name)
    plot_bipartite_network(df, max_genes=30, max_group=10, name=name)
    plot_bipartite_network(df, max_genes=40, max_group=20, name=name)
    plot_bipartite_network(df, max_genes=50, max_group=20, name=name)
    plot_bipartite_network(df, max_genes=50, max_group=10, name=name)

    plot_bipartite_network(df, max_genes=20, max_group=10, group_col="term", name=name)
    plot_bipartite_network(df, max_genes=30, max_group=10, group_col="term", name=name)
    plot_bipartite_network(df, max_genes=40, max_group=20, group_col="term", name=name)
    plot_bipartite_network(df, max_genes=50, max_group=20, group_col="term", name=name)
    plot_bipartite_network(df, max_genes=50, max_group=10, group_col="term", name=name)


def clean_symbols(look_up_hugo, amigo_path, disgnet_path):
    amigo = pd.read_csv(amigo_path, sep="\t", dtype=str)
    disgnet = pd.read_csv(disgnet_path, sep="\t", dtype=str)
    look_up_hugo = pd.read_csv(look_up_hugo, sep="\t", dtype=str)
    look_up_hugo = look_up_hugo.dropna(subset=["Input", "Approved symbol"])

    luh = {
        input_symbol : gene_symbol
        for input_symbol, gene_symbol in zip(
            look_up_hugo["Input"].tolist(), 
            look_up_hugo["Approved symbol"].tolist()
            )
        if not (check_string(input_symbol) or check_string(gene_symbol))
        }
    amigo["bioentity_label"] = amigo["bioentity_label"].replace(luh)
    disgnet["gene_symbol"] = disgnet["gene_symbol"].replace(luh)
    return amigo, disgnet


def build_combined_table(hgnc_complete, amigo, disgnet):
    hgnc_df_complete = pd.read_csv(hgnc_complete, sep="\t", dtype=str)
    hgnc_part = hgnc_df_complete[["symbol", "name", "ensembl_gene_id", "hgnc_id"]]

    amigo_df = amigo[["bioentity_label", "term", "group"]].rename(columns={"bioentity_label": "symbol"})
    disgnet_df = disgnet[["gene_symbol", "disease_name", "group"]].rename(columns={"gene_symbol" : "symbol", "disease_name": "term"})

    df = pd.concat([amigo_df, disgnet_df], ignore_index=True).drop_duplicates(ignore_index=True)

    return df.merge(hgnc_part, on="symbol", how="left")


# fig is a plotly.graph_objects.Figure or plotly.express figure
def save_plotly(fig, file_name):
    pio.write_html(fig, file=file_name, include_plotlyjs='cdn', auto_open=False)
    send_message(f"Saved: {file_name}", 3)

# # plt.figure() ... plotting code ...
def save_matplotlib(plt, file_name):
    plt.savefig(file_name, dpi=150, bbox_inches="tight")
    plt.close()
    send_message(f"Saved: {file_name}", 3)


# def show_simple_barplot(df, x_name, cut_off, name=False, save=True):
#     print("simple_barplot")
#     data = dut.cut_off(df, cut_off, x_name)
#     plt.figure(figsize=(60,40)) #just to make sure
#     sns.barplot(x=data.index, y=data.values, palette='viridis')
#     plt.xlabel(x_name)
#     plt.ylabel("count")
#     plt.title(f"{x_name} with count > {cut_off}")
#     plt.tight_layout()
#     plt.show()
#     if save:
#         name = name if name else f"barplot_cutoff_{cut_off}"
#         save_matplotlib(plt, name)


# def show_better_barplot(df, x_name, cut_off, name=False, save=True):
#     print("better_barplot")
#     data = dut.cut_off(df, cut_off, x_name)
#     #top_df = cut_off.reset_index()#.rename(columns={'index':'genes','genes':'count'})
#     fig = px.bar(data, x=x_name, y='count', orientation='h',
#                 title=f'genes with count > {cut_off}', height=600)
#     fig.update_layout(yaxis={'categoryorder':'total ascending'})
#     #fig.show()
#     if save:
#         name = name if name else f"barplot_cutoff_{cut_off}"
#         save_plotly(fig, name)


# def plot_incidence_heatmap(df, gene_col, y_name, gene_cutoff=25, top_k_groups=50, name=False, save=True):
#     print("heatmap")
#     # select genes above cutoff and the top_k from another colum
#     data = df[[gene_col, y_name]]
#     genes = dut.cut_off(data, gene_cutoff, gene_col).index.tolist()
#     topky = data[y_name].value_counts().index[:top_k_groups].tolist()

#     data = data[data[gene_col].isin(genes) & data[y_name].isin(topky)]

#     # # pivot to incidence (binary)
#     incidence = pd.crosstab(data[y_name], data[gene_col]).clip(upper=1)
#     # # optionally keep top groups by sum
#     group_sums = incidence.sum(axis=0).sort_values(ascending=False)
#     groups_keep = group_sums.index[:top_k_groups]
#     incidence = incidence[groups_keep]

#     # plot
#     plt.figure(figsize=(60,40))
#     sns.heatmap(incidence, cmap='YlOrBr', cbar=False)
#     plt.xlabel('genes ('+gene_col+')')
#     plt.ylabel('groups ('+y_name+')')
#     plt.title(f'heatmap genes with count > {gene_cutoff} and the top {top_k_groups} most named groups)')
#     plt.tight_layout()
#     plt.show()

#     if save:
#         name = name if name else f"heatmap_cutoff_{gene_cutoff}_topk_{top_k_groups}"
#         save_matplotlib(plt, name)


def plot_bipartite_network(df, max_genes=50, max_groups=10, gene_col="symbol", group_col="group", name="_"):
    send_message("build network graph", 3)

    genes_keep = df_sub[gene_col].value_counts().index[:max_genes].tolist()
    groups_keep = df_sub[group_col].value_counts().index[:max_groups].tolist()
    df_sub = df[df[gene_col].isin(genes_keep) & df[group_col].isin(groups_keep)]

    G = nx.Graph()
    # add nodes with bipartite attribute
    for g in genes_keep:
        G.add_node(('g',g), label=g, bipartite=0)

    for tg in groups_keep:
        G.add_node(('t',tg), label=tg, bipartite=1)

    # add edges
    for _, row in df_sub.iterrows():
        G.add_edge(('g', row[gene_col]), ('t', row[group_col]))

#, width=1, color='pink'
#, width=0.5, color='#888'

    pos = nx.spring_layout(G, k=0.3, seed=42)
    edge_x, edge_y, edge_color, edge_width  = [], [], [], []
    for u,v in G.edges():
        x0,y0 = pos[u]; x1,y1 = pos[v]
        edge_x += [x0, x1]; edge_y += [y0, y1]
        edge_color.append("pink")
        edge_width.append(0.7)

    node_x, node_y, node_text, node_color, node_size = [], [], [], [], []
    for n,data in G.nodes(data=True):
        x,y = pos[n]
        node_x.append(x); node_y.append(y)
        node_text.append(data['label'])
        if data['bipartite']==0:
            node_color.append('blue'); node_size.append(5)
        else:
            node_color.append('red'); node_size.append(10)

    edge_trace = go.Scatter(x=edge_x, y=edge_y, mode='lines', line=dict(width=0.5, color='#888'), hoverinfo='none')
    node_trace = go.Scatter(x=node_x, y=node_y, mode='markers+text', text=node_text,
                            marker=dict(color=node_color, size=node_size), textposition='top center')

    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(showlegend=False, title='Gene-Group network')
    #fig.show()
    time_stamp = datetime.now().strftime('%H%M')
    file_name = f"{time_stamp}_{name}_network_graph_top{max_genes}genes_top{max_groups}groups.html"
    save_plotly(fig, file_name)



def sankey_genes_groups(df, max_genes=30, max_groups=10, max_term=20, gen_column="symbol", group_col="group", term_col="term", name="_"):
    send_message("build sankey groups", 3)

    genes_keep = df_sub[gen_column].value_counts().index[:max_genes].tolist()
    groups_keep = df_sub[group_col].value_counts().index[:max_groups].tolist()
    terms_keep = df_sub[term_col].value_counts().index[:max_term].tolist()
    df_sub = df[df[gen_column].isin(genes_keep) & df[group_col].isin(groups_keep) & df[term_col].isin(terms_keep)]
    df_sub = df[[gen_column, group_col, term_col]]
   

    # build nodes
    nodes = genes_keep + groups_keep + terms_keep
    node_idx = {n:i for i,n in enumerate(nodes)}

    # links terms -> genes
    df_s = df_sub.drop_duplicates([term_col,gen_column]).groupby([term_col,gen_column]).size().reset_index(name='count')

    # links genes -> groups
    df_g = df_sub.drop_duplicates([gen_column,group_col]).groupby([gen_column,group_col]).size().reset_index(name='count')

    source, target, value = [], [], []
    for _,r in df_s.iterrows():
        source.append(node_idx[r[term_col]]); target.append(node_idx[r[gen_column]]); value.append(r['count'])

    for _,r in df_g.iterrows():
        source.append(node_idx[r[gen_column]]); target.append(node_idx[r[group_col]]); value.append(r['count'])

    fig = go.Figure(go.Sankey(node=dict(label=nodes), link=dict(source=source, target=target, value=value)))
    fig.update_layout(title_text="Sankey: specific → genes → general", font_size=20)
    #fig.show()
    time_stamp = datetime.now().strftime('%H%M')
    file_name = f"{time_stamp}_{name}_sankey_top{max_genes}genes_top{max_groups}groups_top{max_term}terms.html"
    save_plotly(fig, file_name)
# usage

# def generate_all(df, amigodf, gen_cut_off=0, top_genes=200, top_groups_general=100, top_groups_specific=100):
#     rename_dict = {"gene_symbol": "gene", "bioentity_label": "gene", 
#                    "group_term" : "term_general", 
#                    "term" : "term_specific", "disease_name": "term_specific"}
#     gen_names = "gene"
#     term_general = "term_general"
#     term_specific = "term_specific"

#     #to make sure we don't accidentally break the table we copy here
#     sub_df = df.rename(mapper=rename_dict, axis=1) #implicit copy

#     # try:
#     #     show_better_barplot(sub_df, gen_names, gen_cut_off)
#     # except Exception as e:
#     #     print(str(e))

#     # try:
#     #     plot_incidence_heatmap(sub_df, gen_names, term_general, gen_cut_off, top_groups_general)
#     # except Exception as e:
#     #     print(str(e))

#     # try:
#     #     plot_incidence_heatmap(sub_df, gen_names, term_specific, gen_cut_off, top_groups_specific)
#     # except Exception as e:
#     #     print(str(e))

#     try:
#         plot_bipartite_network(sub_df, gen_names, term_general, 
#                                gene_cutoff=5, max_genes=top_genes, max_groups=20)
#     except Exception as e:
#         print(str(e))

#     try:
#         plot_bipartite_network(sub_df, gen_names, term_specific, gen_cut_off, top_genes, top_groups_specific)
#     except Exception as e:
#         print(str(e))

#     try:
#         sankey_genes_groups(sub_df, gen_names, term_general, term_specific, gen_cut_off, top_genes, top_groups_general, top_groups_specific)
#     except Exception as e:
#         print(str(e))


# def get_sub_df(df, list_of_columns):
#     # col = list(set(list_of_columns))
#     return df[list_of_columns].drop_duplicates(ignore_index=True)

# def plot_top_10(df, id_col, af_col, score):
#     fig = sp(rows=len(id_col), cols=len(af_col))
#     for i, x in enumerate(id_col):
#         for j, y in enumerate(af_col):
#             data = get_sub_df(df, [x, y, score])
#             if data.empty or data.shape[0]<5:
#                 continue
#             fig.add_trace(
#                 go.Scatter(
#                     x=data[score], y=data[y], hovertext=data[x], 
#                     mode='markers', name=f"{score} vs {y}",
#                 ), 
#                 row=i+1, col=j+1
#             )
#             fig.update_xaxes(title_text=f"{score}", row=i+1, col=j+1)
#             fig.update_yaxes(title_text=f"{y}", row=i+1, col=j+1)
#     fig.update_layout(title=f"{score}")
#     #fig.show()
#     save_plotly(fig, f"simple_top10_{score}")
#     return


# def show_top_10(df, id, af, scores, show=False):
#     df, af = extend(df, af, scores)
#     id_here = ["gene", "hgvs"]
#     if show:
#         sub_df = get_sub_df(df, id_here+af+scores)
#         if sub_df.empty or sub_df.shape[0]<5:
#             return df, id, af, scores

#         for s in scores:
#             data = possible_variants(sub_df, s)
#             if data.empty or data.shape[0]<5:
#                 continue
#             plot_top_10(data, id_here, af, s)

#     return df, id, af, scores

# def simplify(df, show_simple=False):
#     classification_columns = ["gene", "term", "group", "hgvs"]
#     af_columns = ["joint.af_afr", "joint.af_nfe"]
#     score_columns = [
#         "in_silico_predictors.revel_max", "in_silico_predictors.cadd","lof", 
#         "in_silico_predictors.sift_max", "in_silico_predictors.polyphen_max"
#     ]

#     data = get_sub_df(df, classification_columns + af_columns + score_columns)
#     scores = [
#         #scores of deleterious
#         "revel_max",    #0..1 high is bad, >0.5 very likely
#         "cadd",         #0..99 ab 20 eher kritisch ab 30 deleterious
#         "lof",          #HC | LC

#         #score for splice altering
#         #"pangolin_largest_ds",    #0..1 high is splice altering
#         #"spliceai_ds_max"         #0..1 same

#         #score for amino substitution
#         "sift_max",               #0..1 under 0.05 protein function is altered
#         "polyphen_max"            #0..1 high is bad >0.85 dangerous
#     ]
#     af_new = [ af.split(".")[1] for af in af_columns]
#     rename_dicty = {}
#     for score, col in zip(scores, score_columns):
#         rename_dicty[col] = score

#     for af_n, af_o in zip(af_new,af_columns):
#         rename_dicty[af_o] = af_n

#     return show_top_10(
#         data.rename(mapper=rename_dicty, axis=1, errors='ignore'), 
#         classification_columns,
#         af_new,
#         scores, 
#         show_simple
#     )

# def extend(df, af, scores):
#     diff_collum = f"diff_{af[0]}_{af[1]}"
#     df[diff_collum] = np.abs(df[af[0]] - df[af[1]])
#     af += [diff_collum]
#     for s in scores:
#         df[s] = roundabout2(df[s],s)
#     for a in af:
#         df[s].fillna(0.0, inplace=True)
#     return df, af

# def possible_variants(df, score):
#     score_gt = {
#         #scores of deleterious
#         "revel_max" : 0.5,
#         "cadd" : 20,

#         #score for splice altering
#         "pangolin_largest_ds" : 0.8, #0..1 high is splice altering
#         "spliceai_ds_max" : 0.8, #0..1 same

#         "polyphen_max" : 0.85,
#         "sift_max_reverse" : 1 - 0.05
#     }
#     combare_value = score_gt.get(score, False)
#     if combare_value:
#         return df[df[score > combare_value]].copy()
#     else:
#         return df[(df[score] == "HC") | (df[score] == "LC")].copy()



# # def roundabout(arr, s):
# #     arr = roundabout2(arr,s)
# #     match s:
# #         case "revel_max":
# #             return arr.apply(lambda x : np.around(x, 2))
# #         case "cadd":
# #             return arr.apply(lambda x : int(5*(round(float(x)/5))))
# #         case "lof":
# #             return arr
# #         case "sift_max":
# #             return arr.apply(lambda x : np.around(x, 2))
# #         case "polyphen_max":
# #             return arr.apply(lambda x : np.around(x, 2))
# #         case _:
# #             return arr


# # def roundabout2(arr, s):
# #     match s:
# #         case "revel_max":
# #             return arr.fillna(0.0)
# #         case "cadd":
# #             return arr.fillna(0.0)
# #         case "lof":
# #             return arr.fillna("NC")
# #         case "sift_max":
# #             return arr.fillna(1.0)
# #         case "polyphen_max":
# #             return arr.fillna(0.0)
# #         case _:
# #             return arr


#  # Large values = strong population differences.
# def entropy(row):
#     p = row.values
#     p = p / p.sum()
#     return -(p * np.log2(p + 1e-9)).sum()


# def build_sankey_groups(df, id, af, scores):
#     pop = [x for x in af if "diff" in x]
#     id_here = ["gene", "hgvs"]
#     data = get_sub_df(df, id_here+pop+scores).sort_values(by=pop[0],ascending=False, ignore_index=True).head(100)
#     pop = pop[0]
#     var = id_here[1]
#     gene = id_here[0]
#     data[pop] = data[pop].apply(lambda x : np.around(x, 2))
#     largest_diff = data[pop].tolist()
#     largest_group1 = data[var].tolist()
#     largest_group2 = data[gene].tolist()
#     for s in scores:
#         data[s] = roundabout(data[s], s)
#         score_groups = data[s].tolist()
#         largest_diff, largest_group1, largest_group2, score_groups
#         if all(len(x)<5 for x in [largest_diff, largest_group1, largest_group2, score_groups]):
#             continue
#     #data[pop] = data[pop].apply(lambda x : np.around(x, 2))
#         nodes = largest_diff + largest_group1 + largest_group2 + score_groups
#         node_idx = {n:i for i,n in enumerate(nodes)}

#         #links score -> variant
#         df_r = data.drop_duplicates([s,var]).groupby([s,var]).size().reset_index(name='count')

#         # links variant -> diff
#         df_s = data.drop_duplicates([var,pop]).groupby([var,pop]).size().reset_index(name='count')

#         # links diff -> gene
#         df_g = data.drop_duplicates([pop,gene]).groupby([pop,gene]).size().reset_index(name='count')

#         source, target, value = [], [], []


#         for _,r in df_r.iterrows():
#             source.append(node_idx[r[s]])
#             target.append(node_idx[r[var]])
#             value.append(r['count'])

#         for _,r in df_s.iterrows():
#             source.append(node_idx[r[var]])
#             target.append(node_idx[r[pop]])
#             value.append(r['count'])

#         for _,r in df_g.iterrows():
#             source.append(node_idx[r[pop]])
#             target.append(node_idx[r[gene]])
#             value.append(r['count'])

#         fig = go.Figure(go.Sankey(node=dict(label=nodes), link=dict(source=source, target=target, value=value)))
#         fig.update_layout(title_text=f"Sankey: {s} → variant → diff_pop → genes", font_size=20)
#         #fig.show()
#         save_plotly(fig, f"sankey_{s}")
#     return



# def analyze_data(merged_df, show_simple):

#     # differenze netzwerk gewicht
#     # pop differenz im netzwerk
#     # sankey pop dif, score

#     data, id, af, scores = simplify(merged_df, show_simple)

#     data = data[(data["af_afr"]>0.05) | (data["af_nfe"]>0.05)]
#     pop_af_cols = ["af_afr", "af_nfe"]

#     data['pop_dominant'] = data[pop_af_cols].idxmax(axis=1)
#     data["afr_to_nfe_ratio"] = data["af_afr"] / (data["af_nfe"] + 1e-9)
#     data["nfe_to_afr_ratio"] = data["af_nfe"] / (data["af_afr"] + 1e-9)
#     data["af_range"] = data[pop_af_cols].max(axis=1) - data[pop_af_cols].min(axis=1)
#     data["pop_entropy"] = data[pop_af_cols].apply(entropy, axis=1)
#     data["max_af"] = data[pop_af_cols].max(axis=1)
#     data["second_af"] = data[pop_af_cols].apply(
#         lambda r: sorted(r)[-2], axis=1
#     )

#     data["dominance_ratio"] = data["max_af"] / (data["second_af"] + 1e-9)
#     strong_pop_specific = data[
#         (data["max_af"] >= 0.05) &
#         (data["dominance_ratio"] >= 4) &
#         (data["af_range"] >= 0.10)
#     ]



#     pop_specific_variants = data[data["pop_entropy"] < 0.5]
#     pop_uniform_variants = data[data["pop_entropy"] > 1.5]

#     def cluster(data, columns):
#         try:
#             plt.figure(figsize=(60,40))
    
#             sns.clustermap(
#                 data.set_index('hgvs')[columns],
#                 cmap='magma', norm=None
#             )
#             #plt.show()
#             plt.tight_layout()
#             save_matplotlib(plt, f"cluster_{columns}")
#         except Exception as e:
#             print(str(e))
#         return

#     #cluster(data, pop_af_cols)
#     # for s in scores:
#     #     cluster(data, ['pop_dominant',s])

#     # cluster(data, ['pop_dominant', 'pop_entropy'])

#     #ratios
#     # for a in af:
#     #     cluster(data, ['afr_to_nfe_ratio', a])
#     #     cluster(data, ['nfe_to_afr_ratio', a])
#     #Range

#     #cluster(data, ['af_range', 'pop_entropy'])

#     # for s in scores:
#     #     for a in af:
#     #         cluster(data, [a, s])

#     # sub_df = get_sub_df(data, ["hgvs", "gene", 'pop_dominant', "afr_to_nfe_ratio", "nfe_to_afr_ratio", 'af_range', "pop_entropy"]+scores)
#     # for s in scores:
#     #     subdf = devastating_variants(sub_df, s)
#     #     plot_top_10(subdf, ["hgvs", "gene"], ['pop_dominant', "afr_to_nfe_ratio", "nfe_to_afr_ratio", 'af_range', "pop_entropy"], s)
#     #entropy shows population outliers
#     #range,ratio -> inf big diff in population
#     def zzz(data, pop_af_cols):
#         from scipy.stats import zscore
#         try:
#             df_z = data[pop_af_cols].apply(zscore)
#             data[[col + "_z" for col in pop_af_cols]] = df_z
#             return data
#         except Exception as e:
#             return data

#     data = zzz(data, pop_af_cols)
#     # for s in scores:
#     #     plot_top_10(pop_specific_variants, ["hgvs", "gene"], pop_af_cols, s)
#     #     plot_top_10(pop_uniform_variants, ["hgvs", "gene"], pop_af_cols, s)

#     from sklearn.cluster import KMeans

#     X = data[pop_af_cols]
#     kmeans = KMeans(n_clusters=4, random_state=42).fit(X)

#     data["pop_cluster"] = kmeans.labels_
#     #shows variants with similiar pop data

#     # ## ✔️ **Population-dominance barplot**
#     def countplot(data, col):
#         try:
#             plt.figure(figsize=(60,40))
    
#             sns.countplot(data=data, x=col)
#             plt.title(f"col")
#             plt.tight_layout()
#             save_matplotlib(plt, f"countplot_{col}")
#         except Exception as e:
#             print(e)

#     #countplot(data, "pop_cluster")
#     # ## ✔️ **Histogram of AF-range**

#     def hist_plot(data):
#         try:
#             plt.figure(figsize=(60,40))
    
#             sns.histplot(data["af_range"], kde=True)
#             plt.title("Distribution of Population AF Differences")
#             plt.tight_layout()
#             save_matplotlib(plt, f"histplot_af_range")
#         except Exception as e:
#             print(e)

#     #hist_plot(data)
#     # ## ✔️ **Scatterplot (pop_entropy vs score1)**
#     def scatter(data, scores):
#         for s in scores:
#             try:
#                 plt.figure(figsize=(60,40))
    
#                 sns.scatterplot(data=data, x="pop_entropy", y=s, hue="pop_dominant")
#                 plt.title(f"histplot_pop_entropy_{s}")
#                 plt.tight_layout()
#                 save_matplotlib(plt, f"histplot_pop_entropy_{s}")
#             except Exception as e:
#                 print(e)
#     #scatter(data, scores)
#     # This often reveals strong patterns:
#     # * population-specific variants having different score distributions

#     # ## ✔️ **Heatmap of top population-specific variants**

#     def heat(strong_pop_specific, var, pop_af_cols):
#         subset = strong_pop_specific.head(100)
#         for v in var:
#             try:
#                 sns.heatmap(subset.set_index(v)[pop_af_cols], cmap="magma", annot=False)
#                 save_matplotlib(plt, f"heat_map_{v}")
#             except Exception as e:
#                     print(e)

#     #heat(strong_pop_specific, ["gene", "hgvs"], pop_af_cols)
#     #heat(pop_specific_variants, ["gene", "hgvs"], pop_af_cols)
#     #heat(pop_uniform_variants, ["gene", "hgvs"], pop_af_cols)

#     #build_sankey_groups(data, id, af, scores)

#     #build_sankey_groups(strong_pop_specific, ["gene", "hgvs"], af, scores)
#     #build_sankey_groups(pop_specific_variants, ["gene", "hgvs"], af, scores)
#     #build_sankey_groups(pop_uniform_variants, ["gene", "hgvs"], af, scores)

#     #try_networkx2(strong_pop_specific, ["gene", "hgvs"], pop_af_cols, scores)
#     #try_networkx2(pop_specific_variants, ["gene", "hgvs"], pop_af_cols, scores)
#     #try_networkx2(pop_uniform_variants, ["gene", "hgvs"], pop_af_cols, scores)
#     n0 = []
#     n1 = []
#     for score in scores:
#         n0.append(devastating_variants(data, score))
#         n1.append(possible_variants(data, score))
#     # #just_bad_vars = pd.concat(n0,ignore_index=True).drop_duplicates(ignore_index=True)
#     # #maybe_bad_bars = pd.concat(n1,ignore_index=True).drop_duplicates(ignore_index=True)

#     #build_sankey_groups(just_bad_vars, id, af, scores)
#     #build_sankey_groups(maybe_bad_bars, id, af, scores)
#     # for bad_var in n0:
#     #     #build_sankey_groups(bad_var, id, af, scores)
#     #     try_networkx2(bad_var, id, pop_af_cols, scores)
#         #heat(bad_var, ["gene", "hgvs"], pop_af_cols)
#         #scatter(bad_var, scores)
#         #hist_plot(bad_var)
#         #countplot(bad_var, "pop_cluster")

#     # for bad_bar in n1:
#     #     #build_sankey_groups(bad_bar, id, af, scores)
#     #     try_networkx2(bad_bar, id, pop_af_cols, scores)
#     #     #heat(bad_bar, ["gene", "hgvs"], pop_af_cols)
#     #     #scatter(bad_bar, scores)
#     #     #hist_plot(bad_bar)
#     #     #countplot(bad_bar, "pop_cluster")

#     # try_networkx2(data, id, pop_af_cols, scores)
#     return

# def try_networkx2(df, id, af, scores):
#     t = {
#         "HC" : 1,
#         "NC" : 0,
#         "LC" : 0.5
#     }
#     id_here = ["gene", "hgvs"]
#     for s in scores:

#         sub_df = possible_variants(get_sub_df(df, [s]+id_here+af), s)
#         sub_df[s] = roundabout(sub_df[s],s)
#         if s == "lof":
#             sub_df[s] = sub_df[s].apply(lambda x : t[x])
#         af0 = sub_df.nlargest(50, af[0], keep='all')
#         af1 = sub_df.nlargest(50, af[1], keep='all')
#         for i in id_here:
#             ident0 = af0[i].unique().tolist() 
#             ident1 = af1[i].unique().tolist()
#             sub_df_g = sub_df[sub_df[i].isin(ident0+ident1)]
#             if sub_df_g.empty or sub_df_g.shape[0]<5:
#                 continue

#             # G = nx.Graph()

#             # for pop in af:
#             #     G.add_node(('t',pop), label=pop, bipartite=0)

#             # pos_fix = {}
#             # for klm, node in enumerate(G.nodes):
#             #     pos_fix[node] = [(-1)**klm,(-1)**klm]

#             # for g in sub_df_g[i].unique().tolist():
#             #     G.add_node(('g',g), label=g, bipartite=1)

#             # for _, row in sub_df_g.iterrows():
#             #     for pop in af:
#             #         G.add_edge(('g', row[i]), ('t', pop), weight=row[pop])

#             # elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] > 0.3]
#             # esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] <= 0.3]

#             # pos = nx.spring_layout(G, seed=42, pos=pos_fix, fixed=pos_fix.keys())  # positions for all nodes - seed for reproducibility

#             # nx.draw_networkx_nodes(G, pos, node_size=10)
#             # nx.draw_networkx_edges(G, pos, edgelist=elarge, width=4)
#             # nx.draw_networkx_edges(
#             #     G, pos, edgelist=esmall, width=1, alpha=0.5, edge_color="b", style="dashed"
#             # )

#             # nx.draw_networkx_labels(G, pos, font_size=12, font_family="sans-serif")

#             # edge_labels = nx.get_edge_attributes(G, "weight")
#             # nx.draw_networkx_edge_labels(G, pos, edge_labels)

#             # ax = plt.gca()
#             # ax.margins(0.08)
#             # plt.axis("off")
#             # plt.tight_layout()
#             # plt.show()
#             # save_matplotlib(plt, "example_internet")

#             G = nx.Graph()
#             for pop in af:
#                 G.add_node(('t',pop), label=pop, bipartite=0)

#             pos_fix = {}
#             for klm, node in enumerate(G.nodes):
#                 pos_fix[node] = [(-1)**klm*10,(-1)**klm*10]

#             for g in sub_df_g[i].unique().tolist():
#                 G.add_node(('g',g), label=g, bipartite=1)

#             for _, row in sub_df_g.iterrows():
#                 for pop in af:
#                     G.add_edge(('g', row[i]), ('t', pop), weight=row[pop])

#             pos = nx.spring_layout(G, seed=42, pos=pos_fix, fixed=pos_fix.keys())
#             edge_x, edge_y, edge_color, edge_width  = [], [], [], []
#             for u,v in G.edges():
#                 x0,y0 = pos[u]; x1,y1 = pos[v]
#                 edge_x += [x0, x1]; edge_y += [y0, y1]
#                 edge_color.append('#888')
#                 edge_width.append(0.5)

#             node_x, node_y, node_text, node_color, node_size = [], [], [], [], []
#             for n,data in G.nodes(data=True):
#                 x,y = pos[n]
#                 node_x.append(x); node_y.append(y)
#                 node_text.append(data['label'])
#                 if data['bipartite']==0:
#                     node_color.append('blue'); node_size.append(10)
#                 else:
#                     node_color.append('red'); node_size.append(20)

#             edge_trace = go.Scatter(x=edge_x, y=edge_y, mode='lines', line=dict(width=0.5, color='#888'), hoverinfo='none')
#             node_trace = go.Scatter(
#                 x=node_x, y=node_y, mode='markers+text', text=node_text,
#                 marker=dict(color=node_color, size=node_size), textposition='top center'
#             )

#             fig = go.Figure(data=[edge_trace, node_trace])
#             fig.update_layout(showlegend=False, title=f'{i}-network {s}')
#             save_plotly(fig, f"network_graph{i}_{s}_{af}")
#     return

# def try_networkx(df, id, af, scores):
#     af = [a for a in af if "diff" not in a]
#     scores = [s for s in scores if "lof" not in s]
#     features = af+scores
#     X = df[features].values
#     id_here = ["gene", "hgvs"]
#     # correlation between variants
#     corr = np.corrcoef(X)
#     for i in id_here:

#         # create graph
#         G = nx.Graph()

#         # add nodes with labels
#         for i, vid in enumerate(df[i]):
#             G.add_node(i, label=f"{vid}")

#         threshold = 0.75

#         for i in range(len(df)):
#             for j in range(i+1, len(df)):
#                 if corr[i, j] > threshold:
#                     G.add_edge(i, j, weight=corr[i, j])

#         pos = nx.spring_layout(G, k=0.5)

#         plt.figure(figsize=(12, 10))
#         nx.draw_networkx_nodes(G, pos, node_size=30)
#         nx.draw_networkx_edges(G, pos, alpha=0.8)
#         nx.draw_networkx_labels(G, pos, labels=nx.get_node_attributes(G, 'label'), font_size=8)
#         plt.title("Variant Similarity Network (Labeled Variants)")
#         plt.axis('off')
#         save_matplotlib(plt, "graph")

#         edge_x, edge_y, edge_color, edge_width  = [], [], [], []
#         for u,v in G.edges():
#             x0,y0 = pos[u]; x1,y1 = pos[v]
#             edge_x += [x0, x1]; edge_y += [y0, y1]
#             edge_color.append('#888')
#             edge_width.append(0.5)

#         node_x, node_y, node_text, node_color, node_size = [], [], [], [], []
#         for n,data in G.nodes(data=True):
#             x,y = pos[n]
#             node_x.append(x); node_y.append(y)
#             node_text.append(data['label'])
#             node_color.append('pink'); node_size.append(15)

#         edge_trace = go.Scatter(x=edge_x, y=edge_y, mode='lines', line=dict(width=0.5, color='#888'), hoverinfo='none')
#         node_trace = go.Scatter(x=node_x, y=node_y, mode='markers+text', text=node_text, marker=dict(color=node_color, size=node_size), textposition='top center')

#         fig = go.Figure(data=[edge_trace, node_trace])
#         fig.update_layout(showlegend=False, title='Gene-Var bipartite network')
#         save_plotly(fig, "gen_var")
# # def other_clust(df, af):
# #     id = ["hgvs"]
# #     score = ["cadd"]
# #     population = []
# #     for name in af:
# #         if "joint" in name and "diff" not in name:
# #             population.append(name)

# #     sub_df = get_sub_df(df, population+id+score)
# #     sub_df["cadd"] = sub_df["cadd"].apply(cluster_cadd)
# #     print(sub_df)
# #     sub_df = sub_df.nlargest(50, "cadd", keep='all')
# #     corr_matrix = sub_df[population+score].T.corr()

# #     sns.clustermap(corr_matrix, cmap='coolwarm', xticklabels=sub_df[id], yticklabels=sub_df[id])
# #     plt.show()
# #     return



# # filtered by lof/deteriousness, /distance betweeen pop, /highest populations inside