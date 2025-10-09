from AmiGo2 import search_and_download as sad
from disgnet import get_tables as dis

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import networkx as nx
import plotly.graph_objects as go
from sklearn.metrics import jaccard_score
from sklearn.cluster import AgglomerativeClustering
import scipy.cluster.hierarchy as sch
from upsetplot import from_indicators, plot
import time
from datetime import datetime
import plotly.io as pio



this_folder = os.getcwd()
path_to_overview = os.path.join(this_folder, "AmiGo2")
download_path_amigo = os.path.join(path_to_overview, "data")

df_amigo, df_overview = sad.get_data(download_path_amigo, path_to_overview, True)
# print(df_overview.head(10))
print(df_overview["Name"].unique().tolist())
#print(df_amigo.head(10))

path_to_disgnet = os.path.join(this_folder, "disgnet")
df_disgnet = dis.build_tables(path_to_disgnet)
# print(df_disgnet.head(10))
print(df_disgnet["disease_name"].unique().tolist())

def get_col_as_unique_and_count(df, name):
    return df[name].value_counts()

def get_col_as_list(df, name):
    return df[name].tolist()

def get_col_unique(df, name):
    return df[name].unique()

def make_new_table(df_amigo, df_disgnet):
    #this is not efficient for larger datasets but quick to implement
    genes_amigo = get_col_as_list(df_amigo, "bioentity_label")
    genes_disgnet = get_col_as_list(df_disgnet, "gene_symbol")

    term_amigo = get_col_as_list(df_amigo, "annotation_class_label")
    term_amigo_gen = get_col_as_list(df_amigo, "term")
    term_disgnet = get_col_as_list(df_disgnet, "disease_name")

    if (len(genes_amigo) != len(term_amigo) or 
        len(genes_disgnet) != len(term_disgnet)):
        print(len(genes_amigo), len(term_amigo), len(genes_disgnet), len(term_disgnet))
        print("something is off")
        return None
    else:
        return pd.DataFrame(
            {
                "genes" : genes_amigo + genes_disgnet,
                "terme_general" : term_amigo_gen + term_disgnet,
                "terme_specific" : term_amigo + term_disgnet
            }
        )

df_general = make_new_table(df_amigo, df_disgnet)

# print(len(df_general)) #4890 datapoints
# print(get_col_as_unique_and_count(df_general, "genes"))
# print(get_col_as_unique_and_count(df_general, "terme_general"))
# print(get_col_as_unique_and_count(df_general, "terme_specific"))
# print(df_general.count()) #1280 gene

unique_count_genes = get_col_as_unique_and_count(df_general, "genes")
unique_count_terme_general = get_col_as_unique_and_count(df_general, "terme_general")
unique_count_terme_specific = get_col_as_unique_and_count(df_general, "terme_specific")
print(unique_count_genes[unique_count_genes>9])

threshold_genes = unique_count_genes[unique_count_genes>9].sort_values(ascending=False)

def show_simple_barplot(cut_off):
    print("simple_barplot")
    plt.figure(figsize=(20,10))
    sns.barplot(x=cut_off.values, y=cut_off.index, palette='viridis')
    plt.xlabel('count')
    plt.ylabel('gene')
    plt.title('genes with count > 9')
    plt.tight_layout()
    plt.show()




def show_better_barplot(cut_off):
    print("better_barplot")
    top_df = cut_off.reset_index()#.rename(columns={'index':'genes','genes':'count'})
    fig = px.bar(top_df, x='count', y='genes', orientation='h',
                title='genes with count > 9', height=600)
    fig.update_layout(yaxis={'categoryorder':'total ascending'})
    fig.show()



def plot_incidence_heatmap(df, group_col='terme_specific', gene_cutoff=20, top_k_groups=30):
    print("heatmap")
    # select genes above cutoff
    gene_counts = df['genes'].value_counts()
    genes_keep = gene_counts[gene_counts >= gene_cutoff].index.tolist()
    df_sub = df[df['genes'].isin(genes_keep)]

    # pivot to incidence (binary)
    incidence = pd.crosstab(df_sub['genes'], df_sub[group_col]).clip(upper=1)
    # optionally keep top groups by sum
    group_sums = incidence.sum(axis=0).sort_values(ascending=False)
    groups_keep = group_sums.index[:top_k_groups]
    incidence = incidence[groups_keep]

    # plot
    plt.figure(figsize=(12, max(4, len(incidence)*0.2)))
    sns.heatmap(incidence, cmap='YlOrBr', cbar=False)
    plt.xlabel('groups ('+group_col+')')
    plt.ylabel('genes')
    plt.title(f'heatmap genes with count > {gene_cutoff})')
    plt.tight_layout()
    plt.show()





def plot_bipartite_network(df, group_col='terme_specific', gene_cutoff=20, max_genes=100, max_groups=80):
    print("network")
    counts = df['genes'].value_counts()
    genes_keep = counts[counts >= gene_cutoff].index.tolist()[:max_genes]
    df_sub = df[df['genes'].isin(genes_keep)]
    group_counts = df_sub[group_col].value_counts()
    groups_keep = group_counts.index.tolist()[:max_groups]
    df_sub = df_sub[df_sub[group_col].isin(groups_keep)]

    G = nx.Graph()
    # add nodes with bipartite attribute
    for g in genes_keep:
        G.add_node(('g',g), label=g, bipartite=0)
    for tg in groups_keep:
        G.add_node(('t',tg), label=tg, bipartite=1)

    # add edges
    for _, row in df_sub.drop_duplicates(['genes', group_col]).iterrows():
        G.add_edge(('g', row['genes']), ('t', row[group_col]))

    pos = nx.spring_layout(G, k=0.5, seed=42)
    edge_x, edge_y = [], []
    for u,v in G.edges():
        x0,y0 = pos[u]; x1,y1 = pos[v]
        edge_x += [x0, x1, None]; edge_y += [y0, y1, None]

    node_x, node_y, node_text, node_color, node_size = [], [], [], [], []
    for n,data in G.nodes(data=True):
        x,y = pos[n]
        node_x.append(x); node_y.append(y)
        node_text.append(data['label'])
        if data['bipartite']==0:
            node_color.append('blue'); node_size.append(10)
        else:
            node_color.append('red'); node_size.append(6)

    edge_trace = go.Scatter(x=edge_x, y=edge_y, mode='lines', line=dict(width=0.5, color='#888'), hoverinfo='none')
    node_trace = go.Scatter(x=node_x, y=node_y, mode='markers+text', text=node_text,
                            marker=dict(color=node_color, size=node_size), textposition='top center')
    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(showlegend=False, title='Gene-Group bipartite network')
    fig.show()

# usage



def gene_overlap_clustering(df, group_col='terme_specific', gene_cutoff=20, top_n_genes=60):
    print("cluster")
    counts = df['genes'].value_counts()
    genes_keep = counts[counts >= gene_cutoff].index[:top_n_genes].tolist()
    incidence = pd.crosstab(df[df['genes'].isin(genes_keep)]['genes'], df[group_col]).clip(upper=1)
    # compute Jaccard distance matrix
    M = 1 - (incidence.dot(incidence.T) / (incidence.sum(axis=1).values[:,None] + incidence.sum(axis=1).values[None,:] - incidence.dot(incidence.T))).fillna(0)
    # hierarchical clustering dendrogram
    plt.figure(figsize=(10,6))
    sch.dendrogram(sch.linkage(M, method='average'), labels=incidence.index.tolist(), orientation='right')
    plt.title('Hierarchical clustering of genes by shared groups (Jaccard)')
    plt.tight_layout()
    plt.show()

# usage


def gpt_makes_stuff(df_general, group_col):
# choose a few top groups
    print("gpt")
    top_groups = df_general[group_col].value_counts().index[:10].tolist()
    inc = pd.crosstab(df_general['genes'], df_general[group_col]).loc[:, top_groups].clip(upper=1)
    up = from_indicators(inc.columns.tolist(), inc.reset_index(drop=True).values)
    plot(up)
    plt.title('UpSet plot: intersections among top groups')
    plt.show()




def sankey_genes_groups(df, gene_cutoff=30, top_genes=50, top_general=20, top_specific=30):
    print("samley groups")
    gene_counts = df['genes'].value_counts()
    genes = gene_counts[gene_counts >= gene_cutoff].index[:top_genes].tolist()
    df_sub = df[df['genes'].isin(genes)]

    general_top = df_sub['terme_general'].value_counts().index[:top_general].tolist()
    specific_top = df_sub['terme_specific'].value_counts().index[:top_specific].tolist()
    df_sub = df_sub[df_sub['terme_general'].isin(general_top) & df_sub['terme_specific'].isin(specific_top)]

    # build nodes
    nodes = list(genes) + general_top + specific_top
    node_idx = {n:i for i,n in enumerate(nodes)}

    # links genes -> general
    df_g = df_sub.drop_duplicates(['genes','terme_general']).groupby(['genes','terme_general']).size().reset_index(name='count')
    # links general -> specific
    df_s = df_sub.drop_duplicates(['terme_general','terme_specific']).groupby(['terme_general','terme_specific']).size().reset_index(name='count')

    source, target, value = [], [], []
    for _,r in df_g.iterrows():
        source.append(node_idx[r['genes']]); target.append(node_idx[r['terme_general']]); value.append(r['count'])
    for _,r in df_s.iterrows():
        source.append(node_idx[r['terme_general']]); target.append(node_idx[r['terme_specific']]); value.append(r['count'])

    fig = go.Figure(go.Sankey(node=dict(label=nodes), link=dict(source=source, target=target, value=value)))
    fig.update_layout(title_text="Sankey: genes → general → specific", font_size=10)
    fig.show()

# usage


def unique_fname(prefix, ext, outdir="figures"):
    os.makedirs(outdir, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%dT%H%M%SZ")
    return os.path.join(outdir, f"{prefix}_{ts}.{ext}")

# # fig is a plotly.graph_objects.Figure or plotly.express figure
# html_path = unique_fname("barplot_genes_above20", "html")
# pio.write_html(fig, file=html_path, include_plotlyjs='cdn', auto_open=False)
# print("Saved:", html_path)


# # plt.figure() ... plotting code ...
# png_path = unique_fname("incidence_heatmap", "png")
# plt.savefig(png_path, dpi=150, bbox_inches="tight")
# plt.close()
# print("Saved:", png_path)


# # pip install -U kaleido
# png_path = unique_fname("sankey", "png")
# fig.write_image(png_path, width=1400, height=800, scale=2)  # uses kaleido

#show_simple_barplot(threshold_genes)
#show_better_barplot(threshold_genes)
#plot_incidence_heatmap(df_general, group_col='terme_specific', gene_cutoff=20, top_k_groups=40)
#plot_incidence_heatmap(df_general, group_col='terme_general', gene_cutoff=20, top_k_groups=40)
#plot_bipartite_network(df_general, group_col='terme_specific', gene_cutoff=20)
#plot_bipartite_network(df_general, group_col='terme_general', gene_cutoff=20)

#gene_overlap_clustering(df_general, group_col='terme_specific', gene_cutoff=20, top_n_genes=50)
#gene_overlap_clustering(df_general, group_col='terme_general', gene_cutoff=20, top_n_genes=50)

# gpt_makes_stuff(df_general, 'terme_specific')
# gpt_makes_stuff(df_general, 'terme_general')

#sankey_genes_groups(df_general, gene_cutoff=20)


# print(unique_count_terme_general[unique_count_terme_general>9])
# print(unique_count_terme_specific[unique_count_terme_specific>9])
# amigo2 zu disgnet matchen
# gnom ad runterladen komplett auf curie (gen namen)
# gnomad data genetic ancestry group mit laden 
# gnom ad . vcf nach genetic ancestry group vorfiltern 
# alls über 0.05 
# größten unterschiede finden


# interessante gene -> HGNC ID approved symbol -> gprofiler -> biomart ensemble 
# README schritte schreiben

#build_data_table(download_path)
