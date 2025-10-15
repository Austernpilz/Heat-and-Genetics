
import os
from datetime import datetime

import data_utility as dut

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import networkx as nx


def unique_fname(name, ext, outdir="figures"):
    os.makedirs(outdir, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%dT%H%M%SZ")
    return os.path.join(outdir, f"{name}_{ts}.{ext}")


# fig is a plotly.graph_objects.Figure or plotly.express figure
def save_plotly(fig, name):
    html_path = unique_fname(name, "html")
    pio.write_html(fig, file=html_path, include_plotlyjs='cdn', auto_open=True)
    print("Saved:", html_path)

# # plt.figure() ... plotting code ...
def save_matplotlib(plt, name):
    png_path = unique_fname(name, "png")
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved:", png_path)


def show_simple_barplot(df, x_name, cut_off, name=False, save=True):
    print("simple_barplot")
    data = dut.cut_off(df, cut_off, x_name)
    plt.figure(figsize=(60,40)) #just to make sure
    sns.barplot(x=data.index, y=data.values, palette='viridis')
    plt.xlabel(x_name)
    plt.ylabel("count")
    plt.title(f"{x_name} with count > {cut_off}")
    plt.tight_layout()
    plt.show()
    if save:
        name = name if name else f"barplot_cutoff_{cut_off}"
        save_matplotlib(plt, name)


def show_better_barplot(df, x_name, cut_off, name=False, save=True):
    print("better_barplot")
    data = dut.cut_off(df, cut_off, x_name)
    #top_df = cut_off.reset_index()#.rename(columns={'index':'genes','genes':'count'})
    fig = px.bar(data, x=x_name, y='count', orientation='h',
                title=f'genes with count > {cut_off}', height=600)
    fig.update_layout(yaxis={'categoryorder':'total ascending'})
    #fig.show()
    if save:
        name = name if name else f"barplot_cutoff_{cut_off}"
        save_plotly(fig, name)


def plot_incidence_heatmap(df, gene_col, y_name, gene_cutoff=25, top_k_groups=50, name=False, save=True):
    print("heatmap")
    # select genes above cutoff and the top_k from another colum
    data = df[[gene_col, y_name]]
    genes = dut.cut_off(data, gene_cutoff, gene_col).index.tolist()
    topky = data[y_name].value_counts().index[:top_k_groups].tolist()

    data = data[data[gene_col].isin(genes) & data[y_name].isin(topky)]

    # # pivot to incidence (binary)
    incidence = pd.crosstab(data[y_name], data[gene_col]).clip(upper=1)
    # # optionally keep top groups by sum
    group_sums = incidence.sum(axis=0).sort_values(ascending=False)
    groups_keep = group_sums.index[:top_k_groups]
    incidence = incidence[groups_keep]

    # plot
    plt.figure(figsize=(60,40))
    sns.heatmap(incidence, cmap='YlOrBr', cbar=False)
    plt.xlabel('genes ('+gene_col+')')
    plt.ylabel('groups ('+y_name+')')
    plt.title(f'heatmap genes with count > {gene_cutoff} and the top {top_k_groups} most named groups)')
    plt.tight_layout()
    plt.show()

    if save:
        name = name if name else f"heatmap_cutoff_{gene_cutoff}_topk_{top_k_groups}"
        save_matplotlib(plt, name)


def plot_bipartite_network(df, gene_col, group_col, special_genes=[], gene_cutoff=0, max_genes=50, max_groups=50, name=False, save=True):
    print("network")

    df_sub = df[[gene_col, group_col]]
    genes_keep = dut.cut_off(df_sub, gene_cutoff, gene_col, True).index.tolist()[:max_genes]
    groups_keep = df_sub[group_col].value_counts().index[:max_groups].tolist()

    df_sub = df_sub[df_sub[gene_col].isin(genes_keep) & df_sub[group_col].isin(groups_keep)]

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
        if u[1] in special_genes:
            edge_color.append('pink')
            edge_width.append(0.7)
        else:
            edge_color.append('#888')
            edge_width.append(0.5)

    node_x, node_y, node_text, node_color, node_size = [], [], [], [], []
    for n,data in G.nodes(data=True):
        x,y = pos[n]
        node_x.append(x); node_y.append(y)
        node_text.append(data['label'])
        if data['bipartite']==0:
            if data['label'] in special_genes:
                node_color.append('pink'); node_size.append(15)
            else:
                node_color.append('blue'); node_size.append(5)
        else:
            node_color.append('red'); node_size.append(10)

    edge_trace = go.Scatter(x=edge_x, y=edge_y, mode='lines', line=dict(width=0.5, color='#888'), hoverinfo='none')
    node_trace = go.Scatter(x=node_x, y=node_y, mode='markers+text', text=node_text,
                            marker=dict(color=node_color, size=node_size), textposition='top center')

    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(showlegend=False, title='Gene-Group bipartite network')
    #fig.show()
    if save:
        name = name if name else f"network_bipartit_graph_{gene_cutoff}_{max_genes}_{max_groups}"
        save_plotly(fig, name)



def sankey_genes_groups(df, gen_column, terme_general, terme_specific, gene_cutoff=10, top_genes=100, top_general=50, top_specific=50, name=False, save=True):
    print("samley groups")

    df_sub = df[[gen_column, terme_general, terme_specific]]
    genes_keep = dut.cut_off(df_sub, gene_cutoff, gen_column, True).index.tolist()[:top_genes]
    df_sub = df_sub[df_sub[gen_column].isin(genes_keep)]

    general_top = df_sub[terme_general].value_counts().index[:top_general].tolist()
    specific_top = df_sub[terme_specific].value_counts().index[:top_specific].tolist()
    df_sub = df_sub[df_sub[terme_general].isin(general_top) & df_sub[terme_specific].isin(specific_top)]

    # build nodes
    nodes = genes_keep + general_top + specific_top
    node_idx = {n:i for i,n in enumerate(nodes)}

    # links specific -> genes
    df_s = df_sub.drop_duplicates([terme_specific,gen_column]).groupby([terme_specific,gen_column]).size().reset_index(name='count')

    # links genes -> general
    df_g = df_sub.drop_duplicates([gen_column,terme_general]).groupby([gen_column,terme_general]).size().reset_index(name='count')

    source, target, value = [], [], []
    for _,r in df_s.iterrows():
        source.append(node_idx[r[terme_specific]]); target.append(node_idx[r[gen_column]]); value.append(r['count'])

    for _,r in df_g.iterrows():
        source.append(node_idx[r[gen_column]]); target.append(node_idx[r[terme_general]]); value.append(r['count'])

    fig = go.Figure(go.Sankey(node=dict(label=nodes), link=dict(source=source, target=target, value=value)))
    fig.update_layout(title_text="Sankey: specific → genes → general", font_size=20)
    #fig.show()
    if save:
        name = name if name else f"sankey_{gene_cutoff}_{top_genes}_{top_general}_{top_specific}"
        save_plotly(fig, name)
# usage

def generate_all(df, amigodf, gen_cut_off=0, top_genes=200, top_groups_general=100, top_groups_specific=100):
    rename_dict = {"gene_symbol": "gene", "bioentity_label": "gene", 
                   "group_term" : "term_general", 
                   "term" : "term_specific", "disease_name": "term_specific"}
    gen_names = "gene"
    term_general = "term_general"
    term_specific = "term_specific"

    #to make sure we don't accidentally break the table we copy here
    sub_df = df.rename(mapper=rename_dict, axis=1) #implicit copy

    # try:
    #     show_better_barplot(sub_df, gen_names, gen_cut_off)
    # except Exception as e:
    #     print(str(e))

    # try:
    #     plot_incidence_heatmap(sub_df, gen_names, term_general, gen_cut_off, top_groups_general)
    # except Exception as e:
    #     print(str(e))

    # try:
    #     plot_incidence_heatmap(sub_df, gen_names, term_specific, gen_cut_off, top_groups_specific)
    # except Exception as e:
    #     print(str(e))

    try:
        plot_bipartite_network(sub_df, gen_names, term_general, 
                               gene_cutoff=5, max_genes=top_genes, max_groups=20)
    except Exception as e:
        print(str(e))

    try:
        plot_bipartite_network(sub_df, gen_names, term_specific, gen_cut_off, top_genes, top_groups_specific)
    except Exception as e:
        print(str(e))

    try:
        sankey_genes_groups(sub_df, gen_names, term_general, term_specific, gen_cut_off, top_genes, top_groups_general, top_groups_specific)
    except Exception as e:
        print(str(e))

#show_better_barplot(threshold_genes)
#plot_incidence_heatmap(df_general, group_col='terme_specific', gene_cutoff=20, top_k_groups=40)
#plot_incidence_heatmap(df_general, group_col='terme_general', gene_cutoff=20, top_k_groups=40)
#plot_bipartite_network(df_general, group_col='terme_specific', gene_cutoff=20)
#plot_bipartite_network(df_general, group_col='terme_general', gene_cutoff=20)


#sankey_genes_groups(df_general, gene_cutoff=20)


##not used
# from sklearn.metrics import jaccard_score
# from sklearn.cluster import AgglomerativeClustering
# import scipy.cluster.hierarchy as sch
# def gene_overlap_clustering(df, group_col='terme_specific', gene_cutoff=20, top_n_genes=60):
#     print("cluster")
#     counts = df[gen_column].value_counts()
#     genes_keep = counts[counts >= gene_cutoff].index[:top_n_genes].tolist()
#     incidence = pd.crosstab(df[df[gen_column].isin(genes_keep)][gen_column], df[group_col]).clip(upper=1)
#     # compute Jaccard distance matrix
#     M = 1 - (incidence.dot(incidence.T) / (incidence.sum(axis=1).values[:,None] + incidence.sum(axis=1).values[None,:] - incidence.dot(incidence.T))).fillna(0)
#     # hierarchical clustering dendrogram
#     plt.figure(figsize=(10,6))
#     sch.dendrogram(sch.linkage(M, method='average'), labels=incidence.index.tolist(), orientation='right')
#     plt.title('Hierarchical clustering of genes by shared groups (Jaccard)')
#     plt.tight_layout()
#     plt.show()

# usage

# from upsetplot import from_indicators, plot
# def gpt_makes_stuff(df_general, group_col):
# # choose a few top groups
#     print("gpt")
#     top_groups = df_general[group_col].value_counts().index[:10].tolist()
#     inc = pd.crosstab(df_general['genes'], df_general[group_col]).loc[:, top_groups].clip(upper=1)
#     up = from_indicators(inc.columns.tolist(), inc.reset_index(drop=True).values)
#     plot(up)
#     plt.title('UpSet plot: intersections among top groups')
#     plt.show()


#gene_overlap_clustering(df_general, group_col='terme_specific', gene_cutoff=20, top_n_genes=50)
#gene_overlap_clustering(df_general, group_col='terme_general', gene_cutoff=20, top_n_genes=50)

# gpt_makes_stuff(df_general, 'terme_specific')
# gpt_makes_stuff(df_general, 'terme_general')