import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.io as pio
from datetime import datetime
import threading
#from src.helpers.figures import figures as fig
#utiliti_functions

from src.helpers.folder_magic import search_for_files
from src.HGNC.search_and_fetch import get_tables
#from src.disgnet.get_tables import build_tables


def unique_fname(name, ext, outdir="results/figures"):
    outdir= os.path.join(outdir, "data")
    os.makedirs(outdir, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%dT%H%M%SZ")
    return os.path.join(outdir, f"{name}_{ts}.{ext}")

def save_figure(fig, name, outdir="results/figures"):
    """
    Universal save function. Detects if the figure is Matplotlib or Plotly 
    and routes it to the correct saving mechanism.
    """
    # Get the base module name of the figure object (e.g., 'matplotlib' or 'plotly')
    fig_module = type(fig).__module__.split('.')[0]
    
    if fig_module == 'matplotlib':
        # --- Handle Matplotlib (saves as PNG) ---
        import matplotlib.pyplot as plt
        png_path = unique_fname(name, "png")
        fig.savefig(png_path, dpi=300, bbox_inches="tight")
        plt.close(fig)  # Free memory (critical for multithreading)
        print(f"Saved Matplotlib: {png_path}")
        
    elif fig_module == 'plotly':
        # --- Handle Plotly (saves as interactive HTML) ---
        import plotly.io as pio
        html_path = unique_fname(name, "html")
        
        # NOTE: auto_open=False is crucial here. If it is True in a threaded 
        # environment, your script will aggressively open 20+ browser tabs at once.
        pio.write_html(fig, file=html_path, include_plotlyjs='cdn', auto_open=False)
        print(f"Saved Plotly: {html_path}")

def get_config(config):
    abs = config.get("absolute_file_paths")
    rel = config.get("relative_file_paths")

    data = abs.get("data")
    data_amigo = os.path.join(data, "AmiGo2", "genes")
    data_disgnet = os.path.join(data, "disgnet", "genes")
    data_hgnc = os.path.join(data, "HGNC", "genes")

    config_path = abs.get("config")
    in_ex_group_amigo = os.path.join(config_path, rel.get("Amigo2_inclue_exclude"))
    in_ex_group_disgnet = os.path.join(config_path, rel.get("disgnet_include_exclude"))

    results = abs.get("results")

    return data_amigo, data_disgnet, data_hgnc, in_ex_group_amigo, in_ex_group_disgnet, results

def get_disgnet_df(data_disgnet, in_ex_group_disgnet):
    ieg_disgnet = load_include_exclude_txt(in_ex_group_disgnet)
    resulttsv = search_for_files(data_disgnet, "search_result", "tsv")
    df_list = []
    for tsv in resulttsv:
        try:
            df = pd.read_csv(tsv, sep="\t", dtype=str)
            if df is None:
                continue
            if df.empty:
                continue
            df["hgnc_id"] = ""
            df["name"] = ""
            df["symbol"] = df["gene_symbol"].astype(str).str.strip().str.upper()
            df["ensembl_gene_id"] = df["geneEnsemblIDs"].astype(str, errors='ignore')
            df.rename(columns={"disease_name": "term"}, inplace=True)
            df["group"] = df["term"].map(lambda x : ieg_disgnet["group"].get(x, "NO GROUP"))
            sub_df = df[["symbol", "name", "term", "group", "hgnc_id", "ensembl_gene_id"]].copy()
            df_list.append(sub_df)
        except Exception as e:
            print(str(e))
    df = pd.concat(df_list)
    return df

def get_amigo_df(data_amigo, in_ex_group_amigo):
    ieg_amigo = load_include_exclude_txt(in_ex_group_amigo)
    amigo_tsv = search_for_files(data_amigo, "data", "tsv")
    amigo_list = []
    for table in amigo_tsv:
        term_name = os.path.basename(os.path.basename(table))
        if term_name in ieg_amigo["exclude"] or term_name in ieg_amigo["plusplus"]:
            continue
        try:
            df = pd.read_csv(table, sep="\t", dtype=str)
            if df is None:
                continue
            if df.empty:
                continue
            df["hgnc_id"] = ""
            df["ensembl_gene_id"] = ""
            df["symbol"] = df["bioentity_label"].astype(str).str.strip().str.upper()
            df["name"] = df["bioentity_name"].astype(str).str.strip().str.lower()
            try:
                sub_df = df[["symbol", "name", "term", "group", "hgnc_id", "ensembl_gene_id"]].copy()
                amigo_list.append(sub_df)
            except Exception as _:
                df["term"] = term_name
                df["group"] = ieg_amigo["group"].get(term_name, "NO GROUP")
                sub_df = df[["symbol", "name", "term", "group", "hgnc_id", "ensembl_gene_id"]].copy()
                amigo_list.append(sub_df)
        except Exception as _:
            continue
    return pd.concat(amigo_list).drop_duplicates()

def get_hgnc_look_up(data_hgnc):
    df = get_tables(data_hgnc)

    # Standardize formats
    df["symbol_norm"] = df["symbol"].astype(str).str.strip().str.upper()
    df["name_norm"] = df["name"].astype(str).str.strip().str.lower()
    df["ensembl_gene_id_norm"] = df["ensembl_gene_id"].astype(str, errors='ignore')

    # Drop duplicates before setting the index to ensure uniqueness
    hgnc_symbol_look_up = (
        df.drop_duplicates(subset=["symbol_norm"])
        .set_index("symbol_norm")[["name", "hgnc_id", "ensembl_gene_id"]]
        .to_dict(orient="index")
    )

    hgnc_name_look_up = (
        df.drop_duplicates(subset=["name_norm"])
        .set_index("name_norm")[["name", "hgnc_id", "ensembl_gene_id"]]
        .to_dict(orient="index")
    )

    hgnc_ensembl_look_up = (
        df.drop_duplicates(subset=["ensembl_gene_id_norm"])
        .set_index("ensembl_gene_id_norm")[["name", "hgnc_id", "ensembl_gene_id"]]
        .to_dict(orient="index")
    )

    return hgnc_symbol_look_up, hgnc_name_look_up, hgnc_ensembl_look_up

def fill_up(df, hgnc):
    hgnc_symbol_look_up, hgnc_name_look_up, hgnc_ensembl_look_up = hgnc
    df = df.copy()
    def pick(existing, new):
        if pd.isna(existing):
            return new
        if str(existing).strip() == "":
            return new
        return existing

    def fill_row(row):
        key0 = row.get("symbol", "")
        info0 = hgnc_symbol_look_up.get(key0, {})
        if key0 and info0:
            row["name"] = pick(row.get("name", ""), info0.get("name", ""))
            row["hgnc_id"] = pick(row.get("hgnc_id", ""), info0.get("hgnc_id", ""))
            row["ensembl_gene_id"] = pick(row.get("ensembl_gene_id", ""), info0.get("ensembl_gene_id", ""))

        key1 = row.get("name", "")
        info1 = hgnc_name_look_up.get(key1, {})
        if key1 and info1:
            row["symbol"] = pick(row.get("symbol", ""), info1.get("symbol", ""))
            row["hgnc_id"] = pick(row.get("hgnc_id", ""), info1.get("hgnc_id", ""))
            row["ensembl_gene_id"] = pick(row.get("ensembl_gene_id", ""), info1.get("ensembl_gene_id", ""))

        key2 = row.get("ensembl_gene_id", "")
        info2 = hgnc_ensembl_look_up.get(key2, {})
        if key2 and info2:
            row["name"] = pick(info2.get("name", ""), row.get("name", ""))
            row["hgnc_id"] = pick(info2.get("hgnc_id", ""), row.get("hgnc_id", ""))
            row["symbol"] = pick(info2.get("symbol", ""), row.get("symbol", ""))
        return row
    #df[["symbol", "name", "term", "group", "hgnc_id", "ensembl_gene_id"]]
    df = df.apply(fill_row, axis=1)
    return df

def build_look_up(config_file):
    data_amigo, data_disgnet, data_hgnc, in_ex_group_amigo, in_ex_group_disgnet, results = get_config(config_file)
    hgnc_df = get_hgnc_look_up(data_hgnc)
    amigo_df = fill_up(get_amigo_df(data_amigo, in_ex_group_amigo), hgnc_df)
    disgnet_df = fill_up(get_disgnet_df(data_disgnet, in_ex_group_disgnet), hgnc_df)
    print(amigo_df.head(2))
    print(disgnet_df.head(2))
    try:
        df = pd.concat([amigo_df, disgnet_df])
        file_name = os.path.join(results, "genes_terms_groups.tsv")
        df.to_csv(file_name, sep='\t', index=False)

    except Exception as e:
        print(str(e))
        file_name_amigo = os.path.join(results, "amigo_genes_terms_groups.tsv")
        file_name_disgnet = os.path.join(results, "disgnet_genes_terms_groups.tsv")
        amigo_df.to_csv(file_name_amigo, sep='\t', index=False)
        disgnet_df.to_csv(file_name_disgnet, sep='\t', index=False)

def focus_on_variants(config):
    abs = config.get("absolute_file_paths")
    data_variants = os.path.join(abs.get("data"), "gnomAD")
    results = abs.get("results")
    list_list = []
    for look_up in search_for_files(results, "genes_terms_groups", "tsv"):
        df = pd.read_csv(look_up, sep='\t')
        for ensembl_id in df["ensembl_gene_id"].unique():
            if isinstance(ensembl_id, float):
                continue
            genes = search_for_files(os.path.join(data_variants, ensembl_id), "simple_variants", "tsv")
            if len(genes) >= 1:
                df = df[df["ensembl_gene_id"] != ensembl_id]
        df.to_csv(look_up, sep='\t', index=False)
        list_list.append(df)
    df_group = list_list[0] if len(list_list) == 1 else pd.concat(list_list, ignore_index=True)
    save_var = os.path.join(results, "all_variants.tsv")
    df_var = save_variants_to_csv(data_variants, save_var)

    #vcf = os.path.join(results, "all_variants.vcf")
    #build_vcf(data_variants, vcf)

    save_clinvar = os.path.join(results, "all_variants.vcf")
    df_clinvar = save_clinvar_to_csv(data_variants ,save_clinvar)

    run_analysis_pipeline(df_var, df_clinvar, df_group)
    return

def save_results(config_file):
    results = config_file.get("absolute_file_paths").get("results")
    build_look_up(config_file)
    df = focus_on_variants(config_file)



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# 1. DATA MERGING & CLEANING
# ==========================================

def build_merged_dataframe(df_variants, df_clinvar, df_genes):
    """Merges the three tables into a single dataframe based on IDs."""
    # Merge variants with clinvar on variant_id (Left merge preserves all variants)
    # Suffixes handle overlapping columns like hgvsc and hgvsp
    merged = pd.merge(df_variants, df_clinvar, on="variant_id", how="left", suffixes=("", "_clinvar"))

    # Merge with genes on symbol
    merged = pd.merge(merged, df_genes, left_on="gene_symbol", right_on="symbol", how="outer", suffixes=("", "_terms"))
    return merged

def get_any_deterious(df):
    """
    Returns a dataframe of variants that trigger AT LEAST ONE deleterious threshold.
    Adds a 'deleterious_hits' column showing how many thresholds were crossed.
    """
    # Create boolean masks for each dangerous condition
    c1 = df["revel"].fillna(0) >= 0.5
    c2 = df["cadd"].fillna(0) >= 20.0
    c3 = df["pangolin"].fillna(0) >= 0.2
    c4 = df["spliceai"].fillna(0) >= 0.2
    c5 = df["sift"].fillna(1.0) <= 0.05  # SIFT is bad when LOW
    c6 = df["polyphen"].fillna(0) >= 0.85
    c7 = df["gold_stars"].fillna(0) >= 2
    
    sig_mapping = lambda x: 0 if "benign" in str(x).lower() else 1 #{"Benign": 0, "Uncertain significance": 1, "Likely pathogenic": 2, "Pathogenic": 3}
    c8 = df["clinical_sig"].map(sig_mapping).fillna(0) == 1
    
    high_impact = ["stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant"]
    c9 = df["major_consequence"].isin(high_impact)
    
    # Sum the booleans to get a total burden score (1 to 9)
    df_scored = df.copy()
    df_scored["deleterious_hits"] = (
        c1.astype(int) + c2.astype(int) + c3.astype(int) + 
        c4.astype(int) + c5.astype(int) + c6.astype(int) + 
        c7.astype(int) + c8.astype(int) + c9.astype(int)
    )
    
    # Filter for variants that hit at least one threshold
    # Drop rows that don't have a variant_id (e.g., orphan gene rows from the outer merge)
    any_deleterious = df_scored[(df_scored["deleterious_hits"] >= 1) & (df_scored["variant_id"].notna())].copy()
    
    return any_deleterious

def extract_and_rename_subtable(df):
    """Selects the requested columns and renames them for easier access."""
    column_mapping = {
        "symbol": "gene_symbol",
        "term": "term",
        "group": "group",
        "hgvsc": "hgvsc",
        "hgvsp": "hgvsp",
        "variant_id": "variant_id",
        
        # Population Allele Frequencies
        "joint.af_nfe": "joint_nfe",
        "joint.af_afr": "joint_afr",
        "exome.af_nfe": "exome_nfe",
        "exome.af_afr": "exome_afr",
        "genome.af_nfe": "genome_nfe",
        "genome.af_afr": "genome_afr",
        
        # 6 in-silico predictors
        "in_silico_predictors.revel_max": "revel",
        "in_silico_predictors.cadd": "cadd",
        "in_silico_predictors.pangolin_largest_ds": "pangolin",
        "in_silico_predictors.spliceai_ds_max": "spliceai",
        "in_silico_predictors.sift_max": "sift",
        "in_silico_predictors.polyphen_max": "polyphen",
        
        # 3 Clinvar columns
        "major_consequence": "major_consequence",
        "clinical_significance": "clinical_sig",
        "gold_stars": "gold_stars"
    }
    
    # Keep only columns that exist in mapping to avoid KeyError
    existing_cols = [col for col in column_mapping.keys() if col in df.columns]
    sub_df = df[existing_cols].copy()
    
    # Rename
    sub_df.rename(columns=column_mapping, inplace=True)
    return sub_df

# ==========================================
# 2. FILTERING DETERIOUS SUBSETS
# ==========================================

def get_deterious_subsets(df):
    """
    Splits the main df into 9 subsets where the specific score is in the 'deterious' range.
    Returns a dictionary of dataframes.
    """
    subsets = {}

    # 1. REVEL (High is bad, > 0.5 is generally pathogenic)
    subsets["revel"] = df[df["revel"] >= 0.5].copy()

    # 2. CADD (High is bad, > 20 is top 1% deleterious)
    subsets["cadd"] = df[df["cadd"] >= 20.0].copy()

    # 3. Pangolin (High is bad, > 0.2 indicates splice altering)
    subsets["pangolin"] = df[df["pangolin"] >= 0.2].copy()

    # 4. SpliceAI (High is bad, > 0.2 is potential splice altering)
    subsets["spliceai"] = df[df["spliceai"] >= 0.2].copy()

    # 5. SIFT (LOW is bad, < 0.05 is deleterious)
    # We invert SIFT so "higher is darker" works across all plots automatically
    df_sift = df[df["sift"] <= 0.05].copy()
    df_sift["sift_inverted"] = 1.0 - df_sift["sift"]
    subsets["sift"] = df_sift
    
    # 6. PolyPhen (High is bad, > 0.85 is probably damaging)
    subsets["polyphen"] = df[df["polyphen"] >= 0.85].copy()
    
    # 7. Gold Stars (Higher is better clinically, let's say >= 2 stars is reliable)
    subsets["gold_stars"] = df[df["gold_stars"] >= 2].copy()
    
    # 8. Clinical Significance (String: needs mapping for "darker colors")
    sig_mapping = lambda x: 0 if "benign" in str(x).lower() else 1 #{"Benign": 0, "Uncertain significance": 1, "Likely pathogenic": 2, "Pathogenic": 3}
    df_clin = df.dropna(subset=["clinical_sig"]).copy()
    df_clin["clinical_sig_num"] = df_clin["clinical_sig"].map(sig_mapping)
    subsets["clinical_sig"] = df_clin[df_clin["clinical_sig_num"] == 2].copy() # Likely Pathogenic or Pathogenic
    
    # 9. Major Consequence (String: mapping severity)
    # Assuming High impact variants like stop_gained, frameshift, splice site
    high_impact = ["stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant"]
    df_cons = df.dropna(subset=["major_consequence"]).copy()
    # Map them to 1 so they can be colored
    df_cons["major_consequence_num"] = df_cons["major_consequence"].apply(lambda x: 1 if x in high_impact else 0)
    subsets["major_consequence"] = df_cons[df_cons["major_consequence_num"] == 1].copy()
    
    return subsets

# ==========================================
# 3. PLOTTING FUNCTIONS
# ==========================================

def calculate_dynamic_figsize(num_rows):
    """Adjusts plot size based on the number of data points."""
    base_size = 6
    if num_rows < 50:
        return (base_size, base_size)
    elif num_rows < 500:
        return (base_size + 2, base_size + 2)
    else:
        return (base_size + 4, base_size + 4)

def plot_population_scatter(df, score_col, target_color_col, title):
    """
    Plots joint_nfe vs joint_afr. Darker colors = higher severity score.
    """
    if df.empty:
        print(f"Skipping Scatter Plot for {title}: Dataframe is empty.")
        return
    
    fig_size = calculate_dynamic_figsize(len(df))
    plt.figure(figsize=fig_size)
    
    # Set colormap: Blues, but reversed so higher number = darker color
    cmap = "Reds"
    
    scatter = plt.scatter(
        x=df["joint_nfe"], 
        y=df["joint_afr"], 
        c=df[target_color_col], 
        cmap=cmap, 
        alpha=0.7, 
        edgecolors='k',
        s=50 # Point size
    )
    
    plt.colorbar(scatter, label=f"Severity ({target_color_col})")
    
    # Diagonal reference line (equal frequency)
    max_val = max(df["joint_nfe"].max(), df["joint_afr"].max()) * 1.1
    if not pd.isna(max_val):
        plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label="Equal Frequency")

    plt.xlabel("Joint Allele Frequency (NFE)")
    plt.ylabel("Joint Allele Frequency (AFR)")
    plt.title(f"{title} (n={len(df)})")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()
    save_figure(plt, score_col)

def plot_group_majority(df, title):
    """
    Maps gene groups to the population where they are more frequent 
    and visualizes the counts.
    """
    if df.empty or "group" not in df.columns:
        print(f"Skipping Group Plot for {title}: Missing data.")
        return
        
    # Determine majority population per variant
    df_clean = df.dropna(subset=["joint_nfe", "joint_afr", "group"]).copy()
    if df_clean.empty:
        return
        
    df_clean["majority_pop"] = np.where(df_clean["joint_nfe"] > df_clean["joint_afr"], "NFE Majority", "AFR Majority")
    
    # Group and count
    group_counts = df_clean.groupby(["group", "majority_pop"]).size().unstack(fill_value=0)
    
    # Plot
    fig_size = calculate_dynamic_figsize(len(group_counts) * 10) # adjust based on number of groups
    group_counts.plot(kind="bar", stacked=True, figsize=(fig_size[0], 5), color=["#E69F00", "#56B4E9"])
    
    plt.title(f"Gene Group Majority Distribution - {title}")
    plt.xlabel("Gene Group")
    plt.ylabel("Number of Variants")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()
    save_figure(plt, title)


# ==========================================
# 4. MAIN ORCHESTRATOR
# ==========================================

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

def run_analysis_pipeline(df_variants, df_clinvar, df_genes):
    print("1. Merging datasets...")
    merged_df = build_merged_dataframe(df_variants, df_clinvar, df_genes)

    print("2. Extracting and renaming sub-table...")
    clean_df = extract_and_rename_subtable(merged_df)

    print("3. Generating deterious subsets...")
    subsets = get_deterious_subsets(clean_df)
    all_negativ = get_any_deterious(clean_df).copy()
    
    print("4. Plotting...")
    # Mapping the dataframe key to the column that contains the numeric value for coloring
    color_cols = {
        "revel": "revel",
        "cadd": "cadd",
        "pangolin": "pangolin",
        "spliceai": "spliceai",
        "sift": "sift_inverted",           # We use the inverted score so darker = worse
        "polyphen": "polyphen",
        "gold_stars": "gold_stars",
        "clinical_sig": "clinical_sig_num", # Custom numeric map created during filtering
        "major_consequence": "major_consequence_num" # Custom numeric map created during filtering
    }

    for key, sub_df in subsets.items():
        print(f"--- Processing {key.upper()} ---")
        plot_population_scatter(sub_df, key, color_cols[key], f"{key.upper()} deterious Variants")
        plot_group_majority(sub_df, f"{key.upper()}_groups")

    plot_group_majority(all_negativ, "difference_groups")
    plot_deleterious_by_gene_group(all_negativ)
    plot_any_deleterious_burden(all_negativ)
    get_cadd_sankey(all_negativ)
    get_biggest_outliers(all_negativ)

def load_include_exclude_txt(path):
    in_out_txt = {
        "group" : {},
        "exclude" : [],
        "plusplus" : []
    }
    with open(path, 'r') as f:
        group = ""
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            elif line.startswith('<'): 
                line = line[1:-1].strip()
                group = line
                continue
            elif line.startswith('--'): 
                line = line[2:].strip()
                in_out_txt["exclude"].append(line)
            elif line.startswith('++'): 
                line = line[2:].strip()
                in_out_txt["plusplus"].append(line)
                in_out_txt["group"][line] = group
            else:
                in_out_txt["group"][line] = group

    # df = exclude_col(df, colum, in_out_txt["exclude"]).copy()
    # df["group_term"] = df[colum].map(lambda x : in_out_txt["group"].get(x, "NO_GROUP_TERM"))
    # print(df[df["group_term"] == "NO_GROUP_TERM"][colum].unique().tolist())

    # # .copy() to make sure, we don't corrupt the original data, 
    # # when performing later transformations
    # df_reduced = exclude_col(df, colum, in_out_txt["plusplus"]).copy()

    return in_out_txt

def plot_any_deleterious_burden(df, title="Any Deleterious: Burden Score"):
    """
    Scatter plot of population frequencies, colored by how many deleterious 
    thresholds (1-9) the variant crossed.
    """
    if df.empty:
        print("Skipping Burden Plot: Dataframe is empty.")
        return
        
    fig_size = calculate_dynamic_figsize(len(df))
    fig, ax = plt.subplots(figsize=fig_size)
    
    # Sort so the highest burden scores are plotted last (on top)
    df_sorted = df.sort_values("deleterious_hits")
    
    scatter = ax.scatter(
        x=df_sorted["joint_nfe"], 
        y=df_sorted["joint_afr"], 
        c=df_sorted["deleterious_hits"], 
        cmap="inferno_r", # Dark/intense colors for higher hits
        alpha=0.8, 
        edgecolors='k',
        s=df_sorted["deleterious_hits"] * 20 + 20 # Size scales with hits
    )
    
    cbar = fig.colorbar(scatter, ax=ax, ticks=range(1, 10))
    cbar.set_label("Number of Deleterious Thresholds Crossed")
    
    max_val = max(df["joint_nfe"].max(), df["joint_afr"].max()) * 1.1
    if not pd.isna(max_val):
        ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label="Equal Frequency")
        
    ax.set_xlabel("Joint Allele Frequency (NFE)")
    ax.set_ylabel("Joint Allele Frequency (AFR)")
    ax.set_title(f"{title} (n={len(df['variant_id'].unique())} unique variants)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    
    save_figure(fig, "any_deleterious_burden_scatter")


def plot_deleterious_by_gene_group(df, title="Any Deleterious: Gene Group Distribution"):
    """
    Horizontal bar chart showing which gene groups contain the most variants
    that crossed at least one deleterious threshold.
    """
    if df.empty or "group" not in df.columns:
        print("Skipping Gene Group Plot: Missing data.")
        return
    
    # Drop rows with no gene group mapped
    df_clean = df.dropna(subset=["group"]).copy()
    if df_clean.empty:
        return
    
    # Count occurrences of each group
    # Because of the outer merge, a variant is counted in every group it belongs to
    group_counts = df_clean["group"].value_counts().head(20) # Top 20 for readability
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    group_counts.sort_values().plot(kind="barh", ax=ax, color="darkred", edgecolor="black", alpha=0.8)
    
    ax.set_title(f"Top 20 Gene Groups with Deleterious Variants")
    ax.set_xlabel("Total Deleterious Variants in Group")
    ax.set_ylabel("Gene Group")
    ax.grid(axis='x', linestyle='--', alpha=0.7)
    fig.tight_layout()
    
    save_figure(fig, "any_deleterious_gene_groups")

def build_table(input_path):
    returndf = []
    chromosomes = []
    if not os.path.isdir(input_path):
        return pd.DataFrame()

    for entry in os.scandir(input_path):
        if not entry.is_dir():
            continue

        if not entry.name.startswith("ENS"):
            continue
        gene = entry.name
        gene_folder = os.path.join(input_path, gene)
        for file in os.scandir(gene_folder):
            if not file.is_file():
                continue
            if not file.name.startswith("simple_variants") and not file.name.endswith(".tsv"):
                continue

            try:
                cleaned_variant_table = os.path.join(gene_folder, file.name)
                df = pd.read_csv(cleaned_variant_table, sep='\t')
                df["ensemble_gene_id"] = gene
                chromosomes += df["chrom"].astype(str).unique().tolist()
                returndf.append(df)
            except Exception as e:
                print(str(e))

            break

    chr = {str(x): 0 for x in chromosomes}
    return pd.concat(returndf, ignore_index=True), list(chr.keys())


def build_header(contig_chr):
    base = """\
##fileformat=VCFv4.1
##INFO=<ID=gnomadID,Number=1,Type=String,Description="gnomAD_ID to find the variant in"
##INFO=<ID=rsID,Number=1,Type=String,Description="rsID to identify the variant"
##INFO=<ID=hgvsc,Number=1,Type=String,Description="">
##INFO=<ID=AF_afr,Number=1,Type=String,Description="Allel Frequence of African Ancestry">
##INFO=<ID=AF_nfe,Number=1,Type=String,Description="Allel Frequence of European(non-finish) Ancestry">\
"""
    help = {"X": 24, "x": 24, "Y": 25, "y":25}
    contig_chr.sort(key=lambda x: int(x) if str(x).isdigit() else help[x])
    for chr in contig_chr:
        base += "\n##contig=<ID=%s>" % (chr,)
    #"CHROM", 
    columns = ["POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    base += "\n#CHROM"
    for col in columns:
        base += f"\t{col}"

    return base


def float_with_x_precision(number, precision=4):
    number = float(number)
    return f"{number:.{precision}f}"


def build_info(df):
    info_list = []

    gnomadID = df["variant_id"].tolist()
    rsID = df["rsids"].tolist()
    hgvcs = df["hgvsc"].tolist()
    AF_afr_joint = df["joint.af_afr"].tolist()
    AF_afr_exome = df["exome.af_afr"].tolist()
    AF_nfe_genome = df["genome.af_afr"].tolist()
    AF_nfe_joint = df["joint.af_nfe"].tolist()
    AF_nfe_exome = df["exome.af_nfe"].tolist()
    AF_nfe_genome = df["genome.af_nfe"].tolist()

    for i in range(df.shape[0]):
        af_afr = float_with_x_precision(AF_afr_joint[i],4)
        af_nfe = float_with_x_precision(AF_nfe_joint[i],4)
        line = f"\
gnomadID={str(gnomadID[i])};\
rsID={str(rsID[i])};\
hgvcs={str(hgvcs[i])};\
AF_afr={af_afr};\
AF_nfe={af_nfe}\
"
        info_list.append(line)

    return info_list


def get_line(df):
    chromosome = df["chrom"].tolist()
    pos = df["pos"].tolist()
    ref = df["ref"].tolist()
    alt = df["alt"].tolist()
    info = build_info(df)
    for i in range(df.shape[0]):
        line = f"{chromosome[i]}\t\
{pos[i]}\t\
.\t\
{ref[i]}\t\
{alt[i]}\t\
.\t\
.\t\
{info[i]}\t\
."
        yield line


def build_vcf(path_input, path_output):
    df, chromosomes = build_table(path_input)
    help = {"X": 24, "x": 24, "Y": 25, "y":25}
    sub_df = df.drop_duplicates(ignore_index= True)#.sort_values(
        #by=["chrom", "pos"],
        #key=(lambda x: int(x) if str(x).isdigit() else help.get(str(x), 99))
        #)
    head = build_header(chromosomes)
    with open (path_output, "w+") as f:
        f.write(head)
        for line in get_line(sub_df):
            f.write("\n")
            f.write(line)

def save_clinvar_to_csv(data_variants ,save_clinvar):
    list_var = []
    for clinvar_tsv in search_for_files(data_variants, "simple_clinvar", "tsv"):
        list_var.append(pd.read_csv(clinvar_tsv, sep='\t'))
    df = pd.concat(list_var)
    df.to_csv(save_clinvar, index=False)
    return df

def save_variants_to_csv(path_input, path_output):
    df, _ = build_table(path_input)
    df.drop_duplicates(ignore_index= True, inplace=True)
    #print(df["lof"].unique())
    #print(df["lof_curation"].unique())
    df.to_csv(path_output, sep='\t', index=False)
    return df



# def get_biggest_outliers(cute_df):


# def get_cadd_sankey(cute_df):
#     cute_df["cadd"].tolist()
#     cute_df["gene_symbol"].tolist()
#     cute_df["variant_id"].tolist()
#     return

"""
potential gnomad scores
'cadd_phred': {'Description': "Cadd Phred-like scores ('scaled C-scores') ranging from 1 to 99, based on the rank of each variant relative to all possible 8.6 billion substitutions in the human reference genome. Larger values are more deleterious.", 'Number': '1'}, 
'cadd_raw_score': {'Description': "Raw CADD scores are interpretable as the extent to which the annotation profile for a given variant suggests that the variant is likely to be 'observed' (negative values) vs 'simulated' (positive values). Larger values are more deleterious.", 'Number': '1'}, 
'pangolin_largest_ds': {'Description': "Pangolin's largest delta score across 2 splicing consequences, which reflects the probability of the variant being splice-altering", 'Number': '1'}, 
'phylop': {'Description': 'Base-wise conservation score across the 241 placental mammals in the Zoonomia project. Score ranges from -20 to 9.28, and reflects acceleration (faster evolution than expected under neutral drift, assigned negative scores) as well as conservation (slower than expected evolution, assigned positive scores).', 'Number': '1'}, 
'polyphen_max': {'Description': 'Score that predicts the possible impact of an amino acid substitution on the structure and function of a human protein, ranging from 0.0 (tolerated) to 1.0 (deleterious).  We prioritize max scores for MANE Select transcripts where possible and otherwise report a score for the canonical transcript.', 'Number': '1'}, 
'revel_max': {'Description': "The maximum REVEL score at a site's MANE Select or canonical transcript. It's an ensemble score for predicting the pathogenicity of missense variants (based on 13 other variant predictors). Scores ranges from 0 to 1. Variants with higher scores are predicted to be more likely to be deleterious.", 'Number': '1'}, 
'sift_max': {'Description': 'Score reflecting the scaled probability of the amino acid substitution being tolerated, ranging from 0 to 1. Scores below 0.05 are predicted to impact protein function. We prioritize max scores for MANE Select transcripts where possible and otherwise report a score for the canonical transcript.', 'Number': '1'}, 
'spliceai_ds_max': {'Description': "Illumina's SpliceAI max delta score; interpreted as the probability of the variant being splice-altering.", 'Number': '1'}
"""

# ID = gnomAD ID or rsID
# chrom, pos, ref, alt = take from db
#qual, filter = .
##INFO=<ID=NS,Number=1,Type=Float,Description="Number of Samples With Data">
# Info = AF_afr={af_afr};AF_nfe={af_nfe},

# CADD Score in bins 5er Schritte bis 20
# CADD Score ab 10,15,20
# Gene mit höchsten Score,
# major consequence variant data, clinical_significance


# Sankey Gene mit lof in silico predictors CADD Score Bin, Anzahl an Varianten

# Die anderen Spalten (QUAL, FILTER, INFO, FORMAT) können alle einen "." (Punkt) haben und dann passt das.

import plotly.graph_objects as go

def get_biggest_outliers(cute_df, top_n=20):
    """
    Extracts the variants with the highest CADD scores, including
    major consequence and clinical significance data.
    """
    cols_to_check = ["cadd", "variant_id", "gene_symbol", "major_consequence", "clinical_sig"]
    existing_cols = [c for c in cols_to_check if c in cute_df.columns]
    
    # Clean and convert CADD to float for sorting
    df_clean = cute_df.dropna(subset=['cadd']).copy()
    df_clean['cadd'] = pd.to_numeric(df_clean['cadd'], errors='coerce')
    
    # Sort descending
    outliers = df_clean.sort_values(by='cadd', ascending=False).head(top_n)
    
    print(f"\n--- Top {top_n} Outliers (Highest CADD Scores) ---")
    print(outliers[existing_cols].to_string(index=False))
    
    return outliers[existing_cols]

def get_cadd_sankey(cute_df):
    """
    Builds a Plotly Sankey diagram flowing from:
    Gene -> Major Consequence -> CADD Score Bin
    """
    # 1. Clean missing data for the flow path
    df = cute_df.dropna(subset=['gene_symbol', 'major_consequence', 'cadd']).copy()
    df['cadd'] = pd.to_numeric(df['cadd'], errors='coerce')

    # 2. Create CADD Bins (5er Schritte bis 20+)
    bins = [-1, 5, 10, 15, 20, float('inf')]
    labels = ["0-5", "5-10", "10-15", "15-20", ">20"]
    df['cadd_bin'] = pd.cut(df['cadd'], bins=bins, labels=labels)

    # 3. Create unique node names (add prefixes to prevent loops/merges if names overlap)
    df['src_node'] = "Gene: " + df['gene_symbol'].astype(str)
    df['mid_node'] = "Effect: " + df['major_consequence'].astype(str)
    df['tgt_node'] = "CADD Bin: " + df['cadd_bin'].astype(str)

    # 4. Aggregate the two flow steps
    flow1 = df.groupby(['src_node', 'mid_node']).size().reset_index(name='count')
    flow2 = df.groupby(['mid_node', 'tgt_node']).size().reset_index(name='count')

    # 5. Build node dictionary for Plotly
    all_nodes = list(pd.concat([flow1['src_node'], flow1['mid_node'], flow2['tgt_node']]).unique())
    node_indices = {name: i for i, name in enumerate(all_nodes)}

    source, target, value = [], [], []

    # Map Flow 1 (Gene -> Consequence)
    for _, row in flow1.iterrows():
        source.append(node_indices[row['src_node']])
        target.append(node_indices[row['mid_node']])
        value.append(row['count'])

    # Map Flow 2 (Consequence -> CADD Bin)
    for _, row in flow2.iterrows():
        source.append(node_indices[row['mid_node']])
        target.append(node_indices[row['tgt_node']])
        value.append(row['count'])

    # 6. Generate Figure
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=all_nodes,
            color="lightblue"
        ),
        link=dict(source=source, target=target, value=value)
    )])

    fig.update_layout(title_text="Variant Flow: Gene -> Major Consequence -> CADD Score", font_size=12)
    
    # Hooking into your existing Plotly save function
    save_figure(fig, "cadd_sankey_diagram")
    return fig