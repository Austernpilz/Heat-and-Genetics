import os
import pandas as pd

#utiliti_functions
def get_col_as_unique_and_count(df, name):
    return df[name].value_counts()

def get_col_as_list(df, name):
    return df[name].tolist()

def get_col_unique(df, name):
    return df[name].unique().tolist()


def include_col(df, col, include):
    return df[df[col].isin(include)]


def exclude_col(df, col, exclude):
    return df[~df[col].isin(exclude)]


def cut_off(df, cutoff, col, ascend=False):
    new_df = get_col_as_unique_and_count(df, col)
    return new_df[new_df>=cutoff].sort_values(ascending=ascend)


def top_cut(df, cutoff, col):
    top_x = df[col].value_counts().index.tolist()[:cutoff]
    return df[df[col].isin(top_x)]


def new_table(df, col_list, rename_dict):
    filtered_list = [col for col in col_list if col in df.columns]
    return df[filtered_list].rename(mapper=rename_dict, axis=1, errors='ignore') #copy's by default

def make_new_table(df_list, cols_list, rename_dict, axis=0):
    # axis=0 means I assume you condense your dict, 
    # and want to merge different columns into one
    # so it will concat on index (under the other data)

    # this is not efficient for larger datasets but quick to implement
    return_list = []
    for df in df_list:
        return_list.append(new_table(df, cols_list, rename_dict))
    return pd.concat(return_list, axis=axis, copy=False) #it always uses new_tabel, so copy=False is safe here

def apply_include_exclude_txt(path, df, colum):
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

    df = exclude_col(df, colum, in_out_txt["exclude"]).copy()
    df["group_term"] = df[colum].map(lambda x : in_out_txt["group"].get(x, "NO_GROUP_TERM"))
    print(df[df["group_term"] == "NO_GROUP_TERM"][colum].unique().tolist())

    # .copy() to make sure, we don't corrupt the original data, 
    # when performing later transformations
    df_reduced = exclude_col(df, colum, in_out_txt["plusplus"]).copy()

    return df_reduced, df


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

        gene_name = entry.name
        gene_folder = os.path.join(input_path, entry.name)
        for file in os.scandir(gene_folder):
            if not file.is_file():
                continue
            if not file.name.startswith("variants_") and not file.name.endswith(".tsv"):
                continue

            try:
                cleaned_variant_table = os.path.join(gene_folder, file.name)
                df = pd.read_csv(cleaned_variant_table, sep='\t')
                df["gene"] = gene_name
                chromosomes += df["chrom"].unique().tolist()
                returndf.append(df)
            except Exception as e:
                print(str(e))

            break

    return pd.concat(returndf, ignore_index=True), list(set(chromosomes))


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
    #AF_afr_exome = df["exome.af_afr"].tolist()
    # AF_nfe_genome = df["genome.af_afr"].tolist()
    AF_nfe_joint = df["joint.af_nfe"].tolist()
    #AF_nfe_exome = df["exome.af_nfe"].tolist()
    #AF_nfe_genome = df["genome.af_nfe"].tolist()

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


def get_vcf(path_input, path_output):
    df, chromosomes = build_table(path_input)
    sub_df = df.drop_duplicates(ignore_index= True)
    head = build_header(chromosomes)
    with open (path_output, "w+") as f:
        f.write(head)
        for line in get_line(sub_df):
            f.write("\n")
            f.write(line)


# ID = gnomAD ID or rsID
# chrom, pos, ref, alt = take from db
#qual, filter = .
##INFO=<ID=NS,Number=1,Type=Float,Description="Number of Samples With Data">
# Info = AF_afr={af_afr};AF_nfe={af_nfe},

# CADD Score in bins 5er Schritte bis 20
# CADD Score ab 10,15,20
# Gene mit höchsten Score, 

# Sankey Gene mit CADD Score Bin, Anzahl an Varianten

# Die anderen Spalten (QUAL, FILTER, INFO, FORMAT) können alle einen "." (Punkt) haben und dann passt das.

