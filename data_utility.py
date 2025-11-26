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



def build_header(contig_chr):
    base = """
    ##fileformat=VCFv4.1
    ##INFO=<ID=gnomadID,Number=1,Type=String,Description="gnomAD_ID to find the variant in"
    ##INFO=<ID=AF_afr,Number=1,Type=Float,Description="Allel Frequence of African Ancestry">
    ##INFO=<ID=AF_nfe,Number=1,Type=Float,Description="Allel Frequence of European(non-finish) Ancestry">
    """
    contig_chr.sort(key=lambda x: x.isdigit())
    for chr in contig_chr:
        base += "\n##contig=<ID=%s>" % (chr,)
    ##contig=<ID=1>
    ##contig=<ID=2>
    ##contig=<ID=3>
    ##contig=<ID=4>
    ##contig=<ID=5>
    ##contig=<ID=6>
    ##contig=<ID=7>
    ##contig=<ID=8>
    ##contig=<ID=9>
    ##contig=<ID=10>
    ##contig=<ID=11>
    ##contig=<ID=12>
    ##contig=<ID=13>
    ##contig=<ID=14>
    ##contig=<ID=15>
    ##contig=<ID=16>
    ##contig=<ID=17>
    ##contig=<ID=18>
    ##contig=<ID=19>
    ##contig=<ID=20>
    ##contig=<ID=21>
    ##contig=<ID=22>
    ##contig=<ID=X>

["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
# ID = gnomAD ID or rsID
# chrom, pos, ref, alt = take from db
#qual, filter = .
##INFO=<ID=NS,Number=1,Type=Float,Description="Number of Samples With Data">
# Info = AF_afr={af_afr};AF_nfe={af_nfe},

# def get_line(var):
    

# Bei ID kannst Du gern die gnomAD ID nehmen, um da alles nachverfolgen zu können.


# Die anderen Spalten (QUAL, FILTER, INFO, FORMAT) können alle einen "." (Punkt) haben und dann passt das.

