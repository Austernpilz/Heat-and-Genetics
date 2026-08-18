import os
import pandas as pd

from src.AmiGo2.fetch_Amigo2 import get_col, download_table_from_go_id

from src.helpers.folder_magic import search_for_files, search_for_file, check_string
from src.helpers.std_out import send_message


def get_single_df_from_path(table_path, columns = get_col() + ["term"]):
    if check_string(table_path):
        return None
    if not os.path.isfile(table_path):
        return None

    try:
        df = pd.read_csv(table_path, sep="\t", dtype=str)
        for col in columns:
            if col in df.columns:
                continue
            elif col == "term":
                term = os.path.basename(os.path.dirname(table_path))
                df[col] = term
            else:
                df[col] = ""
        new_df = df[columns]
        return new_df.drop_duplicates(ignore_index=True)

    except Exception as e:
        send_message(f" - coulnd't get table from path {table_path}\n{str(e)}\n")
        return None


def get_all_genes_from_path(amigo_data):
    files = search_for_files(amigo_data, "data", "tsv")
    columns = get_col() + ["term"]
    df_list = []
    for data_tsv in files:
        df = get_single_df_from_path(data_tsv, columns)
        if df is not None:
            df_list.append(df)

    if not df_list:
        send_message(" - couldn't load amigo tables: no tables found")
        return None

    return pd.concat(df_list, ignore_index=True).drop_duplicates(ignore_index =True)

def check_in_ex(term, group_in_ex):
    if term in group_in_ex["include"]:
        return "include"
    elif term in group_in_ex["exclude"]:
        return "exclude"
    elif term in group_in_ex["extra"]:
        return "extra"
    else:
        return "not_found"

def build_amigo_tables(data_path, result_path, group_in_ex, extra):
    df = get_all_genes_from_path(data_path)
    if df is None:
        return None
    df["group"] = df["term"].map(lambda x : group_in_ex["group"].get(x, "NO_GROUP"))
    df["in_ex"] = df["term"].map(lambda x : check_in_ex(x, group_in_ex))
    file_name = os.path.join(result_path, "amigo_df_complete.tsv")
    df.to_csv(file_name, sep="\t", index=False)

    if extra:
        df = df[df["term"].isin(group_in_ex["include"] + group_in_ex["extra"])]
    else:
        df = df[df["term"].isin(group_in_ex["include"])]

    file_name = os.path.join(result_path, "amigo_df_sub.tsv")
    df.to_csv(file_name, sep="\t", index=False)
    return df


def get_table(data_path, go_id, term_name=None, download=False):
    if check_string(term_name):
        term_name = f"GO_ID_{go_id}"
    dir_path = os.path.join(data_path, term_name)
    file = search_for_file(dir_path, "data", "tsv")
    if not file or download:
        file = download_table_from_go_id(data_path, go_id, term_name)
    return get_single_df_from_path(file)

