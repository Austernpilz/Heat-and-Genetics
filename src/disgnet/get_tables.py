import os
import pandas as pd

from src.helpers.folder_magic import search_for_files
from src.helpers.table_magic import load_include_exclude_txt


def get_config(config):
    storage = config.get("absolute_file_paths")
    data_path = os.path.join(storage.get("data"), "disgnet")
    config_path = storage.get("config")

    disgnet_data = config.get("relative_file_paths").get("disgnet_data")
    data_path = os.path.join(data_path, disgnet_data)
    os.makedirs(data_path, exist_ok=True)

    include_exclude = config.get("relative_file_paths").get("disgnet_include_exclude")
    include_exclude_file = os.path.join(config_path, include_exclude)

    download = config.get("flags").get("download_data")

    return (data_path, include_exclude_file, download)

def build_tables(disgnet_config):
    data_path, include_exclude_file, download = disgnet_config
    in_ex_group = load_include_exclude_txt(include_exclude_file)
    resulttsv = search_for_files(data_path, "search_result", "tsv")

    df_list = []
    for tsv in resulttsv:
        try:
            df = pd.read_csv(tsv, sep="\t", dtype=str)
            if df.empty:
                continue
            df_list.append(df)
        except Exception as e:
            print(str(e))
    df = pd.concat(df_list).drop_duplicates()

    df.rename(columns={"disease_name": "term"}, inplace=True)
    df["group_term"] = df["term"].map(lambda x : in_ex_group["group"].get(x, "NO GROUP"))
    df["geneEnsemblIDs"] = df["geneEnsemblIDs"].fillna("NO ID")
    return df


# empty functions, they are here to use the result_info_txt
# def look_for_extra_information(path_to_disgnet):
#     for path, _, files in os.walk(path_to_disgnet):
#         if "extend_data.txt" in files:
#             return os.path.join(path, "extend_data.txt")
#     return None

# def extend_data(df, path_to_disgnet):
#     extra = look_for_extra_information(path_to_disgnet)
#     if extra is None:
#         return df

#     return df