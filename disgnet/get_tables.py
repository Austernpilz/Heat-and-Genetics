import os
import pandas as pd


def get_paths_to_disgnet_data(path_disgnet = None):
    path_disgnet = os.path.join(os.getcwd(), "disgnet") if path_disgnet is None else path_disgnet

    if not os.path.isdir(path_disgnet):
        return []

    #infotxt = []
    resulttsv = []
    dir_to_visit = [path_disgnet]
    while dir_to_visit:
        current = dir_to_visit.pop(0)
        try:
            for entry in os.scandir(current):
                if entry.name in ["bin", "include", "lib", "overview.txt", "data.tsv", "include_exclude.txt"]:
                    continue
                elif entry.is_dir():
                    dir_to_visit.append(os.path.join(current, entry.name))
                    continue
                elif not entry.is_file():
                    continue

                # if entry.name.endswith(".txt"):
                #     infotxt.append(os.path.join(current, entry.name))
                if entry.name.endswith(".tsv"):
                    resulttsv.append(os.path.join(current, entry.name))

        except Exception as _:
            continue

    return resulttsv #infotxt, 

# print(get_paths_to_disgnet_data(os.getcwd()))

def build_tables(path_to_disgnet = None):
    resulttsv = get_paths_to_disgnet_data(path_to_disgnet)

    df_list = []
    for tsv in resulttsv:
        df_list.append(
            pd.read_csv(tsv, 
                        sep="\t", 
                        dtype=str))

    return pd.concat(df_list)

# print(build_tables(os.getcwd()))