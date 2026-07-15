import os

from src.AmiGo2.fetch_Amigo2 import download_data, get_all_genes_from_path

def get_config(config):
    storage = config.get("absolute_file_paths")
    data_path = os.path.join(storage.get("data"), "Amigo2", "genes")
    config_path = storage.get("config")
    os.makedirs(data_path, exist_ok=True)

    overview = config.get("relative_file_paths").get("AmiGo2_overview")
    overview_file = os.path.join(config_path, overview)

    include_exclude = config.get("relative_file_paths").get("Amigo2_inclue_exclude")
    include_exclude_file = os.path.join(config_path, include_exclude)

    download = config.get("flags").get("download_data")

    return (data_path, overview_file, include_exclude_file, download)
