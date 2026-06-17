import os
import json



def search_for_dir(some_path, some_string):
    if not os.path.isdir(some_path):
        return []

    dir_to_visit = [some_path]
    found_something = []
    while dir_to_visit:
        current = dir_to_visit.pop()
        try: #in case i get a permission denied
            for entry in os.scandir(current):
                if entry.is_file():
                    continue
                elif entry.is_dir():
                    dir_path = os.path.join(current, entry.name)
                    dir_to_visit.append(dir_path)
                    if some_string in entry.name:
                        found_something.append(dir_path)

        except Exception as _:
            continue

    return found_something


def search_for_file(some_path, some_string, some_suffix):
    if os.path.isfile(some_path):
        file_name = os.path.basename(some_string)
        if some_string in file_name and file_name.endswith(some_suffix):
            return [some_path]

    if not os.path.isdir(some_path):
        return []

    dir_to_visit = [some_path]
    found_something = []
    while dir_to_visit:
        current = dir_to_visit.pop()
        try: #in case i get a permission denied
            for entry in os.scandir(current):
                if entry.is_dir():
                    dir_to_visit.append(os.path.join(current, entry.name))
                elif (
                    entry.is_file() and 
                    some_string in entry.name and 
                    entry.name.endswith(some_suffix)
                ):
                    found_something.append(os.path.join(current, entry.name))

        except Exception as _:
            continue

    return found_something


def get_config(pathpath = False):
    # this excpects you to start from the Heat_and_Genetics_Folder 
    # or give the config path as argument
    # parser = argparse.ArgumentParser(
    #                 prog='heat_and_genetics',
    #                 description='pipeline for downloading and plotting gene data for heatrelated healthriscs')
    # parser.add_argument('--folder_path', '-p', default=None)
    # inputs = parser.parse_args()
    if pathpath:
        config = {}
        try:
            with open(pathpath, 'r') as f:
                config = json.load(f)
            return config
        except Exception as e:
            print(str(e))
            #no config, try to build it myself?!?
    this_folder = os.getcwd()
    config = os.path.join(this_folder, "config")
    data = os.path.join(this_folder, "data")
    return {
        "flags": {
            "download_data": False,
            "force_data_update": False,
            "threads": 10
        },
        "relative_file_paths": {
            "AmiGo2_overview": search_for_file(config, "AmiGo2_overview", "txt"),
            "Amigo2_inclue_exclude": search_for_file(config, "AmiGo2_include_exclude", "txt"),
            "disgnet_data": search_for_dir(data, "disgnet"),
            "disgnet_include_exclude": search_for_file(config, "disgnet_include_exclude", "txt"),
            "hgnc_symbol_check": search_for_file(data, "hgnc-symbol-check", "csv"),
            "data_storage": data
        },
        "populations": {
            "afr": "African/African American",
            "nfe": "European (non-Finnish)"
        }
    }