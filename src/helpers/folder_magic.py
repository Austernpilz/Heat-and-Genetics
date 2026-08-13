
import os
import numpy as np


def check_string(some_string):
    if some_string is None:
        return True
    if not isinstance(some_string, str):
        return True
    if some_string in ["", "nan", "Nan", "NaN", "NAN", "empty", " ", "_", "None", "-", np.nan, "NO_GROUP", "NO_SYMBOL", "NO_NAME", "NO_ID", "time_out"]:
        return True
    return False

def save_table(df, some_path, name):
    os.makedirs(some_path, exist_ok=True)
    file_path = os.path.join(some_path, f"{name}.tsv")
    i=0
    while os.path.isfile(file_path):
        file_path = os.path.join(some_path, f"{name}_{i}.tsv")
        i+=1
    df.to_csv(file_path, index=False, sep="\t")


def search_for_dir(some_path, some_string):
    if not os.path.isdir(some_path):
        return False

    dir_to_visit = [some_path]
    while dir_to_visit:
        current = dir_to_visit.pop()
        try: #in case i get a permission denied
            for entry in os.scandir(current):
                if entry.is_dir():
                    dir_path = os.path.join(current, entry.name)
                    if some_string in entry.name:
                        return dir_path
                    else:
                        dir_to_visit.append(dir_path)

        except Exception as _:
            continue

    return False


def search_for_dirs(some_path, some_string):
    if not os.path.isdir(some_path):
        return []

    dir_to_visit = [some_path]
    found_something = []
    while dir_to_visit:
        current = dir_to_visit.pop()
        try: #in case i get a permission denied
            for entry in os.scandir(current):
                if entry.is_dir():
                    dir_path = os.path.join(current, entry.name)
                    dir_to_visit.append(dir_path)
                    if some_string in entry.name:
                        found_something.append(dir_path)

        except Exception as _:
            continue

    return found_something


def search_for_file(some_path, some_string, some_suffix):
    if os.path.isfile(some_path):
        file_name = os.path.basename(some_path)
        if some_string in file_name and file_name.endswith(some_suffix):
            return some_path

    if not os.path.isdir(some_path):
        return None

    dir_to_visit = [some_path]
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
                    return os.path.join(current, entry.name)

        except Exception as _:
            continue

    return None


def search_for_files(some_path, some_string, some_suffix):
    if os.path.isfile(some_path):
        file_name = os.path.basename(some_path)
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