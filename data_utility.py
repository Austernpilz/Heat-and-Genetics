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
    return new_df[new_df>cutoff].sort_values(ascending=ascend)


def top_cut(df, cutoff, col):
    top_x = df[col].value_counts().index.tolist()[:cutoff]
    return df[df[col].isin(top_x)]


def new_table(df, col_list, rename_dict):
    filtered_list = [col for col in col_list if col in df.columns]
    return df[col_list].rename(columns=rename_dict)

def make_new_table(df_list, cols_list, rename_dict, axis=0):
    # axis=0 means I assume you condense your dict, 
    # and want to merge different columns into one
    # so it will concat on index (under the other data)

    # this is not efficient for larger datasets but quick to implement
    # it will copy the dataframe twice, once for new_table and once for concat
    return_list = []
    for df in df_list:
        return_list.append(new_table(df, cols_list, rename_dict))
    return pd.concat(return_list, axis)

