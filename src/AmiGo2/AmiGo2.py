from src.AmiGo2.mod_Amigo2 import get_table, build_amigo_tables
from src.helpers.std_out import send_message

def get_data(amigo_send, amigo_config, go_id_list=None):
    send_message("started", 0, "amigo")
    term_names = []
    if go_id_list is None:
        go_id_list = []
    else:
        term_names = [None] * len(go_id_list)
    overview_df, group_in_ex, data_path, result_path, download, extra = amigo_config
    #offline_data = get_all_genes_from_path(data_path)
    term_names += overview_df["Name"].tolist()
    go_id_list += overview_df["Accession"].tolist()
    send_message(len(go_id_list), 2, "amigo")

    for go_id, term in zip(go_id_list, term_names):
        if term in group_in_ex["exclude"]:
            continue
        if not extra and term in group_in_ex["extra"]:
            continue
        df = get_table(data_path, go_id, term, download)
        if df is None:
            send_message(f"{term} missing", 0, "amigo")
            continue

        name_symbol = df[["bioentity_name", "bioentity_label"]].drop_duplicates(ignore_index=True)
        amigo_send.put(name_symbol)

        term = df["term"].iat[0]
        send_message(1, 1, "amigo")
        send_message(f"got {term}", 0, "amigo")

    amigo_send.put("finished")
    send_message("finished", 0, "amigo")
    df = build_amigo_tables(data_path, result_path, group_in_ex, extra)
    if df is None:
        send_message("Final Amigo table not build")


