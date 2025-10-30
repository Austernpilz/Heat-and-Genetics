import requests as r
import asyncio #because the loading times are so long, i use parallel and async processes
import os
from datetime import datetime
import json
from collections import deque
import pandas as pd
import numpy as np
from time import sleep
import threading
import time

# broad_information
# this is used as the ground truth
# the key values should be corresponding to the ids in gnomAD
GNOMAD_POPULATION_NAMES = {
    'afr': 'African/African American',
    'ami': 'Amish',
    'amr': 'Admixed American',
    'asj': 'Ashkenazi Jewish',
    'eas': 'East Asian',
    'mid': 'Middle Eastern',
    'eur': 'European',
    'nfe': 'European (non-Finnish)',
    'fin': 'European (Finnish)',
    'oth': 'Remaining individuals',
    'sas': 'South Asian',
    'rmi': 'Remaining',
    'remaining': 'Remaining',

    # EAS subpopulations
    'eas_jpn': 'Japanese',
    'eas_kor': 'Korean',
    'eas_oea': 'Other East Asian',

    # NFE subpopulations
    'nfe_bgr': 'Bulgarian',
    'nfe_est': 'Estonian',
    'nfe_nwe': 'North-western European',
    'nfe_onf': 'Other non-Finnish European',
    'nfe_seu': 'Southern European',
    'nfe_swe': 'Swedish',
}

def query_from_ensemble_id(ensemble_id):
    """
    very long query string; 
    look up https://gnomad.broadinstitute.org/api to test out querys Strings
    there is also some aka little to none documentation
    # <- mark comments
    """
    return (f''' 
query VariantsInGene {{
    # this is what we are looking for, in ref_g: GRCh38 look for ensemble_id
    gene(gene_id: "{ensemble_id}", reference_genome: GRCh38) {{
    # these are all the datafields we want to get back
    # some will give back lists
        reference_genome
        gene_id
        gene_version
        symbol
        gencode_symbol
        hgnc_id
        ncbi_id
        omim_id
        name
        chrom
        start
        stop
        strand
        canonical_transcript_id

        mane_select_transcript {{
            ensembl_id
            ensembl_version
            refseq_id
            refseq_version
        }}

        variants(dataset: gnomad_r4) {{
            variant_id
            reference_genome
            chrom
            pos
            ref
            alt
            rsids
            gene_id
            gene_symbol
            transcript_id
            transcript_version
            lof
            hgvsc
            hgvsp

            exome {{
                ac
                an

                populations {{
                    id
                    ac
                    an
                    homozygote_count
                    hemizygote_count
                }}
            }}

            genome {{
                ac
                an

                populations {{
                id
                ac
                an
                homozygote_count
                hemizygote_count
                }}
            }}

            joint {{
                ac
                an

                populations {{
                    id
                    ac
                    an
                    homozygote_count
                    hemizygote_count
                }}
            }}

            lof_curation {{
                verdict
                flags
            }}

            in_silico_predictors {{
                id
                value
                flags
            }}
        }}

        exons {{
            feature_type
            start
            stop
        }}

        transcripts {{
            transcript_id
            start
            stop
            exons {{
                feature_type
                start
                stop
      	    }}
        }}

        clinvar_variants {{
            variant_id
            reference_genome
            chrom
            pos
            ref
            alt
            clinical_significance
            clinvar_variation_id
            gold_stars
            hgvsc
            hgvsp
            major_consequence
            review_status
            transcript_id
        }}
    }}
}} ''')


def query_by_region(chromosom, start, stop):
    """
    creating a very long, f string to use as a function argument,
    maybe there could be a better way, 
    but this is straigt out of the documentation how a query in python should look like
    and this was the easiest to make a function out of it
    """
    return (f'''
query VariantsInGene {{

    region(reference_genome: GRCh38, chrom: "{chromosom}", start: {start}, stop: {stop}) {{

        variants(dataset: gnomad_r4) {{
            gene_id
            variant_id
            gene_symbol
            reference_genome
            lof
            domains
            alt
            hgvsc
            hgvs
            hgvsp

            exome {{
                ac
                an

                populations {{
                    id
                    ac
                    an
                    homozygote_count
                    hemizygote_count
                    ac_hemi
                    ac_hom
                    }}
                }}

            genome {{
                ac
                an

                populations {{
                    id
                    ac
                    an
                    homozygote_count
                    hemizygote_count
                    ac_hemi
                    ac_hom
                }}
            }}

            joint {{
                ac
                an

                populations {{
                    id
                    ac
                    an
                    homozygote_count
                    hemizygote_count
                    ac_hemi
                    ac_hom
                }}
            }}
        }}
    }}
}} ''')


#unique names in case we don't get an gene id or symbol
def get_unique_name(path_to_data, name):
    ts = datetime.now().strftime("%Y%m%dT%H%M%SZ")
    return os.path.join(path_to_data, f"{name}_{ts}.json")

"""
to overload the function, we use the same thing and let the python interpreter choose
"""

def fetch_data_as_json(ensemble_id : str):
    try:
        response = r.post("https://gnomad.broadinstitute.org/api",
                        json={"query": query_from_ensemble_id(ensemble_id)}
                        )
        return response.json()
    except Exception as e:
        print(str(e))
        return {}

# def fetch_data_as_json(chromosome : str, pos_start : int, pos_end : int):
#     #chromosome = f'\"{chromosome}\"'
#     try:
#         response = r.post("https://gnomad.broadinstitute.org/api",
#                         json={"query": query_from(chromosome, pos_start, pos_end)}
#                         )
#         return response.json()
#     except Exception as e:
#         print(str(e))
#         return {}


def save_json(data, path_to_gnomAD_data, name=False):
    if data:
        os.makedirs(path_to_gnomAD_data, exist_ok=True)
        file_name = get_unique_name(path_to_gnomAD_data, name) if name else get_unique_name(path_to_gnomAD_data, "data")
        with open(file_name, 'w') as file:
            json.dump(data, file)
        return True
    else:
        print(f"data was empty, and so wasn't saved name was: {name if name else 'NO_NAME_GIVEN'}")
        return False

def get_paths(path_to_gnomAD, ENSG_ids):
    return_paths = []

    dir_walk = deque()
    dir_walk.append(path_to_gnomAD)

    while dir_walk:
        dir = dir_walk.pop()
        for entry in os.scandir(dir):
            if not entry.is_dir():
                continue

            if entry.name in ["bin", "lib", "include", "share", "HGNC", "AmiGo2", "disgnet", "ensembl", "figures"]:
                continue

            if entry.name in ENSG_ids:
                return_paths.append(os.path.join(dir, entry.name))
            else:
                dir_walk.append(os.path.join(dir, entry.name))

    return return_paths

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def get_data(ENSG_ids, path_to_gnomAD, ancestry_list, cutoff, t=1,  download=True):
    paths = []

    if download:
        paths = download_data(ENSG_ids, path_to_gnomAD)
    else: 
        paths = get_paths(path_to_gnomAD, ENSG_ids)
    print(paths)
    for chunk in chunks(paths, t):
        threads = []
        for ensg_id in chunk:
            # Using `args` to pass positional arguments and `kwargs` for keyword arguments
            t = threading.Thread(target=clean_and_filter, 
                                 kwargs={
                                     "path_to_gene": ensg_id,
                                     "ancestry_list": ancestry_list,
                                     "cutoff" : cutoff
                                     })
            threads.append(t)

        # Start each thread
        for t in threads:
            t.start()

        # Wait for all threads to finish
        for t in threads:
            t.join()


def download_data(ENSG_ids, path_to_gnomAD):
    paths = []
    tasks = deque()
    data_path = os.path.join(path_to_gnomAD, "data")
    os.makedirs(data_path, exist_ok=True)
    for ensemble_id in ENSG_ids:
        tasks.append(ensemble_id)

    while tasks:
        start = datetime.now()
        ensemble_id = tasks.pop()
        path_gen = os.path.join(data_path, ensemble_id)

        res = fetch_data_as_json(ensemble_id)
        gene_ids = res.get("data", {}).get("gene", {})
        variants = gene_ids.pop("variants", [])
        clinvar_variants = gene_ids.pop("clinvar_variants", [])
        exons = gene_ids.pop("exons", [])
        transcripts = gene_ids.pop("transcripts", [])

        check_result = 0

        if save_json(gene_ids, path_gen, f"{ensemble_id}_gene_ids"):
            check_result +=1

        if save_json(variants, path_gen, f"{ensemble_id}_gnomAD_variants"):
            check_result +=1

        if save_json(clinvar_variants, path_gen, f"{ensemble_id}_clinvar_variants"):
            check_result +=1

        if save_json(exons, path_gen, f"{ensemble_id}_exons"):
            check_result +=1

        if save_json(transcripts, path_gen, f"{ensemble_id}_transcripts"):
            check_result +=1

        if check_result < 5 or (path_gen in paths):
            tasks.append(ensemble_id)

        paths.append(path_gen)

        end = 6.0 - (start - datetime.now()).total_seconds()

        if end > 0:
            sleep(end)

    return paths


def found(variant, ancestry, cutoff):
    for item in ["joint", "genome", "exome"]:
        populations_an_ac = variant.get(item, None)
        if populations_an_ac is None:
            continue

        ac = populations_an_ac.get("ac", 0)
        if not ac:
            continue

        populations = populations_an_ac.get("populations", [])
        if not populations:
            continue

        for population in populations:
            pop_id = population.get("id", None)
            if pop_id is None or pop_id not in ancestry:
                continue

            ac = population.get("ac", 0)
            if not ac:
                continue

            an = population.get("an", 0)
            if not an:
                continue

            try:
                af = float(ac) / float(an)
                # the variant exists in our population with the specific cutoff
                if af >= cutoff:
                    return True
            except Exception as _:
                continue

    return False

def extend_and_validate_ancestry_names(ancestry_list):
    new_ancestry_list = set()
    global GNOMAD_POPULATION_NAMES
    for short_id, full_name in GNOMAD_POPULATION_NAMES.items():
        for anc in ancestry_list:
            if anc in full_name:
                new_ancestry_list.add(short_id)
        if short_id in ancestry_list:
            new_ancestry_list.add(short_id)
        if full_name in ancestry_list:
            new_ancestry_list.add(short_id)
    return list(new_ancestry_list)


def get_an_ac_af(some_dict):
    if not isinstance(some_dict, dict):
        return None, None, None

    ac = some_dict.get("ac", None)
    an = some_dict.get("an", None)
    af = None

    try:
        af = float(ac) / float(an)
    except Exception as _:
        af = None

    return ac, an, af

def get_ancestry_p_and_reduce(populations_in_gen, ancestry_list):
    new_ac_an_af = {
        "ac" : None,
        "an" : None,
        "af" : None
    }

    for anc in ancestry_list:
        new_ac_an_af[f"ac_{anc}"] = None
        new_ac_an_af[f"an_{anc}"] = None
        new_ac_an_af[f"af_{anc}"] = None

    if not isinstance(populations_in_gen, dict):
        return new_ac_an_af

    ac, an, af = get_an_ac_af(populations_in_gen)
    new_ac_an_af["ac"] = ac
    new_ac_an_af["an"] = an
    new_ac_an_af["af"] = af

    population_list = populations_in_gen.get("populations", [])
    for population_dict in population_list:
        anc_id = population_dict.get("id", "totaly not a name or population")
        if anc_id not in ancestry_list:
            continue

        ac, an, af = get_an_ac_af(population_dict)
        new_ac_an_af[f"ac_{anc_id}"] = ac
        new_ac_an_af[f"an_{anc_id}"] = an
        new_ac_an_af[f"af_{anc_id}"] = af

    return new_ac_an_af


def filter_variants_by_ancestry(variant_path, ancestry_list, cut_off):
    if variant_path is None:
        return [], []

    variant_list = []
    with open(variant_path, 'r') as f:
        variant_list = json.load(f)

    variants_to_keep = []
    variants_sorted_out = []
    ancestry_list = extend_and_validate_ancestry_names(ancestry_list)
    while variant_list:
        variant = variant_list.pop()
        if not isinstance(variant, dict):
            continue

        if not found(variant, ancestry_list, cut_off):
            variants_sorted_out.append(variant.get("variant_id", "NO_ID"))
            continue

        for item in ["joint", "genome", "exome"]:
            population_data = variant.pop(item, None)
            variant[item] = get_ancestry_p_and_reduce(population_data, ancestry_list)

        variants_to_keep.append(variant)

    return variants_to_keep, variants_sorted_out


def get_gnomAD_datapaths_single(path_to_gene):
    gnomAD_data = {
        "gene_ids" : None,
        "gnomAD_variants" : None,
        "clinvar_variants" : None,
        "exons" : None,
        "transcripts" : None 
    }

    for entry in os.scandir(path_to_gene):
        if not entry.is_file():
            continue

        for names in gnomAD_data.keys():
            if names in entry.name:
                gnomAD_data[names] = os.path.join(path_to_gene, entry.name)

    return gnomAD_data


def clean_and_filter(path_to_gene, ancestry_list, cutoff):
    gnomAD_data_dict = get_gnomAD_datapaths_single(path_to_gene)

    variant_table, variants_sorted_out = filter_variants_by_ancestry(gnomAD_data_dict["gnomAD_variants"], ancestry_list, cutoff)

    variant_table, clinvar_variants_sorted_out = match_variants_to_clinvar(variant_table, gnomAD_data_dict["clinvar_variants"])
    variants_sorted_out = list(set(variants_sorted_out + clinvar_variants_sorted_out))

    data = {}
    if gnomAD_data_dict["gene_ids"] is not None:
        with open(gnomAD_data_dict["gene_ids"], 'r') as f:
            data = json.load(f)
            data.pop("mane_select_transcript")
    
    for key in data.keys():
        if key in ["reference_genome", "chrom", "gene_id"]:
            variant_table[f"{key}_general"] = data[key]
        else:
            variant_table[key] = data[key]

    specifier = "".join( [f"{anc}_" for anc in ancestry_list] + [f"cutoff_{cutoff}.tsv"] )
    filename = os.path.join(path_to_gene, specifier)
    variant_table.to_csv(filename, index=False, sep="\t")


def build_entry(gv, cv):
    new_entry = {
        "variant_id" : None,
        "reference_genome" : None,
        "chrom" : None,
        "pos" : None,
        "ref" : None,
        "alt" : None,
        "transcript_id" : None,
        "hgvsc" : None,
        "hgvsp" : None,
        "no_match" : False
    }

    for item in new_entry.keys():
        gnomad_variation = gv.get(item, None)
        clinvar_variation = cv.get(item, None)
        if gnomad_variation == clinvar_variation:
            new_entry[item] = gnomad_variation
            continue

        if (gnomad_variation is None or 
            clinvar_variation and not gnomad_variation):
            new_entry[item] = clinvar_variation
            continue

        if (clinvar_variation is None or
            gnomad_variation and not clinvar_variation):
            new_entry[item] = gnomad_variation
            continue

        new_entry["no_match"] = True
        break

    if new_entry["no_match"]:
        return gv

    for key in gv:
        if key in new_entry.keys():
            continue
        new_entry[key] = gv[key]

    for key in cv.keys():
        if key in new_entry.keys():
            continue
        new_entry[key] = cv[key]

    return new_entry


def match_variants_to_clinvar(variants, path_to_clinvar_data):
    clinvars = []
    if path_to_clinvar_data is not None:
        with open(path_to_clinvar_data, 'r') as f:
            clinvars = json.load(f)

    if not variants and not clinvars:
        return pd.DataFrame(), []

    if not clinvars:
        return pd.json_normalize(variants), []

    if not variants:
        return pd.DataFrame(), pd.json_normalize(clinvars)["variant_id"].unique().tolist()

    full_variant_table = []
    clinvar_variants_sorted_out = set()
    for clinvar_variant in clinvars:
        var_id = clinvar_variant.get("variant_id", "clinvar_var")
        not_found = True

        for gnomAD_variant in variants:
            if var_id != gnomAD_variant.get("variant_id", "gnomAD_var"):
                continue

            not_found = False
            full_variant_table.append( build_entry(clinvar_variant, gnomAD_variant) )
            break

        if not_found:
            clinvar_variants_sorted_out.add(var_id)

    return pd.json_normalize(variants), list(clinvar_variants_sorted_out)



# def big_loop(big_dict):
#     return_dict = {}
#     found, relevant = 0, 0
#     for gene_id, data in big_dict.items():
#         variants = data.get("variants", [])
#         if not variants:
#             continue

#         d, f =  get_relevant_variants(variants)
#         found += f
#         return_dict[gene_id] = d
#         if not return_dict[gene_id]:
#             return_dict.pop(gene_id)
#             continue

#         clinvar = data.get("clinvar_variants", [])
#         if not clinvar:
#             continue

#         for variant_id, values in return_dict[gene_id]:
#             for cvar in clinvar:
#                 id = cvar.get("variant_id", "None")
#                 if id == "None" or id != variant_id:
#                     continue
#                 else:
#                     for item in ["reference_genome", "chrom", "pos", "ref", "alt", "clinical_significance", "clinvar_variation_id", "gold_stars", "hgvsc", "hgvsp", "major_consequence"]:
#                         key = cvar.get(item, "None")
#                         if key == "None":
#                             continue
#                         if item in return_dict[gene_id][variant_id]["info"].keys():
#                             if key == return_dict[gene_id][variant_id]["info"][item]:
#                                 continue
#                             else:
#                                 key2 = return_dict[gene_id][variant_id]["info"][item]
#                                 return_dict[gene_id][variant_id]["info"][item] = [key, key2]
#                         else:
#                             return_dict[gene_id][variant_id]["info"][item] = key
#     # for gene_id, data in big_dict.items():
#     #     print(gene_id, len(data))

#     # print(len(return_dict), f, r)
#     return return_dict

# anzahl varianten pro gen, 
# länge des gen, anzahl der variant
# table VEP
# get VEP
# 

# def clean(smaller_dict):
#     gene_id, gene_symbol, variant_id = [], [], []
#     clinical_significance, major_consequence = [], []
#     chrom, pos, ref, alt = [], [], [], []
#     hgvsc, hgvsp, reference_genome = [], [], []
#     ancestry_afr, ancestry_nfe = [], []

#     for g_id, gnomAD_data in smaller_dict:
#         for var_id, fields in gnomAD_data:
#             gene_id.append(g_id)
#             variant_id.append(var_id)

#             anc = fields.get("ancestry", [])
#             afr, nfe = 0.0, 0.0
 
#             for name, pop_list in anc:
#                 if name != "exome":
#                     continue
#                 for pop, af in pop_list:
#                     if pop in ['nfe', 'European (non-Finnish)', 
#                                'nfe_bgr', 'Bulgarian', 
#                                'nfe_est', 'Estonian', 
#                                'nfe_nwe', 'North-western European', 
#                                'nfe_onf', 'Other non-Finnish European', 
#                                'nfe_seu', 'Southern European', 
#                                'nfe_swe', 'Swedish']:
#                         nfe = af
#                         continue
#                     elif pop in ['afr', 'African/African American']: 
#                         afr = af
#                         continue
#             ancestry_afr.append(afr)
#             ancestry_nfe.append(nfe)

#             c, p, r, a = 0, 0, "", ""
#             rg, gs = "", ""
#             gs, cs, mc, hc, hp = "", "", "", "", ""
#             information = fields.get("info", {})
#             for key, value in information.items():
#                 if isinstance(value, list):
#                     value = value[0]

#                 if key in ["reference_genome", "reference_genome"]:
#                     rg = value
#                     continue

#                 if key in ["chrom", "chr"]:
#                     c = value
#                     continue
#                 if key == "pos":
#                     p = value
#                     continue

#                 if key == "ref":
#                     r = value
#                     continue

#                 if key =="alt":
#                     a = value
#                     continue

#                 if key == "gene_symbol":
#                     gs = value
#                     continue

#                 if key == "clinical_significance":
#                     cs = value
#                     continue

#                 if key == "major_consequence":
#                     mc = value
#                     continue

#                 if key == "hgvsc":
#                     hc = value
#                     continue

#                 if key == "hgvsp":
#                     hp = value
#             gene_symbol.append(gs)
#             clinical_significance.append(cs)
#             major_consequence.append(mc)
#             chrom.append(c)
#             pos.append(p)
#             ref.append(r)
#             alt.append(a)
#             hgvsc.append(hc)
#             hgvsp.append(hp)
#             reference_genome.append(rg)


#     return pd.DataFrame({
#         "gene_id":gene_id, 
#         "gene_symbol":gene_symbol, 
#         "variant_id":variant_id, 
#         "clinical_significance":clinical_significance, 
#         "major_consequence":major_consequence ,
#         "chr":chrom, 
#         "pos":pos, 
#         "ref": ref,
#         "alt" : alt,
#         "hgvsc": hgvsc, 
#         "hgvsp": hgvsp,
#         "reference_genome": reference_genome,
#         "ancestry_afr" : ancestry_afr,
#         "ancestry_nfe" : ancestry_nfe
#     })

# def get_data(path_to_gnomAD):
#     return_dict = {}
#     path_to_gnomAD = os.path.join(path_to_gnomAD, "data")
#     if not os.path.isdir(path_to_gnomAD):
#         return {}

#     dir_to_visit = deque()
#     dir_to_visit.append(path_to_gnomAD)
#     while dir_to_visit:
#         current = dir_to_visit.pop()

#         for entry in os.scandir(current):
#             if entry.name in ["bin", "include", "lib", "overview.txt", "data.tsv", "include_exclude.txt"]:
#                 continue

#             if entry.is_dir():
#                 dir_to_visit.append(os.path.join(current, entry.name))
#                 continue

#             if not entry.is_file():
#                 continue

#             if not entry.name.endswith(".json"):
#                 continue

#             try:
#                 data = {}
#                 gene_name = os.path.basename(current)
#                 file_path = os.path.join(current, entry.name)

#                 if gene_name not in return_dict:
#                     return_dict[gene_name] = {}

#                 with open(file_path, 'r') as f:
#                     data = json.load(f)

#                 if "gene_ids" in entry.name:
#                     return_dict[gene_name]["gene_ids"] = data
#                 elif "clinvar_variants" in entry.name:
#                     return_dict[gene_name]["clinvar_variants"] = data
#                 elif "variants" in entry.name:
#                     return_dict[gene_name]["variants"] = data
#                 else:
#                     return_dict[gene_name]["NO_NAME"] = data

#             except Exception as _:
#                 continue

#     return return_dict

# def download_data(df: pd.DataFrame, path_to_gnomAD):
#     result_dict = {}
#     tasks = deque()

#     for name, chromosom, start, end in zip(df['genes'].tolist(), df['chr'].tolist(), df['start'].tolist(), df['end'].tolist()):
#         tasks.append((name, chromosom, start, end))

#     data_path = os.path.join(path_to_gnomAD, "data")
#     os.makedirs(data_path, exist_ok=True)

#     while tasks:
#         start = datetime.now()
#         gen_symbol, c, b, e = tasks.pop()
#         path_gen = os.path.join(data_path, gen_symbol)
#         try:
#             res = fetch_data_as_json(c,b,e)
#             variants = res.get("data", {}).get("region", {}).get("variants", [])
#             genes = res.get("data", {}).get("region", {}).get("genes", [])

#             if gen_symbol not in result_dict:
#                 result_dict[gen_symbol] = {}

#             if save_json(variants, path_gen, f"{gen_symbol}_variants"):
#                 result_dict[gen_symbol]["variants"] = variants

#             if save_json(genes, path_gen, f"{gen_symbol}_genes"):
#                 result_dict[gen_symbol]["genes"] = genes

#         except Exception as e:
#             print(str(e))
#             if gen_symbol not in result_dict:
#                 #await generate_task(n,c,b,e, tasks)
#                 result_dict[gen_symbol] = {}
#         # if datetime.now() - start:
#     return result_dict

# def test():
#     gnomAD = os.path.join(os.getcwd(), "gnomAD")

#     # df = pd.DataFrame({
#     #     "genes" : ["ENST00000644486", "ATXN3", "PPP5C"],
#     #     "chr" : ["14", "14","19"],
#     #     "start": [92058552, 92044496, 46347108, ],
#     #     "end": [92106582, 92106622, 46390852],
#     #     })
#     df = ["ENSG00000008869", "ENSG00000066427", "ENSG00000011485"]
#     download_data(df, gnomAD)

#never use print with gnomAD data :'(
#test()

"""
region(reference_genome: GRCh38, chrom: "19", start: 46347108, stop: 46390852) {
    genes {
      gene_id
      symbol
      start
      stop
      exons {
        feature_type
        start
        stop
      }
      transcripts {
        transcript_id
        start
        stop
        exons {
          feature_type
          start
          stop
      	}
      }
    }
    clinvar_variants {
      variant_id
      reference_genome
      chrom
      pos
      ref
      alt
      clinical_significance
      clinvar_variation_id
      gold_stars
      hgvsc
      hgvsp
      major_consequence
      review_status
      transcript_id
    }
    variants(dataset: gnomad_r4) {
     variant_id
     reference_genome
     chrom
     pos
     ref
     alt
     rsids
     gene_id
     gene_symbol
     transcript_id
     transcript_version
     lof
     hgvsc
     hgvsp
     exome {
       ac
       an
        populations {
         id
         ac
          an
         homozygote_count
         hemizygote_count
       }
     }
     genome {
       ac
       an
       populations {
         id
         ac
          an
         homozygote_count
         hemizygote_count
       }
     }
     joint {
       ac
       an
       populations {
         id
         ac
          an
         homozygote_count
         hemizygote_count
       }
     }
      lof_curation {
       verdict
       flags
      }
      in_silico_predictors {
       id
       value
       flags
     }
   }
  }
}
"""


"""
query VariantsInGene {
  gene(gene_id: "ENSG00000008869", reference_genome: GRCh38) {
    reference_genome
    gene_id
    gene_version
    symbol
    gencode_symbol
    hgnc_id
    ncbi_id
    omim_id
    name
    chrom
    start
    stop
    strand
    canonical_transcript_id
    mane_select_transcript {
      ensembl_id
      ensembl_version
      refseq_id
      refseq_version
    }
    variants(dataset: gnomad_r4) {
      variant_id
      reference_genome
      chrom
      pos
      ref
      alt
      rsids
      gene_id
      gene_symbol
      transcript_id
      transcript_version
      lof
      hgvsc
      hgvsp
      exome {
        ac
        an
        populations {
          id
          ac
          an
          homozygote_count
          hemizygote_count
        }
      }
      genome {
        ac
        an
        populations {
          id
          ac
          an
          homozygote_count
          hemizygote_count
        }
      }
      joint {
        ac
        an
        populations {
          id
          ac
          an
          homozygote_count
          hemizygote_count
        }
      }
      lof_curation {
        verdict
        flags
      }
      in_silico_predictors {
        id
        value
        flags
      }
    }
    exons {
      feature_type
      start
      stop
    }
    transcripts {
      transcript_id
      start
      stop
      exons {
        feature_type
        start
        stop
      }
    }
  }
}
"""