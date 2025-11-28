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

def query_gen_ensemble(ensemble_id):
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
    }}
}} ''')

def query_variant_ensemble(ensemble_id):
    return (f''' query VariantsInGene {{
        gene(gene_id: "{ensemble_id}", reference_genome: GRCh38) {{
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
        }}
    }} ''')

def query_clinvar_ensemble(ensemble_id):
    return (f''' query VariantsInGene {{
        gene(gene_id: "{ensemble_id}", reference_genome: GRCh38) {{
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

def query_exons_ensemble(ensemble_id):
    return (f''' query VariantsInGene {{
        gene(gene_id: "{ensemble_id}", reference_genome: GRCh38) {{
            exons {{
                feature_type
                start
                stop
            }}
        }}
    }} ''')

def query_transcripts_ensemble(ensemble_id):
    return (f''' query VariantsInGene {{
        gene(gene_id: "{ensemble_id}", reference_genome: GRCh38) {{
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
        }}
    }} ''')


def query_from_ensemble_id(ensemble_id):
    """
    very long query string; 
    look up https://gnomad.broadinstitute.org/api to test out querys Strings
    there is also some aka little to none documentation
    # <- mark comments
    # too long, every now and then there are mistakes, so I split it in 5 querys
    """
    return [
        query_gen_ensemble(ensemble_id), 
        query_variant_ensemble(ensemble_id),
        query_clinvar_ensemble(ensemble_id),
        query_exons_ensemble(ensemble_id),
        query_transcripts_ensemble(ensemble_id)
    ]


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
    i=0
    for query in query_from_ensemble_id(ensemble_id):
        i+=1
        try:
            response = r.post("https://gnomad.broadinstitute.org/api",
                            json={"query": query}
                            )
            if response.ok:
                print(f"downloaded {ensemble_id} part {i}/5")

            yield response.json(), i
        except Exception as e:
            print(str(e))
            yield {}, i

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
        print(f"saved {file_name}")
        return True
    else:
        #print(f"data was empty, and so wasn't saved name was: {name if name else 'NO_NAME_GIVEN'}")
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
    # print(paths)
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

    return os.path.join(path_to_gnomAD, "data")

def check_for_offline(path_gen, ensemble_id):
    if not os.path.isdir(path_gen):
        return False

    look_for = {
        ensemble_id + names : 0
        for names in ["_gene_ids", "_gnomAD_variants", "_clinvar_variants", "_exons", "_transcripts"]
    }

    for entry in os.scandir(path_gen):
        if not entry.is_file():
            continue
        for key in look_for.keys():
            if key in entry.name:
                look_for[key] = 1

    return (sum(look_for.values()) == 5)

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
        print(f"getting {ensemble_id} data")
        path_gen = os.path.join(data_path, ensemble_id)
        if check_for_offline(path_gen, ensemble_id):
            paths.append(path_gen)
            continue

        check_result = 0
        for res, index in fetch_data_as_json(ensemble_id):
            gene_ids = res.get("data", {}).get("gene", {})

            if not gene_ids:
                tasks.append(ensemble_id)
                continue

            match index:
                case 1:
                    if save_json(gene_ids, path_gen, f"{ensemble_id}_gene_ids"):
                        check_result +=1
                case 2:
                    variants = gene_ids.get("variants", [])
                    if save_json(variants, path_gen, f"{ensemble_id}_gnomAD_variants"):
                        check_result +=1
                case 3:
                    clinvar_variants = gene_ids.get("clinvar_variants", [])
                    if save_json(clinvar_variants, path_gen, f"{ensemble_id}_clinvar_variants"):
                        check_result +=1
                case 4:
                    exons = gene_ids.get("exons", [])
                    if save_json(exons, path_gen, f"{ensemble_id}_exons"):
                        check_result +=1
                case 5:
                    transcripts = gene_ids.get("transcripts", [])
                    if save_json(transcripts, path_gen, f"{ensemble_id}_transcripts"):
                        check_result +=1
                case _:
                    raise("Something is wrong with the requests")

        if check_result < 5 and (path_gen not in paths):
            tasks.append(ensemble_id)

        if check_result == 5:
            print(f"succesfull downloaded {ensemble_id} data")

        if check_result > 0:
            paths.append(path_gen)

        end = 40 - (start - datetime.now()).total_seconds()

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


def filter_variants_by_ancestry(variant_path, ancestry_list, cut_off, gen_folder):
    if not variant_path:
        return [], [], "no_variant"

    variant_list = []
    for variant_json in variant_path:
        try:
            with open(variant_json, 'r') as f:
                variant_list += json.load(f)
        except Exception as e:
            print("filter_variants_by_ancestry1")
            print(str(e))
            continue

    variants_to_keep = []
    variants_sorted_out = set()
    check_for_doubles = []

    while variant_list:
        variant = variant_list.pop()
        if not isinstance(variant, dict):
            continue

        id = variant.get("variant_id", "NO_ID")
        if (id in check_for_doubles or 
            id == "NO_ID"):
            continue

        try:
            if not found(variant, ancestry_list, cut_off):
                variants_sorted_out.add(id)
                continue

            new = variant.copy()
            for item in ["joint", "genome", "exome"]:
                population_data = variant.get(item, None)
                new[item] = get_ancestry_p_and_reduce(population_data, ancestry_list)

            variants_to_keep.append( pd.json_normalize(new) )
            check_for_doubles.append(id)

        except Exception as e:
            print("filter_variants_by_ancestry2")
            print(str(e))
            continue

    if not check_for_doubles:
        return [], list(variants_sorted_out), "no_variant"

    df = pd.concat(variants_to_keep)
    if df.empty:
        return check_for_doubles, list(variants_sorted_out), "no_variant"

    anc_name = "".join( [f"_{anc}" for anc in ancestry_list] )
    file_name = os.path.join(gen_folder, f"variants{anc_name}_{cut_off}.tsv")
    try:
        df = pd.concat(variants_to_keep)
        if not df.empty:
            df.to_csv( file_name, sep='\t')
            print(gen_folder, " successfull")
            return check_for_doubles, list(variants_sorted_out), file_name
    except Exception as e:
        print("filter_variants_by_ancestry3")
        print(str(e))

    return check_for_doubles, list(variants_sorted_out), "no_variant"

def get_gnomAD_datapaths_single(path_to_gene):
    gnomAD_data = {
        "gene_ids" : [],
        "gnomAD_variants" : [],
        "clinvar_variants" : [],
        "exons" : [],
        "transcripts" : [] 
    }

    for entry in os.scandir(path_to_gene):
        if not entry.is_file():
            continue

        for name in gnomAD_data.keys():
            if name in entry.name:
                gnomAD_data[name].append(os.path.join(path_to_gene, entry.name))

    return gnomAD_data


def clean_and_filter(path_to_gene, ancestry_list, cutoff):
    print(path_to_gene, "starting now")
    gnomAD_data_dict = get_gnomAD_datapaths_single(path_to_gene)
    ancestry_list = extend_and_validate_ancestry_names(ancestry_list)

    variants_to_keep, variants_sorted_out, variant_tsv = filter_variants_by_ancestry(gnomAD_data_dict["gnomAD_variants"], ancestry_list, cutoff, path_to_gene)
    clinvar_found, clinvar_variants_sorted_out, clinvar_tsv = match_variants_to_clinvar(variants_to_keep, gnomAD_data_dict["clinvar_variants"], path_to_gene)

    variants_sorted_out = list(set(variants_sorted_out + clinvar_variants_sorted_out))

    file_txt = os.path.join(path_to_gene, "overview.txt")
    file_tsv = os.path.join(path_to_gene, "gene_information.tsv")
    with open (file_txt, "w") as f:
        f.write(f"> gene_information: {str(file_tsv)}\n")
        f.write(f"> variants_found: {str(variant_tsv)}\n")
        for var in variants_to_keep:
            f.write(f"{str(var)}, ")
        f.write(f"\n> clinvar_variants_found: {str(clinvar_tsv)}\n")
        for var in clinvar_found:
            f.write(f"{str(var)}, ")
        f.write("\n> variants_sorted_out\n")
        for var in variants_sorted_out:
            f.write(f"{str(var)}, ")

    return

def match_variants_to_clinvar(gnomAD_variants, path_to_clinvar_data, gen_folder):
    if not path_to_clinvar_data:
        return [], [], "no clinvar variant"

    clinvars = []
    for variant_json in path_to_clinvar_data:
        try:
            with open(variant_json, 'r') as f:
                clinvars += json.load(f)
        except Exception as e:
            print("match_variants_to_clinvar1")
            print(str(e))
            continue

    if not gnomAD_variants:
        return [], pd.json_normalize(clinvars)["variant_id"].unique().tolist(), None

    clinvar_variants_to_keep = []
    clinvar_variants_sorted_out = set()
    doubles = []

    while clinvars:
        variant = clinvars.pop()
        if not isinstance(variant, dict):
            continue

        var_id = variant.get("variant_id", "clinvar_var")

        if (var_id in doubles or 
            var_id == "clinvar_var"):
            continue

        if var_id not in gnomAD_variants:
            clinvar_variants_sorted_out.add(var_id)
            continue

        try:
            clinvar_variants_to_keep.append(pd.json_normalize(variant))
            doubles.append(var_id)
        except Exception as e:
            print("match_variants_to_clinvar2")
            print(str(e))
            continue

    if not doubles:
        return [], list(clinvar_variants_sorted_out), "no clinvar variant"

    file_name = os.path.join(gen_folder, f"matched_clinvar_variants.tsv")

    try:
        df = pd.concat(clinvar_variants_to_keep)
        if df.empty:
            return doubles, list(clinvar_variants_sorted_out), "no clinvar variant"
        df.to_csv( file_name, sep='\t')
        return doubles, list(clinvar_variants_sorted_out), file_name

    except Exception as e:
        print("match_variants_to_clinvar3")
        print(str(e))

    try:
        with open(file_name, 'w') as file:
            json.dump(clinvar_variants_to_keep, file)
        return doubles, list(clinvar_variants_sorted_out), file_name
    except Exception as e:
        print("match_variants_to_clinvar4")
        print(str(e))

    return doubles, list(clinvar_variants_sorted_out), "no clinvar variant"

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
      hgvs
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

