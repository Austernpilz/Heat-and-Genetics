
from datetime import datetime
import os
import json
from time import sleep
import requests as r
import numpy as np

from src.helpers.folder_magic import search_for_files, search_for_file, check_string
from src.helpers.std_out import send_message


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

def extend_and_validate_ancestry_names(ancestry_list):
    new_ancestry_list = set()
    global GNOMAD_POPULATION_NAMES
    for short_id, full_name in GNOMAD_POPULATION_NAMES.items():
        if short_id in ancestry_list:
            new_ancestry_list.add(short_id)
        if full_name in ancestry_list:
            new_ancestry_list.add(short_id)
    return new_ancestry_list


def query_gen_ensemble(ensembl_id):
    return (f''' 
query VariantsInGene {{
    # this is what we are looking for, in ref_g: GRCh38 look for ensembl_id
    gene(gene_id: "{ensembl_id}", reference_genome: GRCh38) {{
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

def query_variant_ensemble(ensembl_id):
    return (f''' query VariantsInGene {{
        gene(gene_id: "{ensembl_id}", reference_genome: GRCh38) {{
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
                hgvs
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

def query_clinvar_ensemble(ensembl_id):
    return (f''' query VariantsInGene {{
        gene(gene_id: "{ensembl_id}", reference_genome: GRCh38) {{
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

def query_exons_ensemble(ensembl_id):
    return (f''' query VariantsInGene {{
        gene(gene_id: "{ensembl_id}", reference_genome: GRCh38) {{
            exons {{
                feature_type
                start
                stop
            }}
        }}
    }} ''')

def query_transcripts_ensemble(ensembl_id):
    return (f''' query VariantsInGene {{
        gene(gene_id: "{ensembl_id}", reference_genome: GRCh38) {{
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
def get_unique_name(path_to_data, name="data"):
    ts = datetime.now().strftime("%Y%m%dT%H%M%SZ")
    return os.path.join(path_to_data, f"{name}_{ts}.json")

def fetch_from_ensembl_id_as_json(query, gene_path, name):
    """
    very long query string; 
    look up https://gnomad.broadinstitute.org/api to test out querys Strings
    there is also some aka little to none documentation
    # <- mark comments
    # too long, every now and then there are mistakes, so I split it in 5 querys
    """
    try:
        response = r.post("https://gnomad.broadinstitute.org/api",
                        json={"query": query},
                        timeout=180
                        )
        if response.status_code != 200:
            return None

        data = response.json().get("data", {}).get("gene", {})
        if data:
            file_name = get_unique_name(gene_path, name)
            with open(file_name, 'w') as file:
                json.dump(data, file)
            return file_name
        else:
            return False

    except Exception as e:
        send_message(f" - fetch {gene_path}/{name} failed\n{str(e)}\n")

    return None


def try_fetching_(query, gene_path, name, first=False):
    counter_for_fail = 0
    while counter_for_fail < 4:
        file = fetch_from_ensembl_id_as_json(query, gene_path, name)
        if file is None:
            counter_for_fail += 1 #download failed
            sleep(6)
        elif file:
            counter_for_fail = 99
        else:
            counter_for_fail += 2 #data was empty
            sleep(6)
    if not first:
        sleep(6)


def load_variants(files):
    variants = []
    for variant_json in files:
        if not variant_json.endswith("json") or "gnomAD_variants" not in variant_json:
            continue
        try:
            with open(variant_json, 'r') as f:
                variants += json.load(f).get("variants", None)
        except Exception as e:
            send_message(f" - couldn't load variant_file {variant_json}\n{str(e)}\n")
            continue
    return variants


def find_diff_cutoff(variant, ancestry, cutoff=0.05):
    for item in ["joint", "genome", "exome"]:
        populations_an_ac = variant.get(item, None)
        if populations_an_ac is None:
            continue

        populations = populations_an_ac.get("populations", [])
        if not populations:
            continue

        pop_af = []
        for population in populations:
            pop_id = population.get("id", None)
            if pop_id not in ancestry:
                continue

            ac = population.get("ac", None)
            if ac is None:
                pop_af.append(0)
                continue

            an = population.get("an", None)
            if an is None:
                pop_af.append(0)
                continue

            try:
                af = float(ac) / float(an)
                pop_af.append(af)
            except Exception as _:
                pop_af.append(0)
                continue

        if max(pop_af) - min(pop_af) >= cutoff:
            return True

    return False


def filter_variants(variants, population):
    long_list = extend_and_validate_ancestry_names(population)
    found = set()
    not_sure = set()
    for v in variants:
        if not isinstance(v, dict):
            continue
        variant_id = v.get("variant_id", "NO_ID")
        if find_diff_cutoff(v, long_list):
            found.add(variant_id)
        else:
            not_sure.add(variant_id)

    return found, not_sure


def send_to_VEP(gnomAD_send, files, populations, ensemble_id):
    gnomad_variants = load_variants(files)
    found, not_sure = filter_variants(gnomad_variants, populations)
    gnomAD_send.put((files, populations, found, not_sure, ensemble_id))
    send_message(len(found)+len(not_sure), 2, "vep")

def load_or_fetch(query, gene_path, file_name, download):
    if search_for_file(gene_path, file_name, ".json") is None or download:
        try_fetching_(query, gene_path, file_name)

def load_rest(already_checked, gene_path, ensembl_id, download):
    load_or_fetch(query_gen_ensemble(ensembl_id), gene_path, "gnomAD_gene", download)
    load_or_fetch(query_clinvar_ensemble(ensembl_id), gene_path, "gnomAD_clinvar", download)
    load_or_fetch(query_exons_ensemble(ensembl_id), gene_path, "gnomAD_exons", download)
    load_or_fetch(query_transcripts_ensemble(ensembl_id), gene_path, "gnomAD_transcripts", download)
    already_checked.add(ensembl_id)
    return already_checked

def download_data(gnomAD_receive, gnomAD_send, gnomAD_config):
    send_message("starting",0,"gnomad")
    data_path, populations, download = gnomAD_config
    already_checked = set()
    counter = 0
    while (True):
        try:
            if not gnomAD_receive.empty():
                ensembl_id = gnomAD_receive.get()
                if check_string(ensembl_id):
                    continue
                if ensembl_id in already_checked:
                    continue
                if ensembl_id == "finished" or counter > 3600:
                    break

                path_gen = os.path.join(data_path, ensembl_id)
                os.makedirs(path_gen, exist_ok=True)
                if download:
                    try_fetching_(query_variant_ensemble(ensembl_id), path_gen, "gnomAD_variants", True)

                files = search_for_files(path_gen, "gnomAD_variants", "json")

                if not files:
                    continue
                else:
                    send_to_VEP(gnomAD_send, files, populations, ensembl_id)

                already_checked = load_rest(already_checked, path_gen, ensembl_id, download)
                send_message(1,1,"gnomad")

            else:
                sleep(1)
                counter += 1

        except Exception as _:
            counter += 1
            sleep(1)

    gnomAD_send.put("finished")
    print("gnomAD thread done")
