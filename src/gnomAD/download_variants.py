
from datetime import datetime
import os
import json
from time import sleep
import requests as r
import numpy as np

from src.helpers.folder_magic import search_for_files

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
    #ts = datetime.now().strftime("%Y%m%dT%H%M%SZ")
    return os.path.join(path_to_data, f"{name}.json")

def fetch_from_ensembl_id_as_json(query, ensembl_id, data_path, name, download=False):
    """
    very long query string; 
    look up https://gnomad.broadinstitute.org/api to test out querys Strings
    there is also some aka little to none documentation
    # <- mark comments
    # too long, every now and then there are mistakes, so I split it in 5 querys
    """
    print("gnomad hängt?")
    gnomAD_data = os.path.join(data_path, ensembl_id)
    os.makedirs(gnomAD_data, exist_ok=True)
    file_name = get_unique_name(gnomAD_data, name)
    if not os.path.isfile(file_name) or download:
        try:
            #time_start = datetime.now()
            response = r.post("https://gnomad.broadinstitute.org/api",
                            json={"query": query},
                            timeout=180
                            )
            #time_elapsed = (time_start - datetime.now()).total_seconds()
            #if time_elapsed < 6:
            #    sleep(6 - time_elapsed)
            sleep(6)
            print("sleep was succesful")
            data = response.json().get("data", {}).get("gene", {})

            if data:
                with open(file_name, 'w') as file:
                    json.dump(data, file)
                return file_name
            else:
                return False

        except Exception as e:
            print(f"\n\n fetch {name} failed")
            print(str(e))

        return None


def try_fetching_(query, id, data_path, name, files=[]):
    counter_for_fail = 0
    while counter_for_fail < 4:
        file = fetch_from_ensembl_id_as_json(query, id, data_path, name)
        if file is None:
            counter_for_fail += 1 #download failed
        elif file:
            counter_for_fail = 99
            files.append(file)
        else:
            counter_for_fail += 2 #data was empty
    return files

def found(variant_path, ancestry, cutoff):
    variant_list = []
    try:
        with open(variant_path, 'r') as f:
            variant_list = json.load(f).get("variants", None)
            if variant_list is None:
                return False
    except Exception as e:
        print(str(e))
        return False

    for variant in variant_list:
        for item in ["joint", "genome", "exome"]:
            populations_an_ac = variant.get(item, None)
            if populations_an_ac is None:
                continue

            populations = populations_an_ac.get("populations", [])
            for population in populations:
                pop_id = population.get("id", None)
                if pop_id not in ancestry:
                    continue

                ac = population.get("ac", 0)
                if ac == 0:
                    continue

                an = population.get("an", 0)
                if an == 0:
                    continue

                try:
                    af = float(ac) / float(an)
                    # the variant exists in our population with the specific cutoff
                    if af >= cutoff:
                        return True
                except Exception as _:
                    continue

    return False


def download_data(gnomAD_receive, gnomAD_send, gnomAD_config):
    ensembl_id = ""
    data_path, populations, download = gnomAD_config
    #ancestry = list(populations.keys()) + list(populations.values())
    already_checked = set()
    print("starting gnomAD")
    counter = 0
    while (True):
        try:
            if gnomAD_receive.poll(timeout=60):
                ensembl_id = gnomAD_receive.recv()
            else:
                counter += 1
                ensembl_id = "NO ID"

            if ensembl_id == "finished":
                gnomAD_receive.close()
                break
            elif (
                ensembl_id == "NO ID" or 
                ensembl_id == np.nan or 
                ensembl_id == "nan" or 
                not isinstance(ensembl_id, str) or
                not "ENS" in ensembl_id or
                ensembl_id in already_checked
                ):
                if counter > 10:
                    gnomAD_receive.close()
                    break
                else:
                    continue
            if not download:
                path_gen = os.path.join(data_path, ensembl_id)
                files = search_for_files(path_gen, "", "json")
                if len(files) > 1:
                    gnomAD_send.send(files)
                    already_checked.add(ensembl_id)
                    print(f"{datetime.now().strftime('%H%M')} gnomAD already got {ensembl_id}")
                    continue
        except Exception as _:
            sleep(1)

        files = try_fetching_(query_variant_ensemble(ensembl_id), ensembl_id, data_path, "gnomAD_variants", [])
        if not files:
            #no variants, no population data, no variants to look at
            continue

        files = try_fetching_(query_gen_ensemble(ensembl_id), ensembl_id, data_path, "gnomAD_gene", files)
        files = try_fetching_(query_clinvar_ensemble(ensembl_id), ensembl_id, data_path, "gnomAD_clinvar", files)
        files = try_fetching_(query_exons_ensemble(ensembl_id), ensembl_id, data_path, "gnomAD_exons", files)
        files = try_fetching_(query_transcripts_ensemble(ensembl_id), ensembl_id, data_path, "gnomAD_transcripts", files)
        gnomAD_send.send(files)
        already_checked.add(ensembl_id)
        print(f"{datetime.now().strftime('%H%M')} gnomAD got {ensembl_id}")

    gnomAD_send.send("finished")
    gnomAD_send.close()
    print("gnomAD thread done")
