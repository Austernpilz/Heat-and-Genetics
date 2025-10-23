import requests as r
import asyncio #because the loading times are so long, i use parallel and async processes
import os
from datetime import datetime
import json
from collections import deque
import pandas as pd
from time import sleep

# broad_information
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


def query_chr_pos(chromosom, start, stop):
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
#                         json={"query": query_chr_pos(chromosome, pos_start, pos_end)}
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


def download_data(df: list, path_to_gnomAD):
    result_dict = {}
    tasks = deque()
    data_path = os.path.join(path_to_gnomAD, "data")
    os.makedirs(data_path, exist_ok=True)
    for ensemble_id in df:
        tasks.append(ensemble_id)

    while tasks:
        start = datetime.now()
        ensemble_id = tasks.pop()
        path_gen = os.path.join(data_path, ensemble_id)

        res = fetch_data_as_json(ensemble_id)
        gene_ids = res.get("data", {}).get("gene", {})
        variants = gene_ids.pop("variants", [])
        clinvar_variants = gene_ids.pop("clinvar_variants", [])

        if ensemble_id not in result_dict:
            result_dict[ensemble_id] = {}

        if save_json(variants, path_gen, f"{ensemble_id}_variants"):
            result_dict[ensemble_id]["variants"] = variants

        if save_json(gene_ids, path_gen, f"{ensemble_id}_gene_ids"):
            result_dict[ensemble_id]["genes"] = gene_ids

        if save_json(clinvar_variants, path_gen, f"{ensemble_id}_clinvar_variants"):
            result_dict[ensemble_id]["clinvar_variants"] = clinvar_variants

        if any(item not in result_dict[ensemble_id]
               for item in ["clinvar_variants", "genes", "variants"]):
            # try twice
            tasks.append(ensemble_id)
            result_dict[ensemble_id]["clinvar_variants"] = []
            result_dict[ensemble_id]["variants"] = []
            result_dict[ensemble_id]["genes"] = {}
        end = 6 - (start - datetime.now()).total_seconds()
        if end > 0:
            sleep(end)

    return result_dict

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

def test():
    gnomAD = os.path.join(os.getcwd(), "gnomAD")

    # df = pd.DataFrame({
    #     "genes" : ["ENST00000644486", "ATXN3", "PPP5C"],
    #     "chr" : ["14", "14","19"],
    #     "start": [92058552, 92044496, 46347108, ],
    #     "end": [92106582, 92106622, 46390852],
    #     })
    df = ["ENSG00000008869", "ENSG00000066427", "ENSG00000011485"]
    download_data(df, gnomAD)

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
