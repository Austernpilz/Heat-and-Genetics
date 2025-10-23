import requests as r
import asyncio #because the loading times are so long, i use parallel and async processes
import os
from datetime import datetime
import json
from collections import deque
import pandas as pd
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
def query_ensemble_id(ensemble_id):
    return (f'''
query VariantsInGene {{
  gene(gene_id: "{ensemble_id}", reference_genome: GRCh38){{
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
    
  }}
}}
''')
def query_chr_pos(chromosom, start, stop):
    # creating a very long, f string to use as a function argument,
    # maybe there could be a better way, 
    # but this is straigt out of the documentation how a query in python should look like
    # and this was the easiest to make a function out of it
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
}}
''')

#unique names in case we don't get an gene id or symbol
def get_unique_name(path_to_data, name):
    ts = datetime.now().strftime("%Y%m%dT%H%M%SZ")
    return os.path.join(path_to_data, f"{name}_{ts}.json")

def fetch_data_as_json(chromosome, pos_start, pos_end):
    #chromosome = f'\"{chromosome}\"'
    try:
        response = r.post("https://gnomad.broadinstitute.org/api",
                        json={"query": query_chr_pos(chromosome, pos_start, pos_end)}
                        )
        return response.json()
    except Exception as e:
        print(str(e))
        return {}

def fetch_data_as_json(chromosome, pos_start, pos_end):
    #chromosome = f'\"{chromosome}\"'
    try:
        response = r.post("https://gnomad.broadinstitute.org/api",
                        json={"query": query_chr_pos(chromosome, pos_start, pos_end)}
                        )
        return response.json()
    except Exception as e:
        print(str(e))
        return {}

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

# async def generate_task(name, chromosom, start, end, tasks=deque()):
#     tasks.append((
#         (name, chromosom, start, end),
#         asyncio.create_task(fetch_data_as_json(chromosom, start, end))
#         ))
#     await asyncio.sleep(6)

#async 
def download_data(df, path_to_gnomAD):
    result_dict = {}
    tasks = deque()

    #generate fetch requests
    for name, chromosom, start, end in zip(df['genes'].tolist(), df['chr'].tolist(), df['start'].tolist(), df['end'].tolist()):
        tasks.append((name, chromosom, start, end))

    data_path = os.path.join(path_to_gnomAD, "data")
    os.makedirs(data_path, exist_ok=True)

    while tasks:
        start = datetime.now()
        gen_symbol, c, b, e = tasks.pop()
        path_gen = os.path.join(data_path, gen_symbol)
        try:
            res = fetch_data_as_json(c,b,e)
            variants = res.get("data", {}).get("region", {}).get("variants", [])
            genes = res.get("data", {}).get("region", {}).get("genes", [])

            if gen_symbol not in result_dict:
                result_dict[gen_symbol] = {}

            if save_json(variants, path_gen, f"{gen_symbol}_variants"):
                result_dict[gen_symbol]["variants"] = variants

            if save_json(genes, path_gen, f"{gen_symbol}_genes"):
                result_dict[gen_symbol]["genes"] = genes

        except Exception as e:
            print(str(e))
            if gen_symbol not in result_dict:
                #await generate_task(n,c,b,e, tasks)
                result_dict[gen_symbol] = {}
        # if datetime.now() - start:
    return result_dict

def test():
    gnomAD = os.path.join(os.getcwd(), "gnomAD")

    df = pd.DataFrame({
        "genes" : ["ENST00000644486", "ATXN3", "PPP5C"],
        "chr" : ["14", "14","19"],
        "start": [92058552, 92044496, 46347108, ],
        "end": [92106582, 92106622, 46390852],
        })

    download_data(df, gnomAD)

#never use print with gnomAD data :'(
test()


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