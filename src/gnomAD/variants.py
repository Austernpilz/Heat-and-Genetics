from collections import deque
import os
from time import sleep
import json
from datetime import datetime

import pandas as pd

from src.gnomAD.download_variants import download_data
from src.helpers.folder_magic import search_for_file

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


def found(variant, ancestry, cutoff):
    for item in ["joint", "genome", "exome"]:
        populations_an_ac = variant.get(item, None)
        if populations_an_ac is None:
            continue

        ac = populations_an_ac.get("ac", None)
        if ac is None:
            continue

        populations = populations_an_ac.get("populations", [])
        if not populations:
            continue

        for population in populations:
            pop_id = population.get("id", None)
            if pop_id not in ancestry:
                continue

            ac = population.get("ac", None)
            if ac is None:
                continue

            an = population.get("an", None)
            if an is None:
                continue

            try:
                af = float(ac) / float(an)
                # the variant exists in our population with the specific cutoff
                if af >= cutoff:
                    return True
            except Exception as _:
                continue

    return False


def get_config(config):
    storage = config.get("absolute_file_paths")
    data_path = os.path.join(storage.get("data"), "gnomAD")
    os.makedirs(data_path, exist_ok=True)
    download = config.get("flags").get("download_data")
    populations = config.get("populations")
    return (data_path, populations, download)

# def find_variant(files, ancestry):
#     variant_list = []
#     for file in files:
#         if "gnomAD_variants" in file:
#             try:
#                 with open(file, "r") as f:
#                     variant_list += json.load(f)
#             except Exception as e:
#                 print(str(e))
#     for variant in variant_list:
#         exome = variant.get("exome", {}).get("populations")
#         genome = variant.get("genome", {}).get("populations")
#         joint = variant.get("joint", {}).get("populations")
#         for pop in exome:
#             pop_id = pop.get("id", None)
#             if pop_id not in ancestry:
#                 continue


def simplify_df(gnomAD_receive, gnomAD_config):
    files = []
    data_path, populations, download = gnomAD_config
    ancestry = list(populations.keys()) + list(populations.values())

    while (True):
        try:
            ensembl_id, variants = gnomAD_receive.recv()
            if ensembl_id == "finished" and variants == "finished":
                gnomAD_receive.close()
                break
            if not ensembl_id:
                continue
            path_to_gene = os.path.join(data_path, ensembl_id)
            clean_and_filter(path_to_gene, ancestry, 0.05)

        except Exception as _:
            sleep(180)


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
        return 0, 0, 0

    ac = some_dict.get("ac", 0)
    an = some_dict.get("an", 0)
    af = 0

    try:
        af = float(ac) / float(an)
    except Exception as _:
        af = 0

    return ac, an, af

def simplified_scores(in_silico_predictors):
    scores = [
        #scores of deleterious (higher number, more devestating effect)
        "revel_max",  #0..1 high is bad
        "cadd",       #0..99 ab 20 eher kritisch

        #score for splice altering
        "pangolin_largest_ds",    #0..1 high is splice altering
        "spliceai_ds_max",        #0..1 same

        #score for amino substitution
        "sift_max",               #0..1 under 0.05 protein function is altered not there
        "polyphen_max"            #0..1 high is bad
    ]

    isp_single = {
            s: None
            for s in scores
        }

    for score in in_silico_predictors:
        id = score.get("id", None)
        if id not in scores:
            continue

        isp_single[id] = score.get("value", 0)

    return isp_single


def get_ancestry_p_and_reduce(populations_in_gen, ancestry_list):
    new_ac_an_af = {
        "ac" : 0,
        "an" : 0,
        "af" : 0
    }

    for anc in ancestry_list:
        new_ac_an_af[f"ac_{anc}"] = 0
        new_ac_an_af[f"an_{anc}"] = 0
        new_ac_an_af[f"af_{anc}"] = 0

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
        return [], "no_variant"

    variant_list = []
    for variant_json in variant_path:
        try:
            with open(variant_json, 'r') as f:
                v = json.load(f).get("variants", [])
                if v:
                    variant_list += v
        except Exception as e:
            print(str(e))
            continue

    variants_to_keep = []
    check_for_doubles = []

    while variant_list:
        variant = variant_list.pop()
        if not isinstance(variant, dict):
            continue

        variant_id = variant.get("variant_id", "NO_ID")
        if (variant_id in check_for_doubles or 
            variant_id == "NO_ID"):
            continue

        try:
            if not found(variant, ancestry_list, cut_off):
                continue

            #making the variant_new
            new = variant.copy()
            for item in ["joint", "genome", "exome"]:
                population_data = variant.get(item, None)
                new[item] = get_ancestry_p_and_reduce(population_data, ancestry_list)

            #cleaning the scores, for what i need
            new["in_silico_predictors"] = simplified_scores(variant.get("in_silico_predictors",[]))

            variants_to_keep.append( pd.json_normalize(new) )
            check_for_doubles.append(variant_id)

        except Exception as e:
            print(str(e))
            continue

    if not check_for_doubles:
        return [], "no_variant"


    try:
        df = pd.concat(variants_to_keep)
        if df.empty:
            return check_for_doubles, "no_variant"

        file_name = os.path.join(gen_folder, "simple_variants.tsv")
        df.to_csv( file_name, sep='\t', index=False)
        #print(gen_folder, " successfull")
        return check_for_doubles, file_name
    except Exception as e:
        print("couldn't save simple variants")
        print(str(e))

    return check_for_doubles, "no_variant"



def get_gnomAD_datapaths_single(path_to_gene):
    gnomAD_data = {
        "variants" : search_for_file(path_to_gene, "gnomAD_variants", "json"),
        "clinvar" : search_for_file(path_to_gene, "clinvar", "json"),
    }
    return gnomAD_data


def clean_and_filter(path_to_gene, ancestry_list, cutoff):
    gnomAD_data_dict = get_gnomAD_datapaths_single(path_to_gene)
    ancestry_list = extend_and_validate_ancestry_names(ancestry_list)

    variants_to_keep, variant_tsv = filter_variants_by_ancestry(gnomAD_data_dict["variants"], ancestry_list, cutoff, path_to_gene)
    clinvar_found, clinvar_tsv = match_variants_to_clinvar(variants_to_keep, gnomAD_data_dict["clinvar"], path_to_gene)
    #variants_sorted_out = list(set(variants_sorted_out + clinvar_variants_sorted_out))

    file_txt = os.path.join(path_to_gene, "overview.txt")
    #file_tsv = os.path.join(path_to_gene, "gene_information.tsv")
    with open (file_txt, "w") as f:
        #f.write(f"> gene_information: {str(file_tsv)}\n")
        f.write(f"> variants_found: {str(variant_tsv)}\n")
        for var in variants_to_keep:
            f.write(f"{str(var)}, ")
        f.write(f"\n> clinvar_variants_found: {str(clinvar_tsv)}\n")
        for var in clinvar_found:
            f.write(f"{str(var)}, ")
        f.write("\n> variants_sorted_out\n")
        # for var in variants_sorted_out:
        #     f.write(f"{str(var)}, ")
    print(f"{datetime.now().strftime('%H%M')} gnomAD filtered {os.path.basename(path_to_gene)}")


def match_variants_to_clinvar(gnomAD_variants, path_to_clinvar_data, gen_folder):
    if not path_to_clinvar_data:
        return [], [], "no clinvar variant"

    clinvars = []
    for variant_json in path_to_clinvar_data:
        try:
            with open(variant_json, 'r') as f:
                cv = json.load(f).get("clinvar_variants", [])
                if cv:
                    clinvars += cv
        except Exception as e:
            print("couln't read clinvar_json: ", variant_json)
            print(str(e))
            continue

    if not gnomAD_variants:
        return [], pd.json_normalize(clinvars)["variant_id"].unique().tolist(), None

    clinvar_variants_to_keep = []
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
            continue

        try:
            df = pd.json_normalize(variant).drop_duplicates()
            clinvar_variants_to_keep.append(df)
            doubles.append(var_id)
        except Exception as e:
            print(str(e))
            continue

    if not doubles:
        return [], "no clinvar variant"

    file_name = os.path.join(gen_folder, f"simple_clinvar.tsv")

    try:
        df = pd.concat(clinvar_variants_to_keep)
        if df.empty:
            return doubles, "no clinvar variant"
        df.to_csv( file_name, sep='\t', index=False)
        return doubles, file_name

    except Exception as e:
        print("couldn't save clinvar_variants")
        print(str(e))

    try:
        with open(file_name, 'w') as file:
            json.dump(clinvar_variants_to_keep, file)
        return doubles, file_name
    except Exception as e:
        print("match_variants_to_clinvar4")
        print(str(e))

    return doubles, "no clinvar variant"


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

