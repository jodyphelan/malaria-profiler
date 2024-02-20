from pathogenprofiler.models import SpeciesPrediction, Variant, BamQC, FastaQC, DrVariant, Variant
from typing import List, Union, Tuple
from .models import ProfileResult, GeoClassificationResult

def split_variants(
    variants: List[Variant]
) -> Tuple[List[DrVariant], List[Variant]]:

    dr_variants = []
    other_variants = []
    fail_variants = []
    for var in variants:
        if var.filter.upper() == "PASS":
            if isinstance(var, DrVariant):
                dr_variants.append(var)
            else:
                other_variants.append(var)
        else:
            fail_variants.append(var)
    return dr_variants,other_variants,fail_variants

def create_resistance_result(
    id: str,
    species: SpeciesPrediction,
    geo_classification: GeoClassificationResult,
    genetic_elements: List[Variant],
    qc: Union[BamQC, FastaQC],
    notes: List[str]
) -> ProfileResult:
    for var in genetic_elements:
        print(var)
        var.convert_to_dr_element()
    dr_variants, other_variants, fail_variants = split_variants(genetic_elements)
    data = {
        'id':id,
        'notes':notes,
        'geo_classification':geo_classification,
        'dr_variants':dr_variants,
        'other_variants':other_variants,
        'fail_variants':fail_variants,
        'species':species,
    }
    return ProfileResult(**data, qc=qc)



# def reformat(results,conf):
#     results["variants"] = [x for x in results["variants"] if len(x["consequences"])>0]
#     results["variants"] = select_csq(results["variants"])
#     results["variants"] = dict_list_add_genes(results["variants"],conf)
#     results = reformat_annotations(results,conf)

#     if "region_qc" in results["qc"]:
#         if conf['amplicon']==True:
#             amplicon2gene = load_bed(conf['bed'],columns=[4,5,7],key1=7)
#             for d in results["qc"]["region_qc"]:
#                 d["gene_id"] = amplicon2gene[d["region"]][0]
#             results["qc"]["region_qc"] = dict_list_add_genes(results["qc"]["region_qc"],conf)
#         else:
#             for d in results["qc"]["region_qc"]:
#                 d["gene_id"] = d['region']
#     if "missing_positions" in results["qc"]:
#         results["qc"]["missing_positions"] = reformat_missing_genome_pos(results["qc"]["missing_positions"],conf)
    
#     return results
