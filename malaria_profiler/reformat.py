from pathogenprofiler.models import SpeciesPrediction, Variant, BamQC, FastaQC, DrVariant, Variant, FastqQC
from typing import List, Union, Tuple
from .models import ProfileResult, GeoClassificationResult, Pipeline, SpeciesResult
from pathogenprofiler.utils import shared_dict
from pathogenprofiler import get_db
import argparse

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

def filter_missing_positions(missing_positions: List[str]) -> List[str]:
    return [ann for ann in missing_positions if len(ann.annotation)>0]

def load_amplicon_target_names(bed_file:str) -> List[str]:
    names = []
    for l in  open(bed_file,'r'):
        row = l.strip().split('\t')
        names.append(row[6])
    return names 

def create_resistance_result(
    args: argparse.Namespace,
    id: str,
    species: SpeciesPrediction,
    geo_classification: GeoClassificationResult,
    genetic_elements: List[Variant],
    qc: Union[BamQC, FastaQC],
    notes: List[str]
) -> ProfileResult:
    for var in genetic_elements:
        var.convert_to_dr_element()
    
    if hasattr(qc, 'missing_positions'):
        qc.missing_positions = filter_missing_positions(qc.missing_positions)

    pipeline = Pipeline(
        software_version=args.version,
        db_version=args.conf['version'],
        software=[{'process':k,'software':v} for k,v in shared_dict.items()]
    )

    if 'amplicon' in args.conf and args.conf['amplicon']==True:
        amplicon_target_names = load_amplicon_target_names(args.conf['bed'])
        print(amplicon_target_names)
        print(len(qc.target_qc))
        for i,name in enumerate(amplicon_target_names):
            qc.target_qc[i].target = name

    dr_variants, other_variants, fail_variants = split_variants(genetic_elements)
    data = {
        'id':id,
        'notes':notes,
        'geo_classification':geo_classification,
        'dr_variants':dr_variants,
        'other_variants':other_variants,
        'fail_variants':fail_variants,
        'species':species,
        'pipeline':pipeline,
    }
    return ProfileResult(**data, qc=qc)


def create_species_result(
    args: argparse.Namespace,
    id: str,
    species: SpeciesPrediction,
    qc: Union[FastqQC, FastaQC]
    
) -> SpeciesResult:
    args.conf = get_db(args.software_name,args.species_db)
    pipeline = Pipeline(
        software_version=args.version,
        db_version=args.conf['version'],
        software=[{'process':k,'software':v} for k,v in shared_dict.items()]
    )
    data = {
        'id':id,
        'species':species,
        'qc':qc,
        'pipeline':pipeline,
    }
    return SpeciesResult(**data)
