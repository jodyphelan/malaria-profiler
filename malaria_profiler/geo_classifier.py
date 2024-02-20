import pickle
from typing import Dict, List, Tuple
from pathogenprofiler import GenomePosition, Bam, Vcf
import argparse
from .models import GeoClassificationResult, GeoClassificationProbability

class Geopredictor:
    def __init__(self,clf,positions: List[Tuple[str, int]], name: str, version: str):
        self.clf = clf
        self.positions = list(positions)
        self.name = name
        self.version = version
    def predict(self,X: List[int]):
        probs = self.clf.predict_proba([X])[0]
        results = []
        for region,prob in zip(self.clf.classes_,probs):
            results.append(
                GeoClassificationProbability(
                    region=region,
                    probability=prob
                )
            )
        result = GeoClassificationResult(
            classifier_name=self.name,
            classifier_version=self.version,
            probabilities=results
        )
        return result
    def get_positions(self):
        return self.positions
    def write_bed(self,outfile: str):
        with open(outfile,'w') as O:
            for chrom,pos,alt in self.get_positions():
                O.write(f'{chrom}\t{pos-1}\t{pos}\t{alt}\n')

def get_barcoding_mutations(args: argparse.Namespace,bed_file:str):
    if args.bam:
        bam = Bam(args.bam, prefix=args.files_prefix, platform=args.platform)
        mutations = bam.get_bed_gt(
            bed_file=bed_file,
            ref_file=args.conf['ref'],
            caller=args.caller,
            platform=args.platform
        )
    elif args.vcf:
        vcf = Vcf(args.vcf)
        mutations = vcf.get_bed_gt(
            bed_file=bed_file,
            ref_file=args.conf['ref'],
        )
    else:
        raise ValueError('Please provide either a bam or vcf file')
    return mutations


def predict_geographic_source(args: argparse.Namespace):

    gp: Geopredictor = pickle.load(open(args.conf['geographic_model'], 'rb'))
    print(gp.positions)
    model_positions = gp.get_positions()
    bed_file = f'{args.files_prefix}.barcode.bed'
    gp.write_bed(bed_file)
    mutations = get_barcoding_mutations(args,bed_file)
    print(mutations)
    genotype_vector = []
    print(list(model_positions))
    for chrom,pos,alt in model_positions:
        p = GenomePosition(chrom=chrom,pos=int(pos))
        print(p)
        if p in mutations:
            if mutations[p].get(alt,0) > 0:
                genotype_vector.append(1)
            else:
                genotype_vector.append(0)
        else:
            genotype_vector.append(0)
    print(genotype_vector)
    return gp.predict(genotype_vector)

