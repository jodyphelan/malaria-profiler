from pathogenprofiler import bam_profiler, get_barcoding_mutations, Bam, Vcf
from copy import deepcopy
import os

def assign_region(snps_report, conf):
    condition = set([(r[3],int(r[5])) for r in [l.strip().split() for l in open(conf['geo_barcode'])]])
    region_result = []
    for item in condition:
        sum_snps = sum([it[5] for it in snps_report if it[0] == item[0]])
        region_result.append([item[0], item[1], sum_snps])
    res = [x[0] for x in region_result if x[2] >= x[1]]

    if "geoclassification_order" in conf:
        order = conf["geoclassification_order"]
        if len(res)>1:
            region = [(sorted(res,key=lambda x:order.index(x)))[0]]
        elif len(res) == 1:
            region = res
        else:
            region = ["Unassigned"]
    else:
        region = res

    return(region)




def add_geobarcode(args):
    
    if "bam" in vars(args):
        bam = Bam(args.bam,prefix=args.files_prefix,platform=args.platform)
        barcode_mutations = bam.get_bed_gt(bed_file=args.conf['geo_barcode'],ref_file=args.conf['ref'],caller=args.caller,platform=args.platform)
    elif args.vcf:
        vcf = Vcf(args.vcf_file,prefix=args.files_prefix)
        barcode_mutations = vcf.get_bed_gt(args.conf['geo_barcode'],args.conf['ref'])
    elif args.fasta:
        vcf = Vcf(f"{args.files_prefix}.vcf.gz",prefix=args.files_prefix)
        barcode_mutations = vcf.get_bed_gt(args.conf['geo_barcode'],args.conf['ref'])
    else:
        raise Exception("No input file specified")
    barcode_support,snps_report = get_barcoding_mutations(barcode_mutations,args.conf["geo_barcode"])
    region_result = assign_region(snps_report, args.conf)
    return {'geoclassification':region_result}