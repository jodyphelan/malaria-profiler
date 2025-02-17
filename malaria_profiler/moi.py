from pathogenprofiler import run_cmd
import argparse
import json

def get_moi(args:argparse.Namespace) -> dict:
    region_vcf = f"{args.files_prefix}.moi_region.vcf"
    tmp_result_file = f"{args.files_prefix}.moi_result.json"
    cmd = f"bcftools mpileup -f {args.conf['ref']} -R {args.conf['moi_regions']} {args.bam} | bcftools call -mv > {region_vcf}"
    run_cmd(cmd)
    cmd = f'pymoi --bam {args.bam} --vcf {region_vcf} --outfile {tmp_result_file}'
    run_cmd(cmd)

    d = json.load(open(tmp_result_file))
    return d
