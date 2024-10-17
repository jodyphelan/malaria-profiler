#! /usr/bin/env python

# Load useful libraries
import json
from collections import defaultdict
import argparse
import os
from tqdm import tqdm
import sys
import pathogenprofiler as pp
import pandas as pd
import subprocess as sp
import logging

logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser()
parser.add_argument('--out',type=str,help='File with samples',required = True)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default=".",type=str,help='Directory containing results')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')

args = parser.parse_args()


if args.samples:
    samples = [x.rstrip() for x in open(args.samples).readlines()]
else:
    samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

if len(samples) == 0:
    sys.stderr.write(f"No samples found in directory '{args.dir}'\n")
    sys.exit(1)

fields = ["sample_id","chrom","pos","gene_name","gene_id","change","nucleotide_change","protein_change","depth","freq","forward_reads","reverse_reads","filter","type","genotype"]
variant_rows = []
variant_info = {}
coverage_info = []

for s in tqdm(samples, desc="Loading reports"):
    # Data has the same structure as the .result.json files
    if not os.path.isfile("%s/%s%s" % (args.dir,s,args.suffix)):
        sys.stderr.write(f"Can't find sample {s}\n")
        continue
    data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
    for region in data['qc']['target_qc']:
        region['sample_id'] = s
        coverage_info.append(region)
    for var in data["dr_variants"] + data['other_variants'] + data['fail_variants']:
        var['sample_id'] = s
        if 'drugs' in var:
            var['drugs'] = ",".join([a["drug"] for a in var["drugs"]])
        else:
            var['drugs'] = ''
        var['genotype'] = None if var['filter'] == "soft_fail" else 1
        variant_rows.append({f:var[f] for f in fields})
        variant_info[var['gene_name'],var['change']] = var

df = pd.DataFrame(variant_rows)

variants = set([(v.chrom,v.pos) for v in list(df[['chrom','pos']].value_counts().reset_index(name='count').itertuples(index=False))])

depth = defaultdict(lambda: defaultdict(int))
for s in tqdm(samples, desc="Calculating depth"):
    bam = pp.filecheck(f"{args.dir}/{s}.bam")
    for l in sp.Popen(f'samtools depth {bam}',shell=True,stdout=sp.PIPE).stdout:
        chrom,pos,dp = l.decode().strip().split()
        if (chrom,int(pos)) in variants:
            depth[(chrom,int(pos))][s] = dp

mutation2genomepos = {}
samples2mutations = defaultdict(set)
for v in df.itertuples(index=False):
    mutation2genomepos[(v.gene_name,v.change)] = (v.chrom,v.pos)
    samples2mutations[v.sample_id].add((v.gene_name,v.change))





for s in samples:
    for m in mutation2genomepos:
        gpos = mutation2genomepos[m]
        if m in samples2mutations[s]:
            pass
        else:
            v = dict(variant_info[m])
            variant_rows.append({
                "sample_id": s,
                "chrom":gpos[0],
                "pos": gpos[1],
                "gene_name":v['gene_name'],
                "gene_id": v['gene_id'],
                "change": v['change'],
                "nucleotide_change": v['nucleotide_change'],
                "protein_change": v['protein_change'],
                "depth": depth[gpos][s],
                "freq": 1,
                "forward_reads": None,
                "reverse_reads": None,
                "filter":'reference',
                "type": 'reference',
                "genotype": 0,
            })

df = pd.DataFrame(variant_rows)

df.to_csv(args.out+'.variants.csv',index=False,na_rep="NA")

covdf = pd.DataFrame(coverage_info)[['sample_id','target','percent_depth_pass','median_depth']]
covdf.to_csv(args.out+'.coverage.csv',index=False,na_rep="NA")