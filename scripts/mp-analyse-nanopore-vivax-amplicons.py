#! /usr/bin/env python
from pathogenprofiler import run_cmd
from glob import glob
import os
import argparse
import pandas as pd

argparser = argparse.ArgumentParser(description='Run multiple samples in parallel')
argparser.add_argument('--input_directory', help='Directory containing fastq files',required = True)
argparser.add_argument('--output_directory', help='Directory containing fastq files',required = True)
argparser.add_argument('--experiment-id', help='Directory containing fastq files')

args = argparser.parse_args()

def collate_fastq_files(directory):
    files = [f for f in glob(directory + '/*.fastq.gz') if f!=f'{directory}/all.fastq.gz']
    if len(files) == 0:
        # write a dummy file
        with open(f'{directory}/all.fastq.gz','w') as fh:
            fh.write('')
    else:
        run_cmd(f'cat {" ".join(files)} > {directory}/all.fastq.gz')

def find_bardode_directories(directory):
    current_dir = os.getcwd()
    directories = [(d,os.path.join(current_dir,directory,d)) for d in os.listdir(directory)]
    directories = [d for d in directories if d[0].startswith('barcode') and os.path.isdir(d[1])]
    return directories

def load_seqkit_stats(file):
    # file	format	type	num_seqs	sum_len	min_len	avg_len	max_len
    # /Users/jody/temp/malaria/vivax/fastq_pass/barcode02/all.fastq.gz	FASTQ	DNA	32000	11105214	86	347.0	176071
    with open(file) as fh:
        header = next(fh).strip().split('\t')
        data = next(fh).strip().split('\t')
        data = dict(zip(header,data))
    return data


if not args.experiment_id:
    args.experiment_id = args.input_directory.replace('/','_')

final_data = []
for barcode, directory in find_bardode_directories(args.input_directory):
    print(barcode, directory)
    run_id = f'{args.experiment_id}_{barcode}'
    collate_fastq_files(directory)
    # run_cmd(f'malaria run {directory}/all.fastq.gz {directory}/results')
    run_cmd(f'seqkit stats -T {directory}/all.fastq.gz > {directory}/stats.txt')
    stats = {
        'barcode': barcode
    }
    stats.update(load_seqkit_stats(f'{directory}/stats.txt'))
    final_data.append(stats)
    run_cmd(f'malaria-profiler profile -1 {directory}/all.fastq.gz --resistance_db vivax_amplicon --dir {args.output_directory} -p {run_id} --platform nanopore --caller bcftools')

df = pd.DataFrame(final_data)
df.to_csv(f'{args.output_directory}/{args.experiment_id}.fastq_stats.csv',index=False)


