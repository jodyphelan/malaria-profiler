#! /usr/bin/env python
from pathogenprofiler import run_cmd
from glob import glob
import os
import argparse
from joblib import Parallel, delayed
from tqdm import tqdm
import re

argparser = argparse.ArgumentParser(description='Run multiple samples in parallel')
argparser.add_argument('--input-dir', help='Directory containing fastq files',required = True)
argparser.add_argument('--output-dir', help='Directory containing fastq files',required = True)
argparser.add_argument('--experiment-id', help='Directory containing fastq files')
argparser.add_argument('--threads', default=8, help='Directory containing fastq files')

args = argparser.parse_args()

def collate_fastq_files(barcode,directory):
    files = [f for f in glob(directory + '/*.fastq.gz') if f!=f'{directory}/all.fastq.gz']
    if len(files) == 0:
        # write a dummy file
        with open(f'{directory}/all.fastq.gz','w') as fh:
            fh.write('')
    else:
        run_cmd(f'cat {" ".join(files)} > {directory}/all.fastq.gz')

    return (barcode,f'{directory}/all.fastq.gz')

def find_bardode_directories(directory):
    current_dir = os.getcwd()
    directories = [(d,os.path.join(current_dir,directory,d)) for d in os.listdir(directory)]
    directories = [d for d in directories if d[0].startswith('barcode') and os.path.isdir(d[1])]
    return directories

def get_file_structure(directory):
    files = os.listdir(directory)
    fastqs_in_directory = any([(f.endswith('.fastq.gz') or f.endswith('.fastq')) and 'barcode' in f for f in files])
    dirs = [os.path.join(directory,f) for f in files]
    dirs_in_directory = any([os.path.isdir(d) and 'barcode' in d for d in dirs])
    if fastqs_in_directory:
        return 'fastq'
    elif dirs_in_directory:
        return 'directory'
    else:
        raise ValueError(f'Could not determine file structure for {directory}')

def get_fastq_files(directory):
    files =  glob(f'{directory}/*.fastq*')
    results = []
    for f in files:
        r = re.search('barcode(\d+)',f)
        if r:
            results.append((r.group(1),f))
    return results

if not args.experiment_id:
    args.experiment_id = args.input_dir.replace('/','_')


def process_sample(barcode, fastq_file):
    run_id = f'{args.experiment_id}_{barcode}'
    run_cmd(f'malaria-profiler profile -1 {fastq_file} --resistance_db vivax_amplicon --dir {args.output_dir} -p {run_id} --platform nanopore --caller bcftools')

file_structure = get_file_structure(args.input_dir)

if file_structure == 'directory':
    print("Found directories with fastq files")
    jobs = find_bardode_directories(args.input_dir)
    parallel = Parallel(n_jobs=args.threads, return_as='generator')
    barcode_fastqs = [r for r in tqdm(parallel(delayed(collate_fastq_files)(barcode,directory) for barcode,directory in jobs),total=len(jobs),desc="Collating fastq files")]
elif file_structure == 'fastq':
    print("Using fastq files with 'barcode' in the name")
    barcode_fastqs = get_fastq_files(args.input_dir)

print("Found %d samples" % len(barcode_fastqs))

parallel = Parallel(n_jobs=args.threads, return_as='generator')
barcode_fastqs = [r for r in tqdm(parallel(delayed(process_sample)(barcode,fastq) for barcode,fastq in barcode_fastqs),total=len(barcode_fastqs),desc="Running malaria-profiler")]
