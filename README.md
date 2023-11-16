# malaria-profiler

## Install

```
# Run the following depending on your platorm

# linux
wget  https://raw.githubusercontent.com/jodyphelan/malaria-profiler/main/conda/linux.env.txt
conda create --name malaria-profiler --file linux.env.txt

# macos
wget  https://raw.githubusercontent.com/jodyphelan/malaria-profiler/main/conda/macos.env.txt
conda create --name malaria-profiler --file macos.env.txt

# Then the following commands on all platforms
conda activate malaria-profiler
pip install git+https://github.com/jodyphelan/malaria-profiler.git
pip install git+https://github.com/jodyphelan/pathogen-profiler.git
malaria-profiler update_db
```

## Updating

```
conda activate malaria-profiler
pip install --force-reinstall git+https://github.com/jodyphelan/malaria-profiler.git
pip install --force-reinstall git+https://github.com/jodyphelan/pathogen-profiler.git
```

## Usage

### Input types

Malaria-Profiler species prediciton is currently available to run on a fastq, bam, cram, fasta or vcf data. The output is a txt file with the species prediction and if there is a resistance database then it will also output a  list of variants and if they have been associated with drug resistance.

#### Fastq data

Raw sequencing data in fastq format can been used as input using the following command. The second read is optional.

```bash
malaria-profiler profile -1 </path/to/reads_1.fq.gz> -2 </path/to/reads_2.fq.gz> -p <sample_name> -t [threads] --txt 
```

#### Bam/Cram 

Aligned data in the form of bam or cram files can be used. Please note that the alignment files **must** have been generated with the same reference genome (even the chromosome names) as those used by `malaria-profiler` database.

```bash
malaria-profiler -a </path/to/bam/cram> -p <sample_name> -t [threads] --txt
```

#### Fasta 

Assembled genomes or gene sequencves in fasta format can been used as input using the following command.

```bash
malaria-profiler -f </path/to/fasta> -p <sample_name> -t [threads] --txt
```

#### VCF 

Varaints stored in VCF format can been used as input using the following command. Again the chromosome names must match those in the species-specific database.

```bash
malaria-profiler -v </path/to/vcf> -p <sample_name> -t [threads] --txt
```

#### General options

If you have used a reference genome with different sequence names that you have used to generate a bam/cram/vcf then it is possible to align the `malaria-profiler` databases to use the same sequence names. Please go to the custom databases section to find out more.

Other useful options arguments include 
* `--threads` - sets the number of parallel threads
* `--platform` - sets the platform that was used to generate the data (default=illumina) 
* `--txt` - outputs a text based report

A full list of arguments can be found by running `malaria-profiler profile -h`

### Collating results across runs

The results from numerous runs can be collated into one table using the following command.

```bash
malaria-profiler collate 
```

## How it works?

### Species prediction
Species prediction is performed by looking for pre-detemined kmers in read files which belong to a specific species. If no species is found using this method, mash is run using a database of all Plasmodium mitochondrial sequences from GTDB to find the top 10 closest genomes.

### Resistance prediction
Resistance prediction is performed by aligning the read data to a species-specific reference genome and looking for resistance associated genes and variants. The reference and resistance database is stored in the [malaria-db github repo](https://github.com/jodyphelan/malaria-db). At the moment resistance prediction is available for:

* _Plasmodium falciparum_

If you would like to suggest another organism please leave a comment in [this thread](https://github.com/jodyphelan/malaria-profiler/).
