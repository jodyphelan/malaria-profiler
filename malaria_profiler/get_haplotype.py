import argparse
import pathogenprofiler as pp
from collections import Counter, OrderedDict
import csv

def get_haplotype(args):
    if "falciparum" in args.conf["species"]:
        positions = OrderedDict()
    # read in bed positions
        for l in open(args.conf['geo_barcode']):
            row = l.strip().split('\t')
            positions[(row[0],int(row[2]))] = 1
        if "bam_file" in vars(args):
            bam_class = pp.bam(args.bam_file,prefix=args.files_prefix,platform=args.platform)
            mutations = bam_class.get_bed_gt(bed_file=args.conf['geo_barcode'],ref_file=args.conf['ref'],caller=args.caller,platform=args.platform)
        elif args.fasta:
            vcf_class = pp.vcf(f"{args.files_prefix}.vcf.gz",prefix=args.files_prefix)
            mutations = vcf_class.get_bed_gt(args.conf['geo_barcode'],args.conf['ref'])
        elif "vcf_file" in vars(args):
            vcf_class = pp.vcf(args.vcf_file,prefix=args.files_prefix)
            mutations = vcf_class.get_bed_gt(args.conf['geo_barcode'],args.conf['ref'])
        
            
    #bam_obj = pp.bam(args.bam,'_','Nanopore')
    #mutations = bam_obj.get_bed_gt(args.bed,args.ref, caller="bcftools",platform="Nanopore")

        barcode = []
        for chrom,pos in positions:
        
            total_cov = sum(mutations[chrom][pos].values())
            count = Counter(mutations[chrom][pos])
            #print(chrom,pos,count)
            most_common = count.most_common()[0]
            if total_cov < 10:
                barcode.append('N')
            elif most_common[1]/total_cov > 0.8:
                barcode.append(most_common[0][0])
            else:
                barcode.append('N')

        seq = ''.join(barcode)


        from itertools import product

        d = { "A":["A"],"C":["C"],"G":["G"],"T":["T"],"N": ["A", "G", "T", "C"] }

        possible_seqs = [ "".join(i) for i in product(*[ d[j] for j in seq ]) ]
        #print(possible_seqs)
        #print(seq)
        
        print(args.conf)
        output_rows = []
        for row in csv.DictReader(open("haplotypes.csv")):
            
            if row['haplotype'] in possible_seqs:
                output_rows.append(row)
        

        
        for haplotype_data in output_rows:
            haplotype = haplotype_data['haplotype']
            regions = []

            for key, value in haplotype_data.items():
                if key not in ['haplotype', 'nsamples']:
                    if float(value) > 0.0:
                        regions.append([key, float(value)])
            
        
        return {'haplotype_geoclassification':regions}
    else:
        return {'haplotype_geoclassification':"non_falciparum"}

        