from .utils import *
from .reformat import *
from .output import *
from .speciation import *
from .geo_classifier import *
from .moi import *
from pathogenprofiler.variant_calling import VariantCaller
__version__="0.0.10"




class MalariaFreebayesCaller(VariantCaller):
    __software__ = "freebayes-illumina-malaria"
    def call_variants(self) -> Vcf:
        # Call variants using Lofreq
        
        self.calling_cmd = """
            samtools view -T %(ref_file)s  -h %(bam_file)s {region} %(samclip_cmd)s | \
            samtools view -b > %(temp_file_prefix)s.{region_safe}.tmp.bam && \
            samtools index %(temp_file_prefix)s.{region_safe}.tmp.bam && \
            freebayes -f %(ref_file)s -r {region} %(calling_params)s  %(temp_file_prefix)s.{region_safe}.tmp.bam | \
            bcftools norm -a -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz
            """ % vars(self)
        
        return self.run_calling(self.calling_cmd)