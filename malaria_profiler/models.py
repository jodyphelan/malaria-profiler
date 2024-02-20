from pydantic import BaseModel
from pathogenprofiler.models import SpeciesPrediction, DrVariant, Variant, BamQC, FastaQC, VcfQC
from typing import List, Union
from pathogenprofiler import object_list2text

class Result(BaseModel):
    id: str

class SpeciesResult(Result):
    species: SpeciesPrediction
    result_type: str = 'Species'

class GeoClassificationProbability(BaseModel):
    region: str
    probability: float 

class GeoClassificationResult(BaseModel):
    classifier_name: str
    classifier_version: str
    probabilities: List[GeoClassificationProbability]

class ProfileResult(SpeciesResult):
    notes: List[str] = []
    resistance_db: dict = {}
    geo_classification: GeoClassificationResult
    dr_variants: List[DrVariant] = []
    other_variants: List[Variant] = []
    fail_variants: List[Variant] = []
    qc: Union[BamQC, FastaQC, VcfQC]
    result_type: str = 'Profile'

    def get_qc(self):
        if isinstance(self.qc, (BamQC, FastaQC)):
            text = object_list2text(l = self.qc.target_qc)
        else:
            text = "Not available for VCF input"
        return text