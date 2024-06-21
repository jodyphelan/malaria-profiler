from pydantic import BaseModel
from pathogenprofiler.models import SpeciesPrediction, DrVariant, Variant, BamQC, FastaQC, VcfQC
from typing import List, Union
from pathogenprofiler import object_list2text

class Pipeline(BaseModel):
    software_version: str
    db_version: dict
    software: List[dict]

class Result(BaseModel):
    id: str
    pipeline: Pipeline

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
    qc: Union[BamQC, FastaQC, VcfQC]
    geo_classification: Union[GeoClassificationResult, None]
    dr_variants: List[DrVariant] = []
    other_variants: List[Variant] = []
    fail_variants: List[Variant] = []
    result_type: str = 'Profile'

    def get_qc(self):
        if isinstance(self.qc, (BamQC, FastaQC)):
            text = object_list2text(l = self.qc.target_qc)
        else:
            text = "Not available for VCF input"
        return text