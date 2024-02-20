import sys
import pathogenprofiler as pp
import json
import argparse
from pathogenprofiler.models import SpeciesPrediction

def process_args(args: argparse.Namespace) -> None:
    args.no_delly = False if args.run_delly else True

def get_conf_dict_with_path(library_path):
    files = {"ref":".fasta","barcode":".barcode.bed","bed":".bed","json_db":".dr.json","version":".version.json","variables":".variables.json"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_path+files[key]))
        if key=="variables":
            tmp = json.load(open(pp.filecheck(library_path+files[key])))
            for k in tmp:
                conf[k] = tmp[k]
        else:
            conf[key] = pp.filecheck(library_path+files[key])
    return conf

def get_conf_dict(library_prefix):
    library_prefix = "%s/share/malaria-profiler/%s" % (sys.base_prefix,library_prefix)
    return get_conf_dict_with_path(library_prefix)


def get_species(args: argparse.Namespace) -> SpeciesPrediction:
    if args.resistance_db:
        return pp.set_species(args)
    else:
        return pp.get_sourmash_species_prediction(args)

