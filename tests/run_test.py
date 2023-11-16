from pathogenprofiler import run_cmd
import pathogenprofiler as pp
import json
import os

def extract_mutations(result):
    return [(d['gene'],d['change']) for d in result["dr_variants"]]

example = extract_mutations(json.load(open("ERR039929.results.json")))

if not os.path.isdir("scratch"):
    os.mkdir("scratch")
os.chdir("scratch")

def get_variants(data,type):
    return [(x["gene"],x["change"]) for x in data[type]]

def test_update_library():
    run_cmd("malaria-profiler update_db ")

# def test_fasta_profile():
#     run_cmd("malaria-profiler profile -f ~/test_data/ERR039929.contigs.fa -p ERR039929_fasta -t 3")
#     result = json.load(open("ERR039929_fasta.results.json"))
#     assert result["resistance_genes"] == example["resistance_genes"]
#     assert get_variants(example,"dr_variants")== get_variants(result,"dr_variants")
#     assert get_variants(example,"other_variants")== get_variants(result,"other_variants")

def test_profile():
    run_cmd("malaria-profiler profile -1 ~/test_data/ERR039929_1.fastq.gz -2 ~/test_data/ERR039929_2.fastq.gz -p ERR039929 -t 3 --txt --csv")
    result = extract_mutations(json.load(open("ERR039929.results.json")))
    assert result == example


def test_clean():
    os.chdir("../")
    run_cmd("rm -r scratch")