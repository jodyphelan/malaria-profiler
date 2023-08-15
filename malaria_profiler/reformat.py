from pathogenprofiler import load_bed, dict_list_add_genes, reformat_annotations, reformat_missing_genome_pos, select_csq

def reformat(results,conf):
    results["variants"] = [x for x in results["variants"] if len(x["consequences"])>0]
    results["variants"] = select_csq(results["variants"])
    results["variants"] = dict_list_add_genes(results["variants"],conf)
    results = reformat_annotations(results,conf)

    if "region_qc" in results["qc"]:
        if conf['amplicon']==True:
            amplicon2gene = load_bed(conf['bed'],columns=[4,5,7],key1=7)
            for d in results["qc"]["region_qc"]:
                d["gene_id"] = amplicon2gene[d["region"]][0]
        results["qc"]["region_qc"] = dict_list_add_genes(results["qc"]["region_qc"],conf)
    if "missing_positions" in results["qc"]:
        results["qc"]["missing_positions"] = reformat_missing_genome_pos(results["qc"]["missing_positions"],conf)
    
    return results