import pandas as pd
# from comet_search import run_comet
# from search_output_reconciliation import reconcile_search_results
# from annotate_psms_via_models import annotate_pin_file
from run_mokapot import concatenate_psm_files, run_mokapot
# from generate_final_report import generate_PSM_report


varmods = [{'mass': 15.994915,
            'residues': 'M',
            'binary': 0,
            'max_mods_per_peptide': 5,
            'term_distance': -1,
            'N/C-term': 0,
            'required': 0},
           {'mass': 42.010565,
            'residues': 'n',
            'binary': 0,
            'max_mods_per_peptide': 2,
            'term_distance': 0,
            'N/C-term': 0,
            'required': 0}]

if __name__ == '__main__':
#     rawfile = '/data/test_data/draftprot_testcase/Fetal_Ovary_bRP_Elite_25_f08.raw'
#     fastafile = '/data/gencode/r43_h38p13/gencode.v43.pc_translations.with_pseudorev.fasta'

#     rawfile = '/data/test_data/draftprot_testcase/Fetal_Ovary_bRP_Elite_25_f08.raw'
#     mgffile = '/data/test_data/draftprot_testcase/Fetal_Ovary_bRP_Elite_25_f08.mgf'

#     sagefile = '/data/test_data/draftprot_testcase/results.sage.pin'

#     # comet_pin, comet_txt = run_comet(mgffile, fastafile, varmods=varmods)
#     # print((comet_pin, comet_txt))
#     comet_pin, comet_txt = '/data/test_data/draftprot_testcase/Fetal_Ovary_bRP_Elite_25_f08.pin', '/data/test_data/draftprot_testcase/Fetal_Ovary_bRP_Elite_25_f08.txt'    
#     # rec_file = reconcile_search_results(comet_pin, comet_txt, sagefile)
#     # rec_file = '/data/test_data/draftprot_testcase/Fetal_Ovary_bRP_Elite_25_f08.reconciled.pin' 
#     # ann_pin = annotate_pin_file(rec_file, rawfile)
#     ann_pin = '/data/test_data/draftprot_testcase/Fetal_Ovary_bRP_Elite_25_f08.reconciled.annotated.pin'
#     # mokpsm = run_mokapot(ann_pin)
#     # print(mokpsm)
#     
#     mokpsm = '/data/test_data/draftprot_testcase/mokapot.psms.txt'

    ann_pin = '/home/max/biostuff/proteomic_pipeline/work/PTRC_exp12_plex_05_G_f02.reconciled.annotated.psm'
    mokpsm = run_mokapot(ann_pin)
    # report = generate_PSM_report(ann_pin, mokpsm)
    
    print(mokpsm)
    foo = pd.read_csv(str(mokpsm), sep='\t')
    print(foo[foo['mokapot q-value'] < 0.01].shape)
    print("Done")



