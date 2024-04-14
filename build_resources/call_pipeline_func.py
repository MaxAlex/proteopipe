import os
import argparse
import json

# Parameters are re-used rather sloppily, and don't have defined meanings across
# functions, which is a bit gross. Eventually will get rid of this and have each
# script do its own arg parsing.
def run_pipeline_func():
    parser = argparse.ArgumentParser(description="Run a proteomics pipeline function")
    parser.add_argument("FUNCTION", type=str, help="Function to run")
    parser.add_argument('-d', '--data', type=str, nargs='+')
    parser.add_argument('-p', '--pin', type=str)
    parser.add_argument('-t', '--txt')
    parser.add_argument('-g', '--sagepin')
    parser.add_argument('-s', '--psm', nargs='+')
    parser.add_argument('-c', '--config')  # Situational config file
    parser.add_argument('-ppp', '--params')  # Comet parameter file
    parser.add_argument('-f', '--fasta')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()

    if args.FUNCTION == 'comet':
        from comet_search import run_comet
        # Will always output two files of form 'FILE.mgf -> FILE.pin + FILE.txt'
        comet_params = json.load(open(args.params))
        comet_varmods = json.load(open(args.config))
        if isinstance(args.data, list):
            mgf_file = args.data[0]
        else:
            mgf_file = args.data
        run_comet(mgf_file, args.fasta, basic_config=True, varmods=comet_varmods, **comet_params)
    elif args.FUNCTION == 'reconcile':
        from search_output_reconciliation import reconcile_search_results
        reconcile_search_results(args.pin, args.txt, args.sagepin, args.output)
    elif args.FUNCTION == 'annotate':
        from annotate_psms_via_models import annotate_pin_file
        if isinstance(args.data, list):
            raw = args.data[0]
        else:
            raw = args.data
        try:
            # Required for thermo reasons - mzAPI fails when called on symlinks
            raw = os.readlink(raw)
        except OSError:
            pass
        
        annotate_pin_file(args.pin, raw, args.output)
    elif args.FUNCTION == 'concatenate':
        from run_mokapot import concatenate_psm_files
        concatenate_psm_files(args.psm, outputfile=args.output)
    elif args.FUNCTION == 'mokapot':
        from run_mokapot import run_mokapot 
        # Mokapot task should publish its directory, since mokapot sprays its files around in the working dir
        if isinstance(args.psm, str):
            run_mokapot(args.psm)
        else:
            [run_mokapot(x) for x in args.psm]
    elif args.FUNCTION == 'dino_annotate':
        from dinosaur_quantitation import annotate_psms_with_features
        annotate_psms_with_features(args.psm, args.data)
    elif args.FUNCTION == 'psm_report':
        from generate_final_report import generate_PSM_report
        generate_PSM_report(args.pin, args.psm[0])
    else:
        raise RuntimeError("Not a valid function: %s" % args.FUNCTION)



if __name__ == "__main__":
    run_pipeline_func()

