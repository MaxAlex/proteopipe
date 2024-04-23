import os
import pandas as pd
import mokapot


def concatenate_psm_files(psm_files, exp_name=None, outputfile=None):
    print("Concatenating %s" % psm_files)
    agg = []
    for psm_file in psm_files:
        try:
            psm = pd.read_csv(psm_file)
            assert(len(psm.columns)>3)
        except (AssertionError, pd.errors.ParserError):
            psm = pd.read_csv(psm_file, sep='\t')

        agg.append(psm)
    
    full = pd.concat(agg)

    # bar.dtypes[bar.dtypes=='object'].index
    # Have to ensure there's no tabs in any text field, otherwise mokapot's
    # homespun parser will barf.
    for col in full.dtypes[full.dtypes=='object'].index:
        full[col] = full[col].apply(lambda x: x.replace('\t', ' '))

    if outputfile is None:
        outputfile = os.path.join(os.path.dirname(psm_files[0]), 'concatenated.pin')
    full.to_csv(outputfile, index=False, sep='\t')
    
    print("Done concatenating")
    return outputfile



# def excerpt_columns_for_mokapot(psm_file, outputfile=None):
#     psm = pd.read_csv(psm_file, sep='\t')
#     psm = psm[['scan', 'peptide', 'charge', 'mz', 'retention_time', 'protein', 'decoy', 'percolator_score']]
#     if outputfile is None:
#         outputfile = os.path.join(os.path.dirname(psm_file), 'mokapot_input.txt')
#     psm.to_csv(outputfile, index=False, sep='\t')
#     return outputfile

def run_mokapot(pin_file, output_stem=None):
    print("Running mokapot")
    # mokamod = mokapot.Model(estimator=GradientBoostingClassifier())
    mokamod = None

    moka = mokapot.read_pin(pin_file, filename_column='filename', calcmass_column='calcmass',
                            expmass_column='expmass', charge_column='charge', rt_column='rt')
    confidences, models = mokapot.brew(moka)
    if output_stem:
        confidences.to_txt(file_root=output_stem, decoys=False)
        print("Mokapot files given output stem %s" % output_stem)
        return
    else:
        files = confidences.to_txt(os.path.dirname(pin_file), decoys=False)
        print(files)
        print("Mokapot completed")
        return files[0]




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser() 
    parser.add_argument('-p', '--psm', nargs='+', required=True)
    args = parser.parse_args()

    psm_files = args.psm

    concat_file = concatenate_psm_files(psm_files)
    run_mokapot(concat_file)
