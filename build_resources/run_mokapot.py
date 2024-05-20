import os
import pandas as pd
import mokapot
import numpy as np

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

    # if outputfile is None:
    #     outputfile = os.path.join(os.path.dirname(psm_files[0]), 'concatenated.pin')
    # full.to_csv(outputfile, index=False, sep='\t')
    # 
    # print("Done concatenating")
    # return outputfile

    print("Done concatenating")
    return full.copy()



# def excerpt_columns_for_mokapot(psm_file, outputfile=None):
#     psm = pd.read_csv(psm_file, sep='\t')
#     psm = psm[['scan', 'peptide', 'charge', 'mz', 'retention_time', 'protein', 'decoy', 'percolator_score']]
#     if outputfile is None:
#         outputfile = os.path.join(os.path.dirname(psm_file), 'mokapot_input.txt')
#     psm.to_csv(outputfile, index=False, sep='\t')
#     return outputfile


def add_imputed_features(psms, outputfile):
    # TODO rationalize where this sort of thing occurs vs the annotatation script
    # psms = pd.read_csv(pin_file, sep='\t')

    # for z in [2, 3, 4, 5, 6]:
    #     psms['charge_%d' % z] = (psms['charge'] == z).astype(int)

    psms['log_frag_count'] = psms['matched_fragments'].apply(lambda x: 0 if x == 0 else np.log(x))
    psms['log_isotope_error'] = psms['isotope_ratio_error'].apply(lambda x: 0 if x == 0 else np.log(abs(x)))
    psms['log_rt_error'] = psms['abs_rt_error'].apply(lambda x: 0 if x == 0 else np.log(x))
    psms['log_fragment_error'] = psms['fragment_error'].apply(lambda x: 0 if x == 0 else np.log(abs(x)))
    psms['log_frag_product'] = (psms['matched_fragments'] * psms['fragment_error']).apply(lambda x: 0 if x == 0 else np.log(abs(x)))
    psms['log_frag_div'] = (psms['fragment_error'] / psms['matched_fragments'].apply(lambda x: x if x else 1)).apply(lambda x: 0 if x == 0 else np.log(abs(x)))
    psms['log_frag_dist'] = psms['fragment_cosine_dist'].apply(lambda x: x if x else 1).apply(lambda x: 0 if x == 0 else np.log(abs(x)))

    psms['is_missed_cleave'] = psms['peptide'].apply(lambda x: 'K' in x[:-1] or 'R' in x[:-1]).astype(int)
    psms['length'] = psms['peptide'].apply(len)
    psms['mod_count'] = psms['peptide'].apply(lambda x: x.count('['))
    
    # psms.to_csv(outputfile, index=False, sep='\t')
    return psms.copy()





def run_mokapot(pin_data, fasta_file, output_dir=None):
    print("Running mokapot")
    # mokamod = mokapot.Model(estimator=GradientBoostingClassifier())

    moka = mokapot.read_pin(pin_data, filename_column='filename', calcmass_column='calcmass',
                            expmass_column='expmass', charge_column='charge', rt_column='rt')
    moka.add_proteins(fasta_file, decoy_prefix='rev_', min_length=5)
    confidences, models = mokapot.brew(moka)

    print(confidences)
    print('\n'.join([f"{x[0]}: {x[1]}" for x in zip(models[0].features, models[0].estimator.coef_.ravel())]))

    confidences.to_txt(dest_dir=output_dir, decoys=False)
    print("Done running mokapot")          




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser() 
    parser.add_argument('-f', '--fasta', required=True)
    parser.add_argument('-p', '--psm', nargs='+', required=True)
    args = parser.parse_args()

    psm_files = args.psm
    fasta_file = args.fasta

    table = concatenate_psm_files(psm_files)
    for psm_file in psm_files:
        os.remove(psm_file)
    table.to_csv('concatenated.pin', index=False, sep='\t')
    print("Done concatenating")
    table = add_imputed_features(table, None)
    run_mokapot(table, fasta_file, output_dir=os.path.dirname(psm_files[0]))
