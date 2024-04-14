import os
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import re




def _comet_mods_to_reconciled(modstr):
    if modstr == '-':
        return ''
    else:
        agg = []
        for mod in modstr.split(','):
            site, _, mass = mod.split('_', 2)
            mass = mass.lower()
            if '_n' in mass:
                mass = mass.replace('_n', '')
                site = '0'
            if '42.' in mass:
                assert(site=='0')
            agg.append((site, float(mass)))
        return ','.join(['%s_%f' % x for x in sorted(agg)])

# peptide_re = re.compile(r'([A-Z]*)([+-][0-9.]+)+')
def _sage_mods_to_reconciled(pepstr):
    agg = []
    loc_count = 0
    for subseq, modmass in re.findall(r'([A-Z]*)\[[+-]([0-9.]+)\]', pepstr):
        loc_count += len(subseq)
        if '42.' in modmass:
             assert(loc_count==0)
        agg.append((loc_count, float(modmass)))
    
    return ','.join(['%s_%f' % x for x in sorted(agg)])



# Neither Sage nor Comet give us the charge state directly. We have to infer it from the one-hot encoding.
def _deconvert_sage_onehot_chg(psm):
    if not pd.isnull(psm['sage_z=2']):
        chg_series = psm[['sage_z=2', 'sage_z=3', 'sage_z=4', 'sage_z=5', 'sage_z=6', 'sage_z=other']]
    else:
        chg_series = psm[['wide_z=2', 'wide_z=3', 'wide_z=4', 'wide_z=5', 'wide_z=6', 'wide_z=other']]
    hot_chg = chg_series.index[chg_series.apply(bool)][0]
    try:
        return int(hot_chg[-1])
    except:
        return 7  # Most likely. Chg 1 is freaky, chg 8 is basically the same as 7, right?
    
def _deconvert_comet_onehot_chg(psm):
    chg_series = psm[['comet_charge1', 'comet_charge2', 'comet_charge3', 'comet_charge4', 'comet_charge5', 'comet_charge6']]
    hot_chg = chg_series.index[chg_series.apply(bool)][0]
    return int(hot_chg[-1])
    # Comet is at least nice enough not to give us that "other" value.





def impute_hyperscore(psms):
    # Imputes hyperscore values, when missing, from comet_xcorr and comet_sp_score. Probably valid in that
    # the hyperscore is only going to be used as a discriminant value; adding imputations into the mix makes it
    # *less* useful in discriminating between good and bad PSMs, but it's better than having no such value at all.

    has_hyperscore_and_other = psms.hyperscore.notnull() & psms.comet_xcorr.notnull() & psms.comet_sp_score.notnull()
    missing_hyperscore = psms.hyperscore.isnull()

    psms['hyperscore_imputed'] = 1
    psms.loc[psms.hyperscore.notnull(), 'hyperscore_type'] = 0

    reg = LinearRegression()
    reg.fit(psms.loc[has_hyperscore_and_other, ['comet_xcorr', 'comet_sp_score']], psms.loc[has_hyperscore_and_other, 'hyperscore'])
    score = reg.score(psms.loc[has_hyperscore_and_other, ['comet_xcorr', 'comet_sp_score']], psms.loc[has_hyperscore_and_other, 'hyperscore'])
    print("Imputing with score of %s" % score)

    psms.loc[missing_hyperscore, 'hyperscore'] = reg.predict(psms.loc[missing_hyperscore, ['comet_xcorr', 'comet_sp_score']])

    return psms




def civilized_pin(pinfile):
    outfile = pinfile.replace('.pin', '.civ.pin')
    assert(outfile!=pinfile)
    agg = []
    with open(pinfile, 'r') as inp:
        header = inp.__next__().strip().split('\t')
        for line in inp:
            values = line.strip().split('\t', len(header)-1)

            agg.append(values)
    foo = pd.DataFrame(agg, columns=header)
    foo.to_csv(outfile, index=False)
    return outfile


def _pref_order(*things):
    try:
        return [x for x in things if not pd.isnull(x)][0]
    except IndexError:
        return np.nan


sage_tsv_columns = ['peptide', 'proteins', 'num_proteins', 'filename', 'scannr', 'rank', 'label', 'expmass', 'calcmass', 'charge', 
                    'peptide_len', 'missed_cleavages', 'isotope_error', 'precursor_ppm', 'fragment_ppm', 'hyperscore', 'delta_next',
                    'delta_best', 'rt', 'aligned_rt', 'predicted_rt', 'delta_rt_model', 'matched_peaks', 'longest_b', 'longest_y',
                    'longest_y_pct', 'matched_intensity_pct', 'scored_candidates', 'poisson', 'sage_discriminant_score',
                    'posterior_error', 'spectrum_q', 'peptide_q', 'protein_q', 'ms1_intensity', 'ms2_intensity']
comet_pin_columns = ['SpecId', 'Label', 'ScanNr', 'ExpMass', 'CalcMass', 'lnrSp', 'deltLCn', 'deltCn', 'lnExpect', 'Xcorr', 'Sp', 
                     'IonFrac', 'Mass', 'PepLen', 'Charge1', 'Charge2', 'Charge3', 'Charge4', 'Charge5', 'Charge6', 'enzN', 'enzC', 
                     'enzInt', 'lnNumSP', 'dM', 'absdM', 'Peptide', 'Proteins']
comet_tsv_columns = ['scan', 'num', 'charge', 'exp_neutral_mass', 'calc_neutral_mass', 'e-value', 'xcorr', 'delta_cn', 'sp_score',
                      'ions_matched', 'ions_total', 'plain_peptide', 'modified_peptide', 'prev_aa', 'next_aa', 'protein', 'protein_count', 
                      'modifications', 'retention_time_sec', 'sp_rank']
def reconcile_search_results(comet_pin_file, comet_txt_file, sage_pin_file, widesage_pin_file, output_file=None):
    print("Reconciling %s" % [comet_pin_file, comet_txt_file, sage_pin_file])

    cometpin = pd.read_csv(civilized_pin(comet_pin_file))
    comettxt = pd.read_csv(comet_txt_file, sep='\t', skiprows=1)
    cometpin['psm_num'] = cometpin.SpecId.apply(lambda x: int(x.split('_')[-1]))
    comet_combined = cometpin.merge(comettxt, left_on=['ScanNr', 'psm_num'], right_on=['scan', 'num'])
    assert(len(comet_combined)==len(cometpin)), len(comet_combined)
    assert(len(comet_combined)==len(comettxt)), len(comet_combined)
    comet_combined.columns = [x.lower() for x in comet_combined.columns]
    comet_combined = comet_combined.loc[:,~comet_combined.columns.duplicated()].copy() # type: ignore

    sagepin = pd.read_csv(sage_pin_file, sep='\t')
    sagepin.columns = [x.lower() for x in sagepin.columns]
    widesagepin = pd.read_csv(widesage_pin_file, sep='\t')
    widesagepin.columns = [x.lower() for x in widesagepin.columns]


    comet_combined = comet_combined.set_index('scannr')
    sagepin = sagepin.set_index('scannr')
    widesagepin = widesagepin.set_index('scannr')

    scan_numbers = set(list(comet_combined.index)+list(sagepin.index)+list(widesagepin.index))

    comet_combined['rec_mods'] = comet_combined.modifications.apply(_comet_mods_to_reconciled)
    sagepin['rec_mods'] = sagepin.peptide.apply(_sage_mods_to_reconciled)
    sagepin['plain_peptide'] = sagepin.peptide.apply(lambda x: ''.join([c for c in x if c.isupper()]))
    widesagepin['rec_mods'] = widesagepin.peptide.apply(_sage_mods_to_reconciled)
    widesagepin['plain_peptide'] = widesagepin.peptide.apply(lambda x: ''.join([c for c in x if c.isupper()]))

    # Used as default values of search rank when the engine in question did not find the match
    sage_default_rank = sagepin['rank'].max()+2
    widesage_default_rank = widesagepin['rank'].max()+2
    comet_default_rank = comettxt.num.max()+2

    sagepin = sagepin.rename(columns={x:'sage_'+x for x in sagepin.columns if x not in {'rec_mods', 'plain_peptide'}})
    widesagepin = widesagepin.rename(columns={x:'wide_'+x for x in widesagepin.columns if x not in {'rec_mods', 'plain_peptide'}})
    comet_combined = comet_combined.rename(columns={x:'comet_'+x for x in comet_combined.columns if x not in {'rec_mods', 'plain_peptide'}})

    file_stem = os.path.basename(comet_pin_file).split('.')[0]

    reconciled_rows = []
    reconciled_columns = ['label', 'filename', 'scannr', 'expmass', 'calcmass', 'peptide', 'rec_mods', 'proteins',
                          'comet_rank', 'sage_rank', 'wide_rank', 'in_both',
                          'charge', 'num_proteins', 'peptide_len', 'missed_cleavages',
                          'delta_mass', 'abs_delta_mass',
                          'hyperscore', 'comet_xcorr', 'comet_sp_score']
    for scan_num in scan_numbers:
        # Comet and Sage both may be supplying multiple results per scan; matching results from each should
        # each be reconciled.
        com_psms = comet_combined.loc[[scan_num]] if scan_num in comet_combined.index else pd.DataFrame(columns=comet_combined.columns)
        sage_psms = sagepin.loc[[scan_num]] if scan_num in sagepin.index else pd.DataFrame(columns=sagepin.columns)
        widesage_psms = widesagepin.loc[[scan_num]] if scan_num in widesagepin.index else pd.DataFrame(columns=widesagepin.columns)

        scan_psms = com_psms.merge(sage_psms, how='outer', 
                                   left_on=['plain_peptide', 'rec_mods'],
                                   right_on=['plain_peptide', 'rec_mods']
                                   ).merge(widesage_psms, how='outer',
                                   left_on=['plain_peptide', 'rec_mods'],
                                   right_on=['plain_peptide', 'rec_mods']
                                           )


        for _, psm in scan_psms.iterrows():
            is_sage = not bool(pd.isnull(psm['sage_peptide']))
            is_comet = not bool(pd.isnull(psm['comet_peptide']))
            is_widesage = not bool(pd.isnull(psm['wide_peptide']))

            if is_sage or is_widesage:
                chg = _deconvert_sage_onehot_chg(psm)
            else:
                chg = _deconvert_comet_onehot_chg(psm)
            
            row = {
                    'sage_rank': psm['sage_rank'] if is_sage else sage_default_rank,
                    'comet_rank': psm['comet_num'] if is_comet else comet_default_rank,
                    'wide_rank': psm['wide_rank'] if is_widesage else widesage_default_rank,

                    'in_both': int(is_sage and is_comet),

                    'filename': file_stem,
                    'scannr': scan_num,
                    'label': _pref_order(psm['sage_label'], psm['comet_label'], psm['wide_label']),
                    'expmass': _pref_order(psm['sage_expmass'], psm['comet_expmass'], psm['wide_expmass']),
                    'calcmass': _pref_order(psm['sage_calcmass'], psm['comet_calcmass'], psm['wide_calcmass']),
                    'peptide': _pref_order(psm['plain_peptide'], psm['plain_peptide'], psm['plain_peptide']),
                    'proteins': _pref_order(psm['sage_proteins'], psm['comet_proteins'], psm['wide_proteins']),
                    'rec_mods': _pref_order(psm['rec_mods'], psm['rec_mods'], psm['rec_mods']),
                    'charge': chg,

                    'num_proteins': _pref_order(psm['comet_protein_count'], str(psm['sage_proteins']).count('\t')+1, str(psm['wide_proteins']).count('\t')+1),
                    'peptide_len': len(psm['plain_peptide']),
                    'missed_cleavages': _pref_order(psm['sage_missed_cleavages'], psm['wide_missed_cleavages'], str(psm['plain_peptide']).count('R')+str(psm['plain_peptide']).count('K')-1),
                    'delta_mass': (_pref_order(psm['sage_expmass'], psm['comet_expmass'], psm['wide_expmass']) - _pref_order(psm['sage_calcmass'], psm['comet_calcmass'], psm['wide_calcmass'])),
                    
                    'hyperscore': _pref_order(np.exp(psm['sage_ln(hyperscore)']), np.exp(psm['wide_ln(hyperscore)'])),  # Could be maximum?
                    'comet_xcorr': psm['comet_xcorr'],
                    'comet_sp_score': psm['comet_sp_score'],
                    }
            row['abs_delta_mass'] = abs(row['delta_mass'])

            reconciled_rows.append(row)

    assert(set(reconciled_rows[0].keys())==set(reconciled_columns)), set(reconciled_rows[0].keys())^set(reconciled_columns)
    reconciled = pd.DataFrame(reconciled_rows, columns=reconciled_columns)

    print(reconciled.shape)
    print(reconciled[reconciled.in_both==1].shape)
    print(reconciled[reconciled.hyperscore.isnull()].shape)
    print(reconciled[reconciled.comet_xcorr.isnull()].shape)
    print(reconciled.sage_rank.value_counts())
    print(reconciled.comet_rank.value_counts())
    print(reconciled.wide_rank.value_counts())
    reconciled = impute_hyperscore(reconciled).drop(['comet_xcorr', 'comet_sp_score'], axis=1)

    # We keep a lot of redundancy for the sake of getting a more robust flase hit distribution, but
    # there's no reason to keep multiple instances of the same PSM.
    reconciled.sort_values(by=['hyperscore'], inplace=True, ascending=False)
    reconciled.drop_duplicates(subset=['filename', 'scannr', 'peptide', 'rec_mods', 'charge'], keep='first', inplace=True)

    for i in range(1, 7):
        reconciled['chg_is_%d' % i] = (reconciled.charge==i).apply(int)
    
    if output_file is None:
        output_file = comet_pin_file.replace('.pin', '.reconciled.pin')
    reconciled.to_csv(output_file, index=False)
    return output_file



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Reconcile Sage and Comet search results')
    parser.add_argument('comet_pin', type=str, help='Comet pin file')
    parser.add_argument('comet_txt', type=str, help='Comet txt file')
    parser.add_argument('sage_pin', type=str, help='Sage pin file')
    parser.add_argument('widesage_pin', type=str, help='Wide-window Sage pin file')
    parser.add_argument('-o', '--output', type=str, help='Output file')
    args = parser.parse_args()

    reconcile_search_results(args.comet_pin, args.comet_txt, args.sage_pin, args.widesage_pin, args.output)
