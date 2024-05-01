import pandas as pd
import numpy as np
from collections import defaultdict
import re


# np.seterr(divide='raise')


def ms1_quant_normalize(pep_table):
    # Corrects for relative absolute abundance of precursor in different runs; originating
    # from different levels of sample in each injection, varying performance of the injection,
    # or etc.
    # This could be moved to the post-search PSM aggregation step

    apex_ints = [x for x in pep_table.columns if 'Sum Apex Int' in x]
    total_ints = [x for x in pep_table.columns if 'Sum Total Int' in x]
    
    apex_medians = pep_table[apex_ints].median()
    total_medians = pep_table[total_ints].median()

    pep_table[apex_ints] = pep_table[apex_ints].divide(apex_medians, axis=1)
    pep_table[total_ints] = pep_table[total_ints].divide(total_medians, axis=1)
   
    print("Apex medians %s" % apex_medians)
    print("Total medians %s" % total_medians)

    return pep_table


# Assuming that all relevant fractions are represented in the output! If somehow some fraction gave
# no PSMs, it won't be in the table (but that would be a super pathological case for other reasons.)
def agg_psms_to_peptides(psm_file, mokapot_peptides, run_parse_regexp, norm_to_isochannel=None):
    psms = pd.read_csv(psm_file)
    moka_pep = pd.read_csv(mokapot_peptides, sep='\t')
    moka_pep = moka_pep.set_index('peptide') 
    assert(moka_pep.label.all())

    # Add a column for the fraction
    psms['Run_ID'] = psms['Filename'].apply(lambda x: re.search(run_parse_regexp, x).group(1))

    runs = sorted(set(psms['Run_ID']))
    # tmt_cols = [x for x in psms.columns if x.startswith('TMT_')]
    tmt_cols = [x for x in psms.columns if x.startswith("1")]
    if tmt_cols and norm_to_isochannel is None:
        norm_to_isochannel = tmt_cols[0]
    pep_rows = []
    for pep_str, pep_psms in psms.groupby('ModifiedPeptide'):
        try:
            mok_pep = moka_pep.loc[pep_str] 
        except KeyError:
            print("Skipping %s, not in mokapot output" % pep_str)
            continue
       
        pep_iso_quants = dict(sum([[("Lg2Ratio_TMT_%s_%s" % (run, channel), np.nan) for channel in tmt_cols] for run in runs], []))
        pep_iso_direct = dict(sum([[("TMT_%s_%s" % (run, channel), np.nan) for channel in tmt_cols] for run in runs], []))
        pep_iso_magnitudes = {} 
        pep_lfree_apexints = {}
        pep_lfree_totalints = {}
        run_pep = {}
        run_qval = {}
        run_score = {}
        for run, run_psms in pep_psms.groupby('Run_ID'):
            run_quant = run_psms[tmt_cols].sum() + 1e-6 
            norm_run_quant = np.log2(run_quant / run_quant[norm_to_isochannel]
                                     ).rename(lambda x: 'Lg2Ratio_TMT_%s_%s' % (run, x))

            if norm_run_quant.max() >= np.inf:
                raise Exception

            pep_iso_magnitudes[run] = run_quant.sum()
            pep_iso_quants.update(norm_run_quant)
            pep_iso_direct.update(run_quant.rename(lambda x: 'TMT_%s_%s' % (run, x)))

            pep_lfree_apexints[run] = run_psms['Feature Apex Int'].sum()
            pep_lfree_totalints[run] = run_psms['Feature Total Int'].sum()

            # Being conscientious here by taking only the minimum PEP/q-value and maximum score, on the basis
            # that the same PSM of the same peptide in a given run isn't an iid observation. Different mod states
            # a more iid-y but that's handled at the protein level
            run_pep[run] = run_psms['Mokapot PEP'].min()
            run_qval[run] = run_psms['Mokapot q-value'].min()
            run_score[run] = run_psms['Mokapot score'].max()

        pep_iso_quants = pd.Series(pep_iso_quants)
        pep_iso_direct = pd.Series(pep_iso_direct)
        pep_lfree_apexints = pd.Series([pep_lfree_apexints.get(k, np.nan) for k in runs],
                                       index=['%s Sum Apex Int' % x for x in runs])
        pep_lfree_totalints = pd.Series([pep_lfree_totalints.get(k, np.nan) for k in runs],
                                        index=['%s Sum Total Int' % x for x in runs])
        pep_lfree_quants = pd.concat([pep_lfree_apexints, pep_lfree_totalints])
        iso_magnitudes = pd.Series([pep_iso_magnitudes.get(k, np.nan) for k in runs],
                                   index=['%s Label Intensity' % k for k in runs])
        pep_peps = [('%s PEP' % k, run_pep.get(k, '-')) for k in runs]
        pep_qvals = [('%s q-value' % k, run_qval.get(k, '-')) for k in runs]
        pep_scores = [('%s score' % k, run_score.get(k, '-')) for k in runs]
        pep_peps, pep_qvals, pep_scores = list(map(lambda xs: pd.Series([x[1] for x in xs], 
                                                                        index=[x[0] for x in xs]),
                                                   [pep_peps, pep_qvals, pep_scores]))

        pep_row = pep_psms.iloc[0][['Peptide', 'Modifications', 'Calculated Mass', 'Proteins']].copy()
        pep_row['Runs'] = ';'.join(sorted(set(pep_psms['Run_ID'])))
        pep_row['Charges'] = ';'.join(sorted(set(pep_psms['Charge'].astype(str))))
        pep_row['PSMs'] = len(pep_psms)

        # Here I'm trusting mokapot's PEP/q-value/score over the PSM-level ones without much justification
        # Likely some guy thought hard about how to do it just right... but possibly not!
        pep_row['PEP'] = mok_pep['mokapot PEP']
        pep_row['q-value'] = mok_pep['mokapot q-value']
        pep_row['Score'] = mok_pep['mokapot score']
        pep_row = pd.concat([pep_row, pep_peps, pep_qvals, pep_scores,
                             iso_magnitudes, pep_lfree_quants, pep_iso_quants, pep_iso_direct])
        pep_rows.append(pep_row)


    # expected_columns = ['Protein', 'Peptides', 'Unique Peptides', 'PSMs', 'PEP', 'q-value', 'Score'] + \
    #         sum([['LgTMT_%s_%s' % (run, channel) for channel in tmt_cols] for run in runs], []) + \
    #         sum([['TMT_%s_%s' % (run, channel) for channel in tmt_cols] for run in runs], []) + \
    #         sum([['%s Sum Apex Int' % run, '%s Sum Total Int' % run, '%s Label Intensity' % run] for run in runs], []) \

    expected_columns = ['Peptide', 'Modifications', 'Calculated Mass', 'Proteins', 'Runs', 'Charges', 'PSMs', 'PEP', 'q-value', 'Score'] + \
            sum([['%s PEP' % run, '%s q-value' % run, '%s score' % run, '%s Label Intensity' % run, '%s Sum Apex Int' % run, '%s Sum Total Int' % run] for run in runs], []) + \
            sum([['Lg2Ratio_TMT_%s_%s' % (run, channel) for channel in tmt_cols] for run in runs], []) + \
            sum([['TMT_%s_%s' % (run, channel) for channel in tmt_cols] for run in runs], [])
    assert(set(pep_rows[0].index).issubset(expected_columns))
    pep_table = pd.DataFrame(pep_rows, columns=expected_columns)

    # pep_table = tmt_channel_correction(pep_table, runs)
    pep_table = ms1_quant_normalize(pep_table)

    pep_table.to_csv(psm_file + '.peptides.csv', index=False)

    return psm_file + '.peptides.csv'


def agg_peptides_to_protein(pep_table):
    peptides = pd.read_csv(pep_table)
    runs = sorted(set([x.split()[0] for x in peptides.columns if len(x.split()) > 1 and x.split()[1] == 'PEP']))
    tmt_cols = sorted(set([x.split('_')[-1] for x in peptides.columns if x.startswith('TMT_')]))
    assert(runs)
    assert(tmt_cols)
    print(runs)
    print(tmt_cols)

    # This will use the ENSP accessions, which aren't directly useful to anyone. Could do [-2] to get the gene name
    # which people would recognize, but the results won't be *protein* per row per se.
    peptides['prot_accessions'] = peptides['Proteins'].apply(lambda x: ([x.split('|')[0] for x in x.split(';')]))
    
    baseTMT_cols = ['TMT_%s_%s' % (run, channel) for run in runs for channel in tmt_cols]
    lgTMT_cols = ['Lg2Ratio_TMT_%s_%s' % (run, channel) for run in runs for channel in tmt_cols]

    lfree_ints = ['%s Label Intensity' % run for run in runs]
    lfree_apexints = ['%s Sum Apex Int' % run for run in runs]
    lfree_totalints = ['%s Sum Total Int' % run for run in runs]

    # assert(not peptides[label_cols].isnull().all().all())

    protein_peptides = defaultdict(list)
    for _, pep in peptides.iterrows():
        for prot in pep['prot_accessions']:
            protein_peptides[prot].append(pep)
    protein_peptides = {k: pd.DataFrame(v) for k, v in protein_peptides.items()}

    prot_rows = []
    for protein, subpeptides in protein_peptides.items():
        prot_row = {'Protein': protein,
                    'Peptides': ', '.join(subpeptides.apply(lambda x: '%s (%s)' % (x['Peptide'], x['Modifications']), axis=1)),
                    'Unique Peptides': ', '.join(subpeptides[subpeptides['prot_accessions'].apply(len)==1].apply(lambda x: '%s (%s)' % (x['Peptide'], x['Modifications']), axis=1)),
                    'PSMs': subpeptides['PSMs'].sum(),
                    'PEP': subpeptides['PEP'].min(),
                    'q-value': subpeptides['q-value'].min(),
                    'Score': subpeptides['Score'].max(),
                    # 'Runs': ';'.join(sorted(set(sum([x['Runs'].split(';') for x in peptides], [])))), # type: ignore
                    }

        tmt_max_vals = subpeptides[baseTMT_cols].max(skipna=True).rename(lambda x: 'Max_' + x)
        tmt_med_vals = subpeptides[baseTMT_cols].median(skipna=True).rename(lambda x: 'Med_' + x)
        tmtrat_med_vals = subpeptides[lgTMT_cols].median(skipna=True).rename(lambda x: 'Med_' + x)

        lfree_max_int = subpeptides[lfree_ints].max(skipna=True).rename(lambda x: 'Max_' + x) 
        lfree_med_int = subpeptides[lfree_ints].median(skipna=True).rename(lambda x: 'Med_' + x)
        lfree_max_apexint = subpeptides[lfree_apexints].max(skipna=True).rename(lambda x: 'Max_' + x)
        lfree_max_totalint = subpeptides[lfree_totalints].max(skipna=True).rename(lambda x: 'Max_' + x)

        quants = pd.concat([tmt_max_vals, tmt_med_vals, tmtrat_med_vals, lfree_max_int, lfree_med_int, lfree_max_apexint, lfree_max_totalint])




        prot_row = pd.Series(prot_row)
        prot_row = pd.concat([prot_row, quants])
        prot_rows.append(prot_row)

    prot_table = pd.DataFrame(prot_rows)
    prot_table.to_csv(pep_table + '.proteins.csv', index=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--psm', required=True)
    parser.add_argument('-m', '--mokapot', required=True)
    parser.add_argument('--mokapot_prots', required=True)
    parser.add_argument('-r', '--run_parse_regexp', required=True)
    parser.add_argument('-n', '--norm_to_isochannel')
    args = parser.parse_args()

    run_parse = args.run_parse_regexp.strip('"')
    
    pep_table = agg_psms_to_peptides(args.psm, args.mokapot, args.run_parse_regexp, args.norm_to_isochannel)
    agg_peptides_to_protein(pep_table)

# if __name__ == '__main__':
#     import sys
#     agg_peptides_to_protein(sys.argv[1])

         
        


        
        


                 
    
