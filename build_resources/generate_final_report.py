import pandas as pd
from mod_mass import identify_mod_by_mass
import re


# def render_nice_peptide_strings(peptide):
#     if '___' not in peptide:
#         return peptide, ''
#     else:
#         peptide, modstr = peptide.split('___')
#         mods = []
#         for mod in modstr.split(','):
#             loc, mass = mod.split('_')
#             loc = int(loc)
#             aa = peptide[loc-1]
#             modname = identify_mod_by_mass(float(mass))
#             mods.append('%s%d: %s' % (aa, loc, modname))
#         return peptide, '; '.join(mods)

def render_nice_peptide_string(peptide):
    offset = 0
    mod_dict = {}
    for start, stop in [(m.start(), m.end()) for m in re.finditer(r'\[\d+\.\d+\]', peptide)]:
        mod_dict[start-offset] = peptide[start:stop].strip('[]')
        offset += stop-start
    
    peptide = re.sub(r'\[\d+\.\d+\]', '', peptide)
    print((len(peptide), mod_dict))
    mod_strs = []
    for loc, name in mod_dict.items():
        if loc==0:
            mod_strs.append(f'N-term: {name}')
        elif loc==len(peptide)+1:
            mod_strs.append(f'C-term: {name}')
        else:
            mod_strs.append('%s%d: %s' % (peptide[loc-1], loc, name))
    return peptide, '; '.join(mod_strs)


COLUMN_ORDER = ['Spectrum_ID', 'Filename', 'Scan', 'Peptide', 'Modifications', 'ModifiedPeptide', 'Charge', 'M/Z',
                'Calculated Mass', 'Experimental Mass', 'Retention Time', 'dMass', 
                'Search Engine', 'Mokapot q-value', 'Mokapot score', 'Mokapot PEP',
                'RT Prediction Error', 'Fragmentation Prediction Error', 'Fragment Match Count', 'Isotope Profile Error',
                'Proteins']

def generate_PSM_report(pin_file, mokapot_psm_file):
    pin = pd.read_csv(pin_file, sep='\t')
    mokap = pd.read_csv(mokapot_psm_file, sep='\t')

    # mokap = mokap.sort_values(by=['mokapot score']).drop_duplicates(subset=['specid'], keep='last').copy() # type: ignore
    moka_best_inds = mokap.groupby('specid')['mokapot score'].idxmax()
    mokap = mokap.loc[moka_best_inds].copy()

    # pin_best_inds = pin.groupby(['res_specid', 'res_peptide'])['res_hyperscore'].idxmax()
    # pin = pin.loc[pin_best_inds].copy()
    
    combined = pd.merge(pin.rename({x:'res_' + x for x in pin.columns}, axis=1),
                        mokap.rename({x:'mok_' + x for x in mokap.columns}, axis=1),
                        how='right', left_on=['res_specid', 'res_peptide'], right_on=['mok_specid', 'mok_peptide'])

    print((pin.shape, mokap.shape, combined.shape))
    
    # Technically it is arbitrary which redundant rows are chosen, the peptide is the same across all
    comb_inds = combined.groupby('res_specid')['res_hyperscore'].idxmax()
    combined = combined.loc[comb_inds].copy()

    assert(mokap['specid'].nunique() == mokap.shape[0])
    assert(combined['res_specid'].nunique() == combined.shape[0])
    assert(combined['mok_specid'].nunique() == combined.shape[0])

    comet_default_rank = combined.res_comet_rank.max()
    sage_default_rank = combined.res_sage_rank.max()
    wide_default_rank = combined.res_wide_rank.max()
    
    new_rows = []
    for _, row in combined.iterrows():
        peptide, modifications = render_nice_peptide_string(row['res_peptide'])
       
        engine = '|'.join([eng for eng, name, default in 
                           [('Comet', 'comet', comet_default_rank),
                            ('Sage', 'sage', sage_default_rank), 
                            ('Sage(Wide-Tolerance)', 'wide', wide_default_rank)]
                           if row['res_%s_rank' % name.lower()] != default])

        new_row = {'Spectrum_ID': row['res_specid'],
                   'Filename': row['res_filename'],
                   'Scan': row['res_scannr'],
                   'Peptide': peptide,
                   'Modifications': modifications,
                   'ModifiedPeptide': row['res_peptide'], 
                   'Charge': row['res_charge'],
                   'M/Z': row['res_mz'],
                   'Calculated Mass': row['res_calcmass'],
                   'Experimental Mass': row['res_expmass'],
                   'Retention Time': row['res_rt'],
                   'Delta Mass': row['res_delta_mass'],
                   'Search Engine': engine,
                   'Mokapot q-value': row['mok_mokapot q-value'],
                   'Mokapot score': row['mok_mokapot score'],
                   'Mokapot PEP': row['mok_mokapot PEP'],
                   'RT Prediction Error': row['res_rt_error'],
                   'Fragmentation Prediction Error': row['res_fragment_error'], 
                   'Fragment Match Count': row['res_matched_fragments'],
                   'Isotope Profile Error': row['res_isotope_ratio_error'],
                   'Proteins': row['res_proteins'],
                   }
        
        new_rows.append(new_row)

    assert(set(COLUMN_ORDER) == set(new_rows[0].keys())), (set(COLUMN_ORDER) ^ set(new_rows[0].keys()))

    report = pd.DataFrame(new_rows, columns=COLUMN_ORDER)
    report.to_csv(pin_file.split('.')[0] + '.report.csv', index=False)


if __name__ == '__main__':
    # import sys
    # generate_PSM_report(sys.argv[1], sys.argv[2])
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pin')
    parser.add_argument('-m', '--mokapot')
    args = parser.parse_args()
    generate_PSM_report(args.pin, args.mokapot)
