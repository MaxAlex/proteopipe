import pandas as pd
from mod_mass import identify_mod_by_mass



def render_nice_peptide_strings(peptide):
    if '___' not in peptide:
        return peptide, ''
    else:
        peptide, modstr = peptide.split('___')
        mods = []
        for mod in modstr.split(','):
            loc, mass = mod.split('_')
            loc = int(loc)
            aa = peptide[loc-1]
            modname = identify_mod_by_mass(float(mass))
            mods.append('%s%d: %s' % (aa, loc, modname))
        return peptide, '; '.join(mods)

def generate_PSM_report(pin_file, mokapot_psm_file):
    pin = pd.read_csv(pin_file, sep='\t')
    mokap = pd.read_csv(mokapot_psm_file, sep='\t')

    mokap = mokap.sort_values(by=['mokapot score']).drop_duplicates(subset=['specid'], keep='last').copy() # type: ignore

    combined = pd.merge(pin.rename({x:'res_' + x for x in pin.columns}, axis=1),
                        mokap.rename({x:'mok_' + x for x in mokap.columns}, axis=1),
                        how='inner', left_on=['res_specid', 'res_peptide'], right_on=['mok_specid', 'mok_peptide'])

    print((pin.shape, mokap.shape, combined.shape))

    comet_default_rank = combined.res_comet_rank.max()
    sage_default_rank = combined.res_sage_rank.max()
    wide_default_rank = combined.res_wide_rank.max()
    
    new_rows = []
    for _, row in combined.iterrows():
        peptide, modifications = render_nice_peptide_strings(row['res_peptide'])
       
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
                   'Peptide_str': row['res_peptide'],  # This should be removed in the, uh, *final* final report
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
                   'Fragmentation Prediction Error': row['res_inc_frag_error'], # TODO: Which one is best?
                   'Isotope Profile Error': row['res_isotope_ratio_error'],
                   'Proteins': row['res_proteins'],
                   }
        
        new_rows.append(new_row)

    report = pd.DataFrame(new_rows, columns=['Spectrum_ID', 'Filename', 'Scan', 'Peptide', 'Modifications', 'Peptide_str', 'Charge', 'M/Z',
                                             'Calculated Mass', 'Experimental Mass', 'Retention Time', 'dMass', 
                                             'Search Engine', 'Mokapot q-value', 'Mokapot score', 'Mokapot PEP',
                                             'RT Prediction Error', 'Fragmentation Prediction Error', 'Isotope Profile Error',
                                             'Proteins'])
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
