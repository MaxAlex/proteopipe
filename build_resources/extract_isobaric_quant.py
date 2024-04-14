import os
import numpy as np
import pandas as pd
from pyteomics import mzml, auxiliary
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


CHANNEL_CONFIGS = {4: list(zip(['114', '115', '116', '117'], 
                               [114.11123, 115.10826, 116.11162, 117.11497])),
                   6: list(zip(['126', '127', '128', '129', '130', '131'],
                         [126.127726, 127.131081, 128.134436,
                          129.137790, 130.141145, 131.138180])),
                   8: list(zip(['113', '114', '115', '116', '117', '118', '119', '121'],
                         [113.11, 114.11, 115.11, 116.11, 117.12, 118.12, 119.12, 121.12])),
                   10: list(zip(['126', '127N', '127C', '128N', '128C', 
                          '129N', '129C', '130N', '130C', '131'],
                         [126.127726, 127.124761, 127.131081, 128.128116, 128.134436,
                          129.131471, 129.137790, 130.134825, 130.141145, 131.138180])),
                   11: list(zip(['126', '127N', '127C', '128N', '128C',
                          '129N', '129C', '130N', '130C', '131N', '131C'],
                         [126.127726, 127.124761, 127.131081, 128.128116, 128.134436,
                          129.131471, 129.137790, 130.134825, 130.141145, 131.138180, 131.144499])),

                    # # PULLED OUT OF COPILOT'S ASS, DO NOT USE
                    # 16: list(zip(['126', '127N', '127C', '128N', '128C',
                    #         '129N', '129C', '130N', '130C', '131N', '131C',
                    #         '132N', '132C', '133N', '133C', '134'],
                    #          [126.127726, 127.124761, 127.131081, 128.128116, 128.134436,
                    #         129.131471, 129.137790, 130.134825, 130.141145, 131.138180, 131.144499,
                    #         132.141534, 132.147853, 133.144888, 133.151207, 134.148242]))

                   # TODO What other isobaric tag setups are there these days???
}


CALIBRANT_MASS = 445.120025
CALIBRANT_TOL = 0.003


def read_mzml_quant_info(mzml_file, tmt_channels, tol = 0.01, inclusive_mode = True):

    with mzml.read(mzml_file) as reader:
        for spectrum in reader:
            if spectrum['ms level'] == 2:
                scan_id = spectrum['id']
                scan_id = int(scan_id.split('scan=')[-1])
                scan_index = spectrum['index']
                scan = np.stack((spectrum['m/z array'], spectrum['intensity array']), axis=1)

                quant = []
                for chan_name, chan_mz in tmt_channels:
                    match_peaks = scan[np.abs(scan[:, 0] - chan_mz) < tol]
                    if inclusive_mode:
                        chan_int = np.sum(match_peaks[:, 1])
                        if chan_int == 0:
                            chan_int = np.nan
                        try:
                            chan_match_mz_err = match_peaks[:, 0].max() - chan_mz
                        except ValueError:
                            chan_match_mz_err = np.nan
                    else:
                        try:
                            chan_int = np.max(match_peaks[:, 1])
                            chan_match_mz_err = match_peaks[np.argmax(match_peaks[:, 1]), 0] - chan_mz
                        except ValueError:
                            chan_int = np.nan 
                            chan_match_mz_err = np.nan 
                    quant.append((chan_name, chan_int, chan_match_mz_err))

                assert(spectrum['precursorList']['count'] == 1)
                assert(spectrum['precursorList']['precursor'][0]['selectedIonList']['count'] == 1)
                precursor = spectrum['precursorList']['precursor'][0]

                source_ms1 = precursor['spectrumRef']
                source_mz = precursor['isolationWindow']['isolation window target m/z']
                source_range = precursor['isolationWindow']['isolation window lower offset'], precursor['isolationWindow']['isolation window upper offset']
                source_charge = precursor['selectedIonList']['selectedIon'][0]['charge state']
                source_int = precursor['selectedIonList']['selectedIon'][0]['peak intensity']

                cal_region = scan[np.abs(scan[:, 0] - CALIBRANT_MASS) < CALIBRANT_TOL]
                if len(cal_region) == 0:
                    cal_peak_int = 0
                    cal_peak_err = 0.005
                else: 
                    cal_peak_int = cal_region[:, 1].max()
                    cal_peak_err = cal_region[np.argmax(cal_region[:, 1]), 0] - CALIBRANT_MASS

                yield scan_id, scan_index, source_ms1, source_mz, source_range, source_charge, source_int, quant, cal_peak_int, cal_peak_err


def normalize_tmt_channels(table, tmt_channels):
    channel_names = [x[0] for x in tmt_channels]

    print("Base median values:")
    print(table[channel_names].median().div(table[channel_names].median().median()))

    channel_values = table[channel_names].dropna()
    # Taking the top half of the rows by total intensity, to avoid noise
    channel_values = channel_values.loc[channel_values.sum(axis=1).sort_values(ascending=False).index[:len(channel_values)//2]]

    scaled_channels = channel_values.div(channel_values.sum(axis=1), axis=0)
    medians = scaled_channels.median()
    
    print("Median values:")
    print(medians)
    print("Pipetting error ratios:")
    pipette_errors = medians.div(medians.median())
    print(medians.div(medians.median()))
   
    table[channel_names] = table[channel_names].div(pipette_errors, axis=1)

    print("Corrected median values:")
    print(table[channel_names].median().div(table[channel_names].median().median()))

    return table


def check_for_skew(table, tmt_channels, tol = 0.01):
    channel_accs = [x[0] + '_mz' for x in tmt_channels] 
    mz_errors = table[channel_accs].abs()
    error_within_bounds = (mz_errors.mean() + (2*mz_errors.std())) < tol

    if not all(error_within_bounds):
        print("Warning: Some channels have skew in their peak matching that will effect quantification accuracy in >1% of cases.")
        print(error_within_bounds)
        print(mz_errors.mean())
        print(mz_errors.std())
   

def classify_valid_quants(table, tmt_channels):
    channel_names = [x[0] for x in tmt_channels]
    known_good = table[table[channel_names].notnull().all(axis=1)].sort_values(by='source_int', ascending=False)[:len(table)//10]

    table['missing_channels'] = table[table[channel_names].isnull().sum(axis=1)]
    table['max_mz_error'] = table[[x[0] + '_mz' for x in tmt_channels]].abs().max(axis=1)
    table['max_int'] = table[channel_names].max(axis=1)
    table['median_int'] = table[channel_names].median(axis=1)
    table['max_int_ratio'] = table[channel_names].apply(lambda x: x.dropna().max()/x.sum(), axis=1)
    
    qual_factors = table[['source_int', 'missing_channels', 'max_mz_error', 'max_int', 'median_int', 'max_int_ratio']]
    scaled_factors = StandardScaler().fit_transform(qual_factors)
    pca = PCA(n_components=2)
    pca.fit(scaled_factors)




    # TODO ???

    

def extract_isobaric_from_mzml(mzml_file, output_file, tmt_channels, tol = 0.01, inclusive_mode = True):
    filestem = os.path.basename(mzml_file).split('.')[0]
    report = []
    for scan_info in read_mzml_quant_info(mzml_file, tmt_channels, tol, inclusive_mode):
        scan_id, scan_index, source_ms1, source_mz, source_range, source_charge, source_int, quant, cal_int, cal_err = scan_info
        row = {'filename': filestem, 'scan_id': scan_id, 'scan_index': scan_index, 'source_ms1': source_ms1, 'source_mz': source_mz, 
               'source_range': source_range, 'source_charge': source_charge, 'source_int': source_int,
               'calibrant_int': cal_int, 'calibrant_err': cal_err}
        row.update({x[0]: x[1] for x in quant})
        row.update({'%s_mz' % x[0]: x[2] for x in quant})
        report.append(row)

    column_order = (['filename', 'scan_id', 'scan_index', 'source_ms1', 'source_mz', 'source_range', 'source_charge', 'source_int', 
                     'calibrant_int', 'calibrant_err'] + [x[0] for x in tmt_channels] + ['%s_mz' % x[0] for x in tmt_channels])
    table = pd.DataFrame(report, columns=column_order)
   
    table.to_csv(output_file + 'before_normalize.csv', index=False)

    check_for_skew(table, tmt_channels, tol)
    table = normalize_tmt_channels(table, tmt_channels)

    table.to_csv(output_file, index=False)

    return output_file


def integrate_quant(psm_file, isobaric_datafiles, output_file):
    psms = pd.read_csv(psm_file)
    isobaric = pd.concat([pd.read_csv(x) for x in isobaric_datafiles])

    # Merge on both scan number and filename
    combined = pd.merge(psms, isobaric, left_on=['Scan', 'Filename'], right_on=['scan_id', 'filename'])
    assert((combined['M/Z'] - combined['source_mz'] < 0.01).all())  # Ensuring that we're matching the right scans

    combined = combined.drop(columns=['scan_id', 'scan_index', 'source_ms1', 'source_mz', 'source_range',
                                      'source_charge', 'source_int', 'calibrant_int', 'calibrant_err'])
    combined = combined.drop([x for x in combined.columns if x.endswith('_mz')], axis=1)

    combined.to_csv(output_file, index=False) 

    return combined


if __name__ == '__main__':
    # table = pd.read_csv('test_output.csvbefore_normalize.csv')
    # tmt_channels = CHANNEL_CONFIGS[11]
    # normalize_tmt_channels(table, tmt_channels)
    import argparse
    parser = argparse.ArgumentParser(description='Generate a report of TMT quantification from an mzML file')
    parser.add_argument('mode', help='Run mode: extract, integrate')
    parser.add_argument('--mzml_file', help='The mzML file to read')
    parser.add_argument('--output_file', help='The output file to write')
    parser.add_argument('--psm_file', help='The PSM file to integrate isobaric quant to')
    parser.add_argument('--isobaric_datafiles', nargs='+', help='The isobaric quant files to integrate into the PSM file')
    parser.add_argument('--channel_mode', help='TMT/iTRAQ channel configuration to use')
    parser.add_argument('--tol', type=float, default=0.01, help='The tolerance to match peaks')
    parser.add_argument('--inclusive', action='store_true', help='Use inclusive mode for peak matching')
    args = parser.parse_args()

    if args.mode == 'extract':
        channels = CHANNEL_CONFIGS[int(args.channel_mode)]
        extract_isobaric_from_mzml(args.mzml_file, args.output_file, channels, args.tol, args.inclusive)
    elif args.mode == 'integrate':
        integrate_quant(args.psm_file, args.isobaric_datafiles, args.output_file)
    else:
        raise ValueError("Invalid mode %s" % args.mode)

