import os
import pandas as pd
import intervaltree as itree
from collections import defaultdict


# Cases exist where the machine is slightly more aggressive than dinosaur, but this is usually rare.
WIDENING_FACTOR = 0.01  


def PPM_bounds(ppm, mass):
    val = (mass/1000000.0)*ppm
    return mass - val, mass + val


def annotate_psms_with_features(psmfile, featurefiles, ppm_tol = 20,
                                mz_col = 'M/Z', rt_col = 'Retention Time', chg_col = 'Charge'): 
    psms = pd.read_csv(psmfile)

    feature_info = {}
    for featurefile in featurefiles:
        features = pd.read_csv(featurefile, sep='\t')

        # Subtle point - rt_intervals has to be a defaultdict, even though the keys are floats, because
        # RT is based on MS1 time, and there's only so many MS1s, so collisiions become a big problem.
        rt_intervals = defaultdict(list)
        for ind, x in features.iterrows():
            width = x['rtEnd'] - x['rtStart']  # type: ignore
            rt_intervals[(x['rtStart'] - (width*WIDENING_FACTOR), x['rtEnd'])].append(ind)  # type: ignore
        mz_intervals = {PPM_bounds(ppm_tol, x['mz']):ind for ind, x in features.iterrows()}
        rt_tree = itree.IntervalTree([itree.Interval(*x) for x in rt_intervals.keys()])
        mz_tree = itree.IntervalTree([itree.Interval(*x) for x in mz_intervals.keys()])

        filestem = os.path.basename(featurefile).split('.')[0]
        feature_info[filestem] = (features, rt_intervals, mz_intervals, mz_tree, rt_tree)

    quant_rows = []
    for ind, row in psms.iterrows():
        mz = float(row[mz_col])  # type: ignore
        rt = row[rt_col]
        chg = int(row[chg_col])  # type: ignore

        features, rt_intervals, mz_intervals, mz_tree, rt_tree = feature_info[row['Filename']]

        rt_matches = set(sum((rt_intervals[tuple(x[:2])] for x in rt_tree[rt]), [])) # type: ignore
        # The feature will always (should always) have the correct monoisotopic MZ listed, but the
        # scan MZ may have been taken from a different isotopic peak. So we'll check for some reasonable
        # number of peaks and take the first one that matches. 
        match_data : pd.DataFrame = None #type: ignore
        iso = None
        for abs_iso in range(0, 7):
            for iso in [-abs_iso, abs_iso]:
                iso_mz = mz - (iso * 1.00335)/chg
                mz_matches = {mz_intervals[tuple(x[:2])] for x in mz_tree[iso_mz]}
                matches = rt_matches.intersection(mz_matches)
    
                match_data = features.iloc[list(matches)]
                match_data = match_data[match_data.charge.apply(int)==int(chg)]
                if len(match_data):
                    break
            if len(match_data):
                break


        if len(match_data) == 0:
            quant_rows.append({'PSM ID': None, 'Feature ID': None, 'Feature RT Start': None,
                               'Feature RT End': None, 'Feature Apex': None,
                               'Feature Intensity': None, 'Feature Apex Int': None,
                               'Feature Apex Int': None, 'Feature MZ': None
                               })
        else:
            best_match = match_data.iloc[(match_data['mz']-mz).argmin()]
            quant_rows.append({'PSM ID': ind,
                               'Feature ID': best_match.name, 
                               'Feature RT Start': best_match['rtStart'],
                               'Feature RT End': best_match['rtEnd'],
                               'Feature RT Apex': best_match['rtApex'],
                               'Feature Total Int': best_match['intensitySum'],
                               'Feature Apex Int': best_match['intensityApex'],
                               'Feature MZ': best_match['mz'],
                               'Feature Isotope Match': iso,
                               'Feature Interval': (rt-best_match['rtStart'])/(best_match['rtEnd']-best_match['rtStart'])
                               })

    quant = pd.DataFrame(quant_rows, columns = ['PSM ID', 'Feature ID', 'Feature RT Start', 'Feature RT End', 'Feature RT Apex', 'Feature Total Int', 
                                                'Feature Apex Int', 'Feature MZ', 'Feature Interval', 'Feature Isotope Match'])
    psms_wfeatures = psms.merge(quant, left_index=True, right_on='PSM ID', how='left').drop('PSM ID', axis=1)
    assert(len(psms_wfeatures) == len(psms))
    psms_wfeatures.to_csv(psmfile + '.w_features.csv', index=False)

        
        




if __name__ == '__main__':
   # annotate_psms_with_features('/data/test_data/draftprot_testcase/Fetal_Ovary_bRP_Elite_25_f08.report.csv',
   #                             '/data/test_data/draftprot_testcase/Fetal_Ovary_bRP_Elite_25_f08.dinosaur_features.tsv',
   #                             ppm_tol = 200)
    import argparse
    parser = argparse.ArgumentParser() 
    parser.add_argument('-p', '--psm', required=True)
    parser.add_argument('-d', '--data', nargs='+', required=True)
    args = parser.parse_args()
    annotate_psms_with_features(args.psm, args.data)
