import os, sys
from collections import defaultdict
import intervaltree as itree
import pandas as pd
import numpy as np
from mod_mass import identify_mod_by_mass

from deeplc import DeepLC



sys.path += '/'
from multiplierz.mzAPI import mzFile
from multiplierz.mass_biochem import fragment
from multiplierz.internalAlgorithms import ProxSeq, PPM_bounds, assign_multiprocess
from multiplierz.mzTools.isoDist import forPeptide
from multiplierz.spectral_process import peak_pick

from peptdeep.model.ms2 import pDeepModel 

try:
    PEPTDEEP_MODEL_PATH = '/data/models/peptdeep/generic/ms2.pth'
    assert(os.path.exists(PEPTDEEP_MODEL_PATH))
except AssertionError:
    PEPTDEEP_MODEL_PATH = '/generic/ms2.pth'
    assert(os.path.exists(PEPTDEEP_MODEL_PATH))



def annotate_with_feature_info(psms, featurefile, ppm_tol=20):
    features = pd.read_csv(featurefile, sep='\t')

    # Subtle point - rt_intervals has to be a defaultdict, even though the keys are floats, because
    # RT is based on MS1 time, and there's only so many MS1s, so collisiions become a big problem.
    rt_intervals = defaultdict(list)
    for ind, x in features.iterrows():
        width = x['rtEnd'] - x['rtStart']  # type: ignore
        rt_intervals[(x['rtStart'] - (width), x['rtEnd'])].append(ind)  # type: ignore
    mz_intervals = {PPM_bounds(ppm_tol, x['mz']):ind for ind, x in features.iterrows()}
    rt_tree = itree.IntervalTree([itree.Interval(*x) for x in rt_intervals.keys()])
    mz_tree = itree.IntervalTree([itree.Interval(*x) for x in mz_intervals.keys()])

    quant_rows = []
    for ind, row in psms.iterrows():
        mz = float(row['mz'])  # type: ignore
        rt = float(row['rt'])
        chg = int(row['charge'])  # type: ignore

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
            quant_rows.append({'_psm_id': ind, 'feature_id': None, 'feature_RT_start': None,
                               'feature_RT_end': None, 'feature_RT_apex': None,
                               'feature_total_int': None, 'feature_apex_int': None,
                               'feature_mz': None})
        else:
            best_match = match_data.iloc[(match_data['mz']-mz).argmin()]
            quant_rows.append({'_psm_id': ind,
                               'feature_id': best_match.name, 
                               'feature_RT_start': best_match['rtStart'],
                               'feature_RT_end': best_match['rtEnd'],
                               'feature_RT_apex': best_match['rtApex'],
                               'feature_total_int': best_match['intensitySum'],
                               'feature_apex_int': best_match['intensityApex'],
                               'feature_mz': best_match['mz'],
                               })


    quant = pd.DataFrame(quant_rows, columns = ['_psm_id', 'feature_id', 'feature_RT_start',
                                                'feature_RT_end', 'feature_RT_apex',
                                                'feature_total_int', 'feature_apex_int',
                                                'feature_mz'])
    assert(len(psms)==len(quant))
    new_psms = psms.merge(quant, left_index=True, right_on='_psm_id', how='left')#.drop('_psm_id', axis=1)
    assert(len(new_psms)==len(quant))
    assert(len(psms)==len(new_psms))
    #assert(not any(psms.feature_id.isnull()))
    #assert(not any(psms.feature_apex_int.isnull()))
    assert(not any(new_psms.index.isnull()))
    return new_psms


def new_annotate_with_feature_info(psms, featurefile, ppm_tol=20):
    all_features = pd.read_csv(featurefile, sep='\t')

    all_charge = all_features.charge.unique()

    quant_rows = []
    for charge in all_charge:
        features = all_features[all_features.charge==charge]

        # Subtle point - rt_intervals has to be a defaultdict, even though the keys are floats, because
        # RT is based on MS1 time, and there's only so many MS1s, so collisiions become a big problem.
        rt_intervals = defaultdict(list)
        for ind, x in features.iterrows():
            width = x['rtEnd'] - x['rtStart']  # type: ignore
            rt_intervals[(x['rtStart'] - (width), x['rtEnd'])].append(ind)  # type: ignore
        mz_intervals = {PPM_bounds(ppm_tol, x['mz']):ind for ind, x in features.iterrows()}
        rt_tree = itree.IntervalTree([itree.Interval(*x) for x in rt_intervals.keys()])
        mz_tree = itree.IntervalTree([itree.Interval(*x) for x in mz_intervals.keys()])

        for ind, row in psms[psms.charge==charge].iterrows():
            mz = float(row['mz'])  # type: ignore
            rt = float(row['rt'])
            chg = int(row['charge'])  # type: ignore

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
        
                    match_data = features.loc[list(matches)]
                    # match_data = match_data[match_data.charge.apply(int)==int(chg)]
                    assert(all(match_data.charge==int(chg)))
                    if len(match_data):
                        break
                if len(match_data):
                    break

            if len(match_data) == 0:
                quant_rows.append({'_psm_id': ind, 'feature_id': None, 'feature_RT_start': None,
                                   'feature_RT_end': None, 'feature_RT_apex': None,
                                   'feature_total_int': None, 'feature_apex_int': None,
                                   'feature_mz': None})
            else:
                best_match = match_data.iloc[(match_data['mz']-mz).argmin()]
                quant_rows.append({'_psm_id': ind,
                                   'feature_id': best_match.name, 
                                   'feature_RT_start': best_match['rtStart'],
                                   'feature_RT_end': best_match['rtEnd'],
                                   'feature_RT_apex': best_match['rtApex'],
                                   'feature_total_int': best_match['intensitySum'],
                                   'feature_apex_int': best_match['intensityApex'],
                                   'feature_mz': best_match['mz'],
                                   })


    quant = pd.DataFrame(quant_rows, columns = ['_psm_id', 'feature_id', 'feature_RT_start',
                                                'feature_RT_end', 'feature_RT_apex',
                                                'feature_total_int', 'feature_apex_int',
                                                'feature_mz'])
    assert(len(psms)==len(quant))
    new_psms = psms.merge(quant, left_index=True, right_on='_psm_id', how='left').drop('_psm_id', axis=1)
    assert(len(new_psms)==len(quant))
    assert(len(psms)==len(new_psms))
    #assert(not any(psms.feature_id.isnull()))
    #assert(not any(psms.feature_apex_int.isnull()))
    assert(not any(new_psms.index.isnull()))
    return new_psms


# This is extremely ad hoc; will eventually be replaced when parameter file autogeneration is built out (TODO)
MASS_TO_MOD = {
        229: 'TMT6plex', # Or 10, 11, etc.
        304: 'TMTpro',
        42: 'Acetyl',
        15: 'Oxidation',
        57: 'Carbamidomethyl',
        79: 'Phospho',
        }

# TODO it seems like deepLC silently throws away 0-indexed mods????? WTF?
def get_rt_pred_deltas(psms, calibration_size = 9000):
        
    def render_deeplc_mod(mod):
        if mod=='-' or pd.isnull(mod):
            return ''
        else:
            # agg = []
            # for thing in mod.split(','):
            #     a, c = thing.split('_')
            #     agg.append(a+'|'+c)
            loc_mods = [x.split('_') for x in mod.split(',')]

            # DeepLC's call into psm_utils to format the peptides will barf specifically on 
            # multiple N-term mods, so this just removes one. To be slightly less hacky, 
            # remove the smaller one, which is likely to make less (?) of a difference
            # to the elution.
            if len([x for x in loc_mods if int(x[0])==0]) > 1:
                extra_zero_mods = sorted([x for x in loc_mods if int(x[0])==0], key=lambda x: float(x[1]))[:-1]
                loc_mods = [x for x in loc_mods if x not in extra_zero_mods]
            loc_mods = [(loc if int(loc) else '1', mod) for loc, mod in loc_mods]

            loc_mods = [(loc, MASS_TO_MOD[int(float(mod))]) for loc, mod in loc_mods]

            return '|'.join(['%s|%s' % (a, c) for a, c in loc_mods])

    psms['deeplc_mods'] = psms.rec_mods.apply(render_deeplc_mod)
    
    # conf_hit_count = min([len(psms[psms['label']>0])/2, calibration_size])
    
    ## This is probably the right way to do it, roughly, but I was sorting in the wrong direction?
    # conf_hit_txt = psms[psms.proteins.apply(lambda x: 'rev_' not in str(x))].sort_values('e-value').iloc[:conf_hit_count]
    ## Well, can at least insert a shuffle, which is probably also the right way to do it.
    # sample(frac=1)
    conf_psms = psms[(psms.label > 0) & (~psms.feature_RT_apex.isnull()) & (psms.in_both > 0)]
    calibration_size = int(min([len(conf_psms)/2, calibration_size]))
    conf_hit_txt = conf_psms.sample(n=calibration_size)

    deeplc_calib = pd.DataFrame({'seq': conf_hit_txt.peptide, 'modifications': conf_hit_txt.deeplc_mods, 'tr': conf_hit_txt.feature_RT_apex})
    deeplc_target = pd.DataFrame({'seq': psms.peptide,
                                 'modifications': psms.deeplc_mods})

    dlc = DeepLC()  # Could add some arguments here, I guess?
    dlc.calibrate_preds(seq_df=deeplc_calib)

    deeplc_preds = dlc.make_preds(seq_df=deeplc_target)

    psms['rt_pred'] = deeplc_preds
    # psms['rt_error'] = psms['feature_RT_apex'] - psms['rt_pred']
    # psms[psms.feature_RT_apex.isnull()]['rt_error'] = psms[psms.feature_RT_apex.isnull()]['rt'] - psms[psms.feature_RT_apex.isnull()]['rt_pred']
    psms['rt_error'] = psms.apply(axis=1, func=lambda x: x['rt'] - x['rt_pred'] if pd.isnull(x['feature_RT_apex']) else x['feature_RT_apex'] - x['rt_pred'])
    psms['abs_rt_error'] = psms.rt_error.abs()

    assert(not any(psms['abs_rt_error'].isna()))

    return psms[['rt', 'rt_error', 'abs_rt_error']]


def calculate_fragment_error(data, pred_lookup, txt_item, ppm_tol=200):
    seq = txt_item['peptide']
    nice_mods = txt_item['alphapept_mods']
    mods = txt_item['rec_mods']
    charge = txt_item['charge']
    scannum = txt_item['scannr']

    scan = ProxSeq(data.scan(scannum, centroid=True), lambda x: x[0])

    pred_info = pred_lookup[seq, mods if mods!='' else np.nan, charge]

    # Redundant! But multiplierz only recognizes a very limited set of mod names.
    frag_mods = []
    for mod in mods.split(','):
        loc, mass = mod.split('_')
        frag_mods.append(f'{seq[int(loc)-1]}{loc}: {float(mass)}')

    frags = fragment(seq.replace('U', 'S'), frag_mods, charges=[1,2])

    obs_pred_pairs = []
    for _, row in pred_info.iterrows():
        index = int(row['cleave_ind']-1)

        for ion in ['b', 'y']:
            for charge_num, charge_plus in [('1', ''), ('2', '++')]:
                pred_int = row[f'{ion}_z{charge_num}']

                mz = frags[ion+charge_plus][index][1]
                zone = scan.returnRange(*PPM_bounds(ppm_tol, mz))
                try:
                    obs_int = min(zone, key = lambda x: abs(x[1]-pred_int))[1]
                except ValueError:
                    obs_int = 0

                # diff = obs_int / pred_int
                # obses.append((obs_int, len(zone)))
                # diffs.append(diff)
                obs_pred_pairs.append((obs_int, pred_int))

    match_pairs = [x for x in obs_pred_pairs if x[0] and x[1]]
    if not match_pairs:
        return len([x for x in obs_pred_pairs if x[0] or x[1]]), 0
    else:
        scale_factor = sum([x[0] for x in match_pairs]) / sum([x[1] for x in match_pairs])

        diffs = [x[0] / (x[1]*scale_factor) for x in match_pairs]
        
        return abs(np.mean(diffs)-1), len(diffs)


# TODO: specific NCE and instrument setting?
def get_frag_pred_deltas(psms, data, nce, instrument, model_path = PEPTDEEP_MODEL_PATH,
                         prediction_batch_size = 100000):
    def render_alphapept_mod(mod, pep):
        if mod=='-' or pd.isnull(mod):
            return ''
        else:
            agg = []
            for thing in mod.split(','):
                loc, mass = thing.split('_')  # Letter feels like it should be AA but it isn't?
                if thing.endswith('_n'):
                    aa = 'Any N-term'
                else:
                    aa = pep[int(loc)-1]
                modname = identify_mod_by_mass(float(mass))
                agg.append(('%s%s: %s' % (aa, loc, modname)))
            return '; '.join(agg)


    def _deeppept_format_table(pep_list, charges=[2,3,4]):
        mods = []
        mod_sites = []

        for item in pep_list:
            peptide = item[0]
            modstr = item[2]
            if pd.isnull(modstr):
                mods.append('')
                mod_sites.append('')
            else:
                submods = []
                submod_sites = []
                for mod in modstr.split(','):
                    loc, mass = mod.split('_')
                    if loc=='0':
                        site = 'Any N-term'
                    else:
                        site = peptide[int(loc)-1]
                    modname = identify_mod_by_mass(float(mass))
                    if modname == 'Oxidation':
                        assert(site=='M')
                    if modname == 'Acetyl':
                        assert(site=='Any N-term')
                    submods.append(modname + '@' + site)
                    submod_sites.append(loc)
                mods.append(';'.join(submods))
                mod_sites.append(';'.join(submod_sites))

        base_tab = pd.DataFrame({
            'pep_name': [x[1] for x in pep_list],
            'sequence': [x[0] for x in pep_list],
            'mods': mods,
            'mod_sites': mod_sites,
            'fmt_mods': [x[2] for x in pep_list],
            'nce': nce,
            'instrument': instrument
        })
        if charges is None:
            base_tab['charge'] = [x[3] for x in pep_list]
            final_tab = base_tab
        else:
            charged_tabs = []
            for charge in charges:
                subtab = base_tab.copy()
                subtab['charge'] = charge
                charged_tabs.append(subtab)
            final_tab = pd.concat(charged_tabs)

        return final_tab


    def _predict_for_specific(peptide_list, model):
        prediction_record = {}
        for p_i in range(0, len(peptide_list), prediction_batch_size):
            sub_list = peptide_list[p_i:p_i+prediction_batch_size]
            pep_tab = _deeppept_format_table(sub_list, charges=None)
            pred_tab = model.predict(pep_tab)

            for ind, row in pep_tab.iterrows():
                pep_pred = pred_tab[pred_tab.pep_ind==ind]
                prediction_record[row.sequence, row.fmt_mods, row.charge] = pep_pred

        return prediction_record

    fragmodel = pDeepModel(device='cpu')
    assert(os.path.exists(model_path))
    fragmodel.load(model_path) # type: ignore

    psms['alphapept_mods'] = psms.apply(axis=1,
                                       func=lambda x: render_alphapept_mod(x['rec_mods'], 
                                                                           x['peptide']))
    psms['name_placeholder'] = list(psms.index)

    # This could be un-rolled (re-rolled?) so that it only calculates the prediction for each peptide
    # once across all files
    peptide_list = list(psms[['peptide', 'name_placeholder', 'rec_mods', 'charge']]
                        .itertuples(index=False, name=None))
    preds = _predict_for_specific(peptide_list, fragmodel)
    frag_errors = psms.apply(axis=1, func=lambda x: calculate_fragment_error(data, preds, x, ppm_tol=20))
    psms['fragment_error'] = [x[0] for x in frag_errors]
    psms['matched_frags'] = [x[1] for x in frag_errors]


    return psms[['fragment_error', 'matched_frags']]


def get_isotopic_envelope_delta(psms, data):

    def _expected_ratio_for_peptide(peptide):
        ratios = forPeptide(peptide.replace('U', 'S'))
        return ratios[0]/ratios[1]


    ms2s_for_ms1 = {}
    agg = []
    for info in data.scan_info():
        _, mz, num, level, _ = info
        if level=='MS2':
            agg.append(info)
        if level=='MS1':
            ms2s_for_ms1[num] = agg
            agg = []
    
    precursor_ratios = []
    for ms1, ms2s in ms2s_for_ms1.items():
        if ms2s:
            ms1scan = data.scan(ms1, centroid=True)
            envelopes, _ = peak_pick(ms1scan)  # type: ignore
            for _, mz, num, _, _ in ms2s:
                exp_chg = data.extra_info(num)['Charge State']
                try:
                    obs_env = min(envelopes[exp_chg], key = lambda x: min(abs(x[0][0]-mz), abs(x[1][0]-mz)))
                    if (abs(obs_env[0][0]-mz) < 0.01 or abs(obs_env[1][0]-mz) < 0.01):
                        obs_ratio = obs_env[0][1]/obs_env[1][1]
                    else:
                        obs_ratio = np.nan
                except KeyError:
                    obs_ratio = np.nan
                precursor_ratios.append((num, obs_ratio))
    precursor_iso_by_scan = dict(precursor_ratios)

    psms['expected_isotope_ratio'] = psms.apply(axis=1,
                                               func=lambda x: _expected_ratio_for_peptide(x['peptide']))
    psms['observed_isotope_ratio'] = psms.apply(axis=1, func=lambda x: precursor_iso_by_scan.get(x['scannr'], 1))
    psms['isotope_ratio_error'] = abs(np.log(psms['observed_isotope_ratio']) -
                                      np.log(psms['expected_isotope_ratio']))

    # Replace nan with median
    psms['isotope_ratio_error'] = psms['isotope_ratio_error'].fillna(psms['isotope_ratio_error'].median())

    return psms['isotope_ratio_error']

import time
def annotate_pin_file(psmfile, rawfile, featurefile, outputfile, coll_e = 32, instrument_name = 'Elite'):
    testout = open(os.path.join(os.path.dirname(outputfile), 'testout.txt'), 'w')
    testout.write(f"Starting annotation, time: {time.time()}\n")
    data = mzFile(rawfile) 
    psms = pd.read_csv(psmfile)

    # print("TEST MODE")
    # psms = psms.iloc[:100]
    
    testout.write(f"Loaded data, time: {time.time()}\n")

     # TODO replace this with numbers from the TMT file?
    mz_rt_lookup = {x[2]:(x[0], x[1]) for x in data.scan_info()} # type: ignore
    # psms[['rt', 'mz']] = psms.scannr.apply(lambda x: mz_rt_lookup[x])
    psms['rt'] = psms.scannr.apply(lambda x: mz_rt_lookup[x][0])
    psms['mz'] = psms.scannr.apply(lambda x: mz_rt_lookup[x][1])

    psms = annotate_with_feature_info(psms, featurefile)

    # Basic imputation to make mokapot work better?  TODO check if taking median is better?
    psms['feature_apex_int'] = psms['feature_apex_int'].fillna(psms['feature_apex_int'].min())
    psms['feature_total_int'] = psms['feature_total_int'].fillna(psms['feature_total_int'].min())
    psms['feature_run_length'] = psms['feature_RT_end'] - psms['feature_RT_start']
    psms['feature_run_length'] = psms['feature_run_length'].fillna(psms['feature_run_length'].min())


    peptide_table = psms.sort_values('feature_apex_int', ascending=False).drop_duplicates(['peptide', 'rec_mods']).copy()

    # TODO multiprocessing/multithreading?
    # psms[['rt', 'rt_error', 'abs_rt_error']] = get_rt_pred_deltas(psms.copy())
    # psms[['fragment_error', 'matched_fragments']] = get_frag_pred_deltas(psms.copy(), data, nce=coll_e, instrument=instrument_name)
    # psms['isotope_ratio_error'] = get_isotopic_envelope_delta(psms.copy(), data)

    peptide_table[['rt', 'rt_error', 'abs_rt_error']] = get_rt_pred_deltas(peptide_table.copy())
    peptide_table[['fragment_error', 'matched_fragments']] = get_frag_pred_deltas(peptide_table.copy(), data, nce=coll_e, instrument=instrument_name)
    peptide_table['isotope_ratio_error'] = get_isotopic_envelope_delta(peptide_table.copy(), data)

    psms = psms.merge(peptide_table[['peptide', 'rec_mods', 'rt_error', 'abs_rt_error', 'fragment_error', 'matched_fragments', 'isotope_ratio_error']],
                      on=['peptide', 'rec_mods'], how='left')

    psms.drop('rec_mods', axis = 1, inplace = True)

    # Rrebuild specid since mokapot will expect it.
    psms['specid'] = psms.apply(axis=1, func=lambda x: '%s_%s' % (x['filename'], x['scannr']))

    psms.to_csv(outputfile, sep='\t', index=False)
    data.close()

    print("Annotated PSMs written to %s" % outputfile)
    return outputfile

def test_feature_thing(psmfile, rawfile, featurefile, etc):
    data = mzFile(rawfile) 
    psms = pd.read_csv(psmfile)

    # print("TEST MODE")
    # psms = psms.iloc[:100]
    

     # TODO replace this with numbers from the TMT file?
    mz_rt_lookup = {x[2]:(x[0], x[1]) for x in data.scan_info()} # type: ignore
    # psms[['rt', 'mz']] = psms.scannr.apply(lambda x: mz_rt_lookup[x])
    psms['rt'] = psms.scannr.apply(lambda x: mz_rt_lookup[x][0])
    psms['mz'] = psms.scannr.apply(lambda x: mz_rt_lookup[x][1])

    start = (time.time())
    bar = new_annotate_with_feature_info(psms.copy(), featurefile)
    print("Elapsed time: ", time.time()-start)
    start = (time.time())
    foo = annotate_with_feature_info(psms.copy(), featurefile)
    print("Elapsed time: ", time.time()-start)

    print(bar.equals(foo))
    foo.to_csv(etc + '.old.csv', sep='\t', index=True)
    bar.to_csv(etc + '.new.csv', sep='\t', index=True)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('pin')
    parser.add_argument('raw')
    parser.add_argument('features')
    parser.add_argument('output')
    parser.add_argument('--coll_e', default = 32)
    parser.add_argument('--instrument_name', default = 'Elite')
    args = parser.parse_args()


    # import cProfile
    # cProfile.run('annotate_pin_file(args.pin, args.raw, args.features, args.output)#, args.coll_e, args.instrument_name)', '/data/annotate_psms_via_models.prof')
    # annotate_pin_file(args.pin, args.raw, args.features, args.output)
    test_feature_thing(args.pin, args.raw, args.features, args.output)
