import os
from collections import defaultdict
import intervaltree as itree
from multiplierz.mzAPI import mzFile
from multiplierz.mass_biochem import fragment
from multiplierz.internalAlgorithms import ProxSeq, PPM_bounds, assign_multiprocess
from multiplierz.mzTools.isoDist import forPeptide
from multiplierz.spectral_process import peak_pick
from multiplierz.mgf import extract
import pandas as pd
import numpy as np
import re
from mod_mass import identify_mod_by_mass

from deeplc import DeepLC

from peptdeep.model.ms2 import pDeepModel 

try:
    PEPTDEEP_MODEL_PATH = '/data/models/peptdeep/generic/ms2.pth'
    assert(os.path.exists(PEPTDEEP_MODEL_PATH))
except AssertionError:
    PEPTDEEP_MODEL_PATH = '/generic/ms2.pth'
    assert(os.path.exists(PEPTDEEP_MODEL_PATH))



def PPM_bounds(ppm, mass):
    val = (mass/1000000.0)*ppm
    return mass - val, mass + val

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
    new_psms = psms.merge(quant, left_index=True, right_on='_psm_id', how='left').drop('_psm_id', axis=1)
    assert(len(new_psms)==len(quant))
    assert(len(psms)==len(new_psms))
    #assert(not any(psms.feature_id.isnull()))
    #assert(not any(psms.feature_apex_int.isnull()))
    assert(not any(new_psms.index.isnull()))
    return new_psms


# This is extremely ad hoc; will eventually be replaced when parameter file autogeneration is built out (TODO)
MASS_TO_MOD = {
        229: 'TMT',
        304: 'TMTpro',
        42: 'Acetyl',
        15: 'Oxidation',
        57: 'Carbamidomethyl',
        79: 'Phospho',
        }


def get_rt_pred_deltas(psms):
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

            loc_mods = [(loc, MASS_TO_MOD[int(float(mod))]) for loc, mod in loc_mods]

            return '|'.join(loc_mods) # type: ignore
        
    psms['deeplc_mods'] = psms.rec_mods.apply(render_deeplc_mod)
    
    conf_hit_count = min([len(psms[psms.proteins.apply(lambda x: 'rev_' not in str(x))])/2, 2000])
    
    ## This is probably the right way to do it, roughly, but I was sorting in the wrong direction?
    # conf_hit_txt = psms[psms.proteins.apply(lambda x: 'rev_' not in str(x))].sort_values('e-value').iloc[:conf_hit_count]
    ## Well, can at least insert a shuffle, which is probably also the right way to do it.
    conf_hit_txt = psms[psms.proteins.apply(lambda x: 'rev_' not in str(x)) & (~psms.feature_RT_apex.isnull())].sample(frac=1).iloc[:conf_hit_count]

    deeplc_calib = pd.DataFrame({'seq': conf_hit_txt.peptide,
                                 'modifications': conf_hit_txt.deeplc_mods,
                                 'tr': conf_hit_txt.feature_RT_apex})
    deeplc_target = pd.DataFrame({'seq': psms.peptide,
                                 'modifications': psms.deeplc_mods})

    dlc = DeepLC()  # Could add some arguments here, I guess?
    dlc.calibrate_preds(seq_df=deeplc_calib)

    deeplc_preds = dlc.make_preds(seq_df=deeplc_target)

    psms['rt_pred'] = deeplc_preds
    psms['rt_error'] = psms['feature_RT_apex'] - psms['rt_pred']
    psms[psms.feature_RT_apex.isnull()]['rt_error'] = psms['rt'] - psms['rt_pred']
    psms['abs_rt_error'] = psms.rt_error.abs()

    assert(psms['abs_rt_error'].any())

    psms = psms.drop('deeplc_mods', axis=1)

    return psms 


def calculate_fragment_error(data, pred_lookup, ppm_tol=200):

    def get_fragmentation_error(txt_item):
        seq = txt_item['peptide']
        nice_mods = txt_item['alphapept_mods']
        mods = txt_item['rec_mods']
        charge = txt_item['charge']
        scannum = txt_item['scannr']

        scan = ProxSeq(data.scan(scannum, centroid=True), lambda x: x[0])

        pred_info = pred_lookup[seq, mods if mods!='' else np.nan, charge]
        frags = fragment(seq, nice_mods, charges=[1,2])
        value_pairs = []
        for _, row in pred_info.iterrows():
            b_i = row['cleave_ind']
            y_i = len(seq)-(b_i)

            b_obs = sum([x[1] for x in scan.returnRange(*PPM_bounds(ppm_tol, frags['b'][b_i]))])
            y_obs = sum([x[1] for x in scan.returnRange(*PPM_bounds(ppm_tol, frags['y'][y_i]))])
            bpp_obs = sum([x[1] for x in scan.returnRange(*PPM_bounds(ppm_tol, frags['b++'][b_i]))])
            ypp_obs = sum([x[1] for x in scan.returnRange(*PPM_bounds(ppm_tol, frags['y++'][y_i]))])
        
            # TODO mod loss values?
            value_pairs += [(b_obs, row['b_z1'], 'b%d' % b_i),
                            (y_obs, row['y_z1'], 'y%d' % y_i),
                            (bpp_obs, row['b_z2'], 'b++%d' % b_i),
                            (ypp_obs, row['y_z2'], 'y++%d' % y_i)]
            
        # validvals = [x for x in value_pairs if x[0] and x[1]]
        # median_scale = np.median([x[1]/x[2] for x in validvals])
        valid_obs_sum = sum([x[0] for x in value_pairs if x[1]])
        if not valid_obs_sum:
            return 1.5, 1.5
        else:
            scaled_pairs = [(x[0]/valid_obs_sum, x[1], x[2]) for x in value_pairs]

            # Square the errors?
            exc_err = sum([abs(x[0]-x[1]) for x in scaled_pairs if x[0] and x[1]]) / len([x[0] for x in value_pairs if x[0] and x[1]])
            inc_err = sum([abs(x[0]-x[1]) for x in scaled_pairs]) / len(scaled_pairs)

            return exc_err, inc_err
    
    return get_fragmentation_error





def pep_table(pep_list, nce, instrument, charges=[2,3,4]):
    mods = []
    mod_sites = []

    # for item in pep_list:
    #     modstr = item[2]
    #     mod_parse = [re.match(r'^([A-Za-z-]+)(\d*): ([A-Za-z0-9-]+)$', x) for x in modstr.split('; ')]
    #     mod_parse = [x.groups() if x is not None else (None, None, None) for x in mod_parse]

    #     submods = []
    #     submod_sites = []
    #     for aa_loc, site, modname in mod_parse:
    #         if (modname == 'Acetyl' or modname == 'TMT6plex') and site=='0':
    #             aa_loc = 'Any N-term'

    #         if modname is not None:
    #             newmod = modname + '@' + aa_loc
    #             assert('Acetyl@N' not in newmod)
    #             submods.append(newmod)
    #             submod_sites.append(site if site != '' else 0)
    #     mods.append(';'.join(submods))
    #     mod_sites.append(';'.join([str(x) for x in submod_sites]))

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


def predict_for_specific(peptide_list, model, nce, instrument, batch_size=10000):
    prediction_record = {}
    for p_i in range(0, len(peptide_list), batch_size):
        sub_list = peptide_list[p_i:p_i+batch_size]
        pep_tab = pep_table(sub_list, nce=nce, instrument=instrument, charges=None)
        pred_tab = model.predict(pep_tab)

        for ind, row in pep_tab.iterrows():
            pep_pred = pred_tab[pred_tab.pep_ind==ind]
            prediction_record[row.sequence, row.fmt_mods, row.charge] = pep_pred

    return prediction_record





# TODO: specific NCE and instrument setting?
def get_frag_pred_deltas(psms, data, nce, instrument,
                         model_path = PEPTDEEP_MODEL_PATH):
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

    fragmodel = pDeepModel(device='gpu')
    assert(os.path.exists(model_path))
    fragmodel.load(model_path)

    psms['alphapept_mods'] = psms.apply(axis=1,
                                       func=lambda x: render_alphapept_mod(x['rec_mods'], 
                                                                           x['peptide']))
    psms['name_placeholder'] = list(psms.index)

    # sub_psms_agg = []
    # for rawfile in rawfiles:
    #     file_psms = psms[psms['source_file'].apply(lambda x: os.path.basename(rawfile.split('.')[0])==x)]
    #     assert(len(file_psms))

    # This could be un-rolled (re-rolled?) so that it only calculates the prediction for each peptide
    # once across all files
    peptide_list = list(psms[['peptide', 'name_placeholder', 'rec_mods', 'charge']]
                        .itertuples(index=False, name=None))
    preds = predict_for_specific(peptide_list, fragmodel, nce=nce, instrument=instrument,
                                batch_size=1000)
    frag_error_calc = calculate_fragment_error(data, preds, ppm_tol=20)
    frag_errors = psms.apply(axis=1, func=frag_error_calc)
    psms['exc_frag_error'] = [x[0] for x in frag_errors]
    psms['inc_frag_error'] = [x[1] for x in frag_errors]

    return psms 


def get_precursor_envelope_ratios(data):
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
            envelopes, etc = peak_pick(ms1scan)  # type: ignore
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
    return precursor_ratios


def get_isotopic_envelope_delta(psms, data):
    def expected_ratio_for_peptide(peptide):
        ratios = forPeptide(peptide)
        return ratios[0]/ratios[1]

    # percursor_iso_by_scan = {}
    # for rawfile in rawfiles:
    #     filename = os.path.basename(rawfile)#.split('.')[0]
    #     rats = get_precursor_envelope_ratios(rawfile)
    #     percursor_iso_by_scan.update([((filename, num), ratio) for num, ratio in rats])
    percursor_iso_by_scan = dict(get_precursor_envelope_ratios(data))

    psms['expected_isotope_ratio'] = psms.apply(axis=1,
                                               func=lambda x: expected_ratio_for_peptide(x['peptide']))
    psms['observed_isotope_ratio'] = psms.apply(axis=1, func=lambda x: percursor_iso_by_scan.get(x['scannr'], 1))
    psms['isotope_ratio_error'] = abs(np.log(psms['observed_isotope_ratio']) -
                                      np.log(psms['expected_isotope_ratio']))

    # Replace nan with median
    psms['isotope_ratio_error'] = psms['isotope_ratio_error'].fillna(psms['isotope_ratio_error'].median())

    psms = psms.drop(['expected_isotope_ratio', 'observed_isotope_ratio'], axis=1)

    return psms 



def merge_features(*tables, merge_on = ['label', 'filename', 'scannr', 'expmass', 'calcmass', 'peptide', 'rec_mods', 'proteins',
                          'comet_rank', 'sage_rank', 'in_both',
                          'charge', 'num_proteins', 'peptide_len', 'missed_cleavages',
                          'delta_mass', 'abs_delta_mass']):
    merged = tables[0]
    for table in tables[1:]:
        merged = merged.merge(table, on=merge_on, how='inner')
    
    assert(all(len(merged)==len(x) for x in tables))
    # assert(set(merged.columns) == set(sum([x[1] for x in feature_files], [])+merge_on))

    return merged


def annotate_pin_file(psmfile, rawfile, featurefile, annotated_psmfile=None, coll_e = 32, instrument_name = 'Elite'):
    print("Annotating %s" % psmfile)

    data = mzFile(rawfile)

    psms = pd.read_csv(psmfile)
    # Neither comet now sage give us the actual MZ! Very impolite
    mz_rt_lookup = {x[2]:(x[0], x[1]) for x in data.scan_info()} # type: ignore
    # psms[['rt', 'mz']] = psms.scannr.apply(lambda x: mz_rt_lookup[x])
    psms['rt'] = psms.scannr.apply(lambda x: mz_rt_lookup[x][0])
    psms['mz'] = psms.scannr.apply(lambda x: mz_rt_lookup[x][1])

    psms = annotate_with_feature_info(psms, featurefile)

    rt_psms = get_rt_pred_deltas(psms.copy())
    frag_psms = get_frag_pred_deltas(psms.copy(), data, nce=coll_e, instrument=instrument_name)
    iso_psms = get_isotopic_envelope_delta(psms.copy(), data)
    print("Model running done.")
    
    merged_psms = merge_features(frag_psms, rt_psms, iso_psms)
    merged_psms = merged_psms[['label', 'filename', 'scannr', 'expmass', 'calcmass', 'mz', 'rec_mods', 
                               'comet_rank', 'sage_rank', 'in_both', 'wide_rank',
                               'charge', 'num_proteins', 'peptide_len', 'missed_cleavages',
                               'delta_mass', 'abs_delta_mass', 
                               'isotope_ratio_error', 'rt', 'rt_error', 'abs_rt_error', 'inc_frag_error', 'exc_frag_error',
                               'hyperscore', 'hyperscore_imputed', 'hyperscore_type', 
                               'chg_is_1', 'chg_is_2', 'chg_is_3', 'chg_is_4', 'chg_is_5', 'chg_is_6',
                               'peptide', 'proteins']]

    merged_psms['peptide'] = merged_psms.apply(axis=1, func=lambda x: x['peptide'] + (("___" + x['rec_mods']) if isinstance(x['rec_mods'], str) else '')) 
    merged_psms = merged_psms.drop('rec_mods', axis=1)

    
    # Check for mz that come up as 0.0; these are invalid scans for some reason.
    if (merged_psms.mz == 0.0).any():
        print("Removing %d scans with invalid MZ" % (merged_psms.mz == 0.0).sum())
        merged_psms = merged_psms[merged_psms.mz != 0.0]

    # Rrebuild specid since mokapot will expect it.
    merged_psms['specid'] = merged_psms.apply(axis=1, func=lambda x: '%s_%s' % (x['filename'], x['scannr']))

    if annotated_psmfile is None:
        annotated_psmfile = psmfile.rsplit('.', 1)[0] + '.annotated.pin'
    merged_psms.to_csv(annotated_psmfile, index=False, sep='\t')

    data.close()
    
    print("Done annotating %s" % annotated_psmfile)
    return annotated_psmfile


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('pin')
    parser.add_argument('raw')
    parser.add_argument('features')
    args = parser.parse_args()

    annotate_pin_file(args.pin, args.raw, args.features)
