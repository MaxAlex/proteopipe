import pandas as pd



def integrate_quant(psm_file, isobaric_datafiles):
    psms = pd.read_csv(psm_file)
    isobaric = pd.concat([pd.read_csv(x) for x in isobaric_datafiles])

    # Merge on both scan number and filename
    combined = pd.merge(psms, isobaric, left_on=['Scan', 'Filename'], right_on=['scan_id', 'Filename'])
    assert((combined['M/Z'] - combined['source_mz'] < 0.01).all())  # Ensuring that we're matching the right scans

    combined = combined.drop(columns=['scan_id', 'scan_index', 'source_ms1', 'source_mz', 'source_range',
                                      'source_charge', 'source_int', 'calibrant_int', 'calibrant_err'])
    combined = combined.drop([x for x in isobaric.columns if x.endswith('_mz')], axis=1)

    combined.to_csv(psm_file + '.w_isobaric.csv', index=False) 

    return combined

