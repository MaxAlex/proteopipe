import os
import pandas as pd



# IsobaricQuant takes in an input hits file with the following columns:
# 
# Search ID: integer value to identify the search where the peptide was detected.
# Peptide ID: integer value that identifies the peptide
# Sequence: peptide sequence
# Reference: protein sequence reference
# Charge: peptide charge
# Start Scan: peptide MS2 Scan
# M/Z: peptide m/z
# 
# Here we're going to convert output from the Comet/Sage search to this format.
# 
# Then we're also going to take the output from IsobaricQuant, which is a table with the columns:
# 
# Peptide ID: peptide id introduced in the input file
# Quant ID: quantification id (unused field at 0)
# Fragment ion charge
# Fragment ion type (a, b…)
# Fragment ion position
# Fragment ion mz
# Fragment ion mz difference
# Fragment ion intensity
# Fragment ion matched (true or false)
# 
# And use that to re-annotate the Comet/Sage results with quantitation info.


def convert_psms_to_isobaricquant(psm_file, raw_file):
    filestem = os.path.basename(raw_file).split('.')[0]
    psms = pd.read_csv(psm_file, keep_default_na=False)

    # isobaricquant supports combined spectra, in which case multiple scans would be labeled with
    # the same peptide ID. We don't support that, so it's one-to-one between scans and "peptides"
    assert(psms['Spectrum_ID'].value_counts().max()==1)

    quant_rows = []
    for ind, psm in psms.iterrows():
        if psm['Filename'] != filestem:
            continue
        row = {'Search ID': str(hash(psm['Filename']))[:7],
               'Peptide ID': ind, 
               'Sequence': psm['Peptide'],
               'Reference': psm['Proteins'],
               'Charge': psm['Charge'],
               'Start Scan': psm['Scan'],
               'M/Z': psm['M/Z']}
        quant_rows.append(row)

    quant = pd.DataFrame(quant_rows, columns = ['Search ID', 'Peptide ID', 'Sequence', 'Reference', 'Charge', 'Start Scan', 'M/Z'])
    quant.to_csv(psm_file.rsplit('.', 1)[0] + '.isobaricquant_in.csv', index=False, header=False)


def _label_name_comp(lbl):
    return float(lbl.replace('C', '.0').replace('N', '.1'))

def convert_isobaricquant_to_psms(psm_file, quant_files):
    psms = pd.read_csv(psm_file)

    psms['peptideID'] = psms.index

    labels = sorted(set(pd.read_csv(quant_files[0])['labelName']), key=_label_name_comp)

    quant_rows = []
    for quant_file in quant_files:
        quant = pd.read_csv(quant_file)
        quant['peptideID'] = quant['peptideID'].astype(int)
    
        quant = quant[quant['labelName']!='zero']  # I think that's a bug?
        
        for (scan_num, pep_id), subquant in quant.groupby(['firstScan', 'peptideID']):
            subquant = subquant.sort_values('labelName', key=lambda x: x.apply(_label_name_comp)) # NB: redundant
            subquant.set_index('labelName', inplace=True)
            new_qrow = subquant.loc[labels]['maxIntensity'].rename(lambda x: 'TMT_' + x).to_dict()
            new_qrow['peptideID'] = pep_id
            new_qrow['Scan'] = scan_num
            # new_qrow['backgroundTMT'] = subquant['basePeakIntensity'].max()
            quant_rows.append(new_qrow)

    quant_table = pd.DataFrame(quant_rows, columns = ['peptideID', 'Scan']+labels) 
    psms_wiso = psms.merge(quant_table, on=['Scan', 'peptideID'], how='inner')
    assert(len(psms)==len(psms_wiso)==len(quant_table)), (len(psms), len(psms_wiso), len(quant_table))
    psms_wiso.to_csv(psm_file + '.w_isobaric.csv', index=False)




if __name__ == '__main__':
    # Run with 'generate_input' to generate the input file for isobaricquant for PSMs matching the
    # specified raw file, or with 'compile' to generate the final output file given the set of quant files
    # and the original PSM file.
    import argparse
    parser = argparse.ArgumentParser() 
    parser.add_argument('-p', '--psm', required=True)
    parser.add_argument('-r', '--raw')
    parser.add_argument('-q', '--quant', nargs='+', help='Quant files to compile into PSM file')
    parser.add_argument('-f', '--function')
    args = parser.parse_args()

    if args.function == 'generate_input':
        convert_psms_to_isobaricquant(args.psm, args.raw)
    elif args.function == 'compile':
        convert_isobaricquant_to_psms(args.psm, args.quant)
    else:
        raise Exception("Must specify either 'generate_input' or 'compile'")

