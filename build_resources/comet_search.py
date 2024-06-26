import os, sys
import re
from multiplierz.mzReport import writer
from subprocess import call
import json
import pandas as pd

__all__ = ['CometSearch', 'perform_comet_search']



parameterPreamble = ['# comet_version 2023.01 rev. 2\n',
                     '# Comet MS/MS search engine parameters file.\n',
                     "# Everything following the '#' symbol is treated as a comment.\n"]

# # For instance, matches 'comet.2015012.win64.exe'
# cometRegex = r'comet\.20[0-9]{5}\.win(64|32)\.exe'
# cometExecutable = re.compile(cometRegex)

# from multiplierz import SettingsFile, myData
# import multiplierz.settings as settings
# cometPath = settings.get_comet()
try:
    cometPath = '/comet.linux.exe'
    assert(os.path.exists(cometPath)), cometPath
    example_par = '/comet.params.new'
    assert(os.path.exists(example_par)), example_par
except AssertionError:
    cometPath = '/home/max/biostuff/Comet/comet.exe'
    assert(os.path.exists(cometPath)), cometPath
    example_par = '/home/max/biostuff/Comet/comet.params.new'
    assert(os.path.exists(example_par)), example_par


# An error will be raised on invoking a CometSearch object if this path is
# missing or invalid.


    
    
enzymeDict = {'Arg_C': ['1', 'R', 'P'],
              'Asp_N': ['0', 'D', '-'],
              'CNBr': ['1', 'M', '-'],
              'Chymotrypsin': ['1', 'FWYL', 'P'],
              'Glu_C': ['1', 'DE', 'P'],
              'Lys_C': ['1', 'K', 'P'],
              'Lys_N': ['0', 'K', '-'],
              'No_enzyme': ['0', '-', '-'],
              'PepsinA': ['1', 'FL', 'P'],
              'Trypsin': ['1', 'KR', 'P'],
              'Trypsin/P': ['1', 'KR', '-']}
unitsDict = {'amu':'0', 'ppm':'2'}


fixedModTypesForward = {'A': 'add_A_alanine',
                        'B': 'add_B_user_amino_acid',
                        'C': 'add_C_cysteine',
                        'D': 'add_D_aspartic_acid',
                        'E': 'add_E_glutamic_acid',
                        'F': 'add_F_phenylalanine',
                        'G': 'add_G_glycine',
                        'H': 'add_H_histidine',
                        'I': 'add_I_isoleucine',
                        'J': 'add_J_user_amino_acid',
                        'K': 'add_K_lysine',
                        'L': 'add_L_leucine',
                        'M': 'add_M_methionine',
                        'N': 'add_N_asparagine',
                        'O': 'add_O_ornithine',
                        'P': 'add_P_proline',
                        'Q': 'add_Q_glutamine',
                        'R': 'add_R_arginine',
                        'S': 'add_S_serine',
                        'T': 'add_T_threonine',
                        'U': 'add_U_user_amino_acid',
                        'V': 'add_V_valine',
                        'W': 'add_W_tryptophan',
                        'X': 'add_X_user_amino_acid',
                        'Y': 'add_Y_tyrosine',
                        'Z': 'add_Z_user_amino_acid',
                        'pep_nterm':'add_Nterm_peptide',
                        'pep_cterm':'add_Cterm_peptide',
                        'prot_nterm':'add_Nterm_protein',
                        'prot_cterm':'add_Cterm_protein'}
fixedModTypesBackward = dict([(x,y) for y, x in list(fixedModTypesForward.items())])


# Does not include varmods args.
defaultParameters = {'activation_method': 'ALL',
                     'add_A_alanine': '0.0000',
                     'add_B_user_amino_acid': '0.0000',
                     'add_C_cysteine': '57.021464',
                     'add_Cterm_peptide': '0.0',
                     'add_Cterm_protein': '0.0',
                     'add_D_aspartic_acid': '0.0000',
                     'add_E_glutamic_acid': '0.0000',
                     'add_F_phenylalanine': '0.0000',
                     'add_G_glycine': '0.0000',
                     'add_H_histidine': '0.0000',
                     'add_I_isoleucine': '0.0000',
                     'add_J_user_amino_acid': '0.0000',
                     'add_K_lysine': '0.0000',
                     'add_L_leucine': '0.0000',
                     'add_M_methionine': '0.0000',
                     'add_N_asparagine': '0.0000',
                     'add_Nterm_peptide': '0.0',
                     'add_Nterm_protein': '0.0',
                     'add_O_ornithine': '0.0000',
                     'add_P_proline': '0.0000',
                     'add_Q_glutamine': '0.0000',
                     'add_R_arginine': '0.0000',
                     'add_S_serine': '0.0000',
                     'add_T_threonine': '0.0000',
                     'add_U_user_amino_acid': '0.0000',
                     'add_V_valine': '0.0000',
                     'add_W_tryptophan': '0.0000',
                     'add_X_user_amino_acid': '0.0000',
                     'add_Y_tyrosine': '0.0000',
                     'add_Z_user_amino_acid': '0.0000',
                     'allowed_missed_cleavage': '2',
                     'clear_mz_range': '0.0 0.0',
                     'clip_nterm_methionine': '0',
                     'database_name': '/some/path/db.fasta',
                     'decoy_prefix': 'DECOY_',
                     'decoy_search': '0',
                     'digest_mass_range': '600.0 5000.0',
                     'fragment_bin_offset': '0.4',
                     'fragment_bin_tol': '1.0005',
                     'isotope_error': '0',
                     'mass_type_fragment': '1',
                     'mass_type_parent': '1',
                     'max_fragment_charge': '3',
                     'max_precursor_charge': '6',
                     'max_variable_mods_in_peptide': '5',
                     'minimum_intensity': '0',
                     'minimum_peaks': '10',
                     'ms_level': '2',
                     'nucleotide_reading_frame': '0',
                     'num_enzyme_termini': '2',
                     'num_output_lines': '5',
                     'num_results': '100',
                     'num_threads': '0',
                     'output_outfiles': '0',
                     'output_pepxmlfile': '1',
                     'output_percolatorfile': '0',
                     'output_sqtfile': '0',
                     'output_sqtstream': '0',
                     'output_suffix': '',
                     'output_txtfile': '0',
                     'override_charge': '0',
                     'peff_format':'0',
                     'peff_obo':'-1',
                     'peptide_mass_tolerance': '3.00',
                     'peptide_mass_units': '0',
                     'precursor_charge': '0 0',
                     'precursor_tolerance_type': '1',
                     'print_expect_score': '1',
                     'remove_precursor_peak': '0',
                     'remove_precursor_tolerance': '1.5',
                     'require_variable_mod': '0',
                     'sample_enzyme_number': '1',
                     'scale_fragmentNL': '0',
                     'scan_range': '0 0',
                     'search_enzyme_number': '1',
                     'search_enzyme2_number': '0', # 0 == No second enzyme
                     'show_fragment_ions': '0',
                     'skip_researching': '1',
                     'spectrum_batch_size': '0',
                     'theoretical_fragment_ions': '1',
                     'use_A_ions': '0',
                     'use_B_ions': '1',
                     'use_C_ions': '0',
                     'use_NL_ions': '1',
                     'use_X_ions': '0',
                     'use_Y_ions': '1',
                     'use_Z_ions': '0',
                     'use_Z1_ions': '0',
                     'use_sparse_matrix': '1'}

# This needs '[COMET ENZYME INFO]' placed in front of it in the param file.
defaultEnzymes = {'0': {'active': 0, 'name': 'No_enzyme', 'p': '-', 'specificity': '-'},
                  '1': {'active': 1, 'name': 'Trypsin', 'p': 'P', 'specificity': 'KR'},
                  '2': {'active': 1, 'name': 'Trypsin/P', 'p': '-', 'specificity': 'KR'},
                  '3': {'active': 1, 'name': 'Lys_C', 'p': 'P', 'specificity': 'K'},
                  '4': {'active': 0, 'name': 'Lys_N', 'p': '-', 'specificity': 'K'},
                  '5': {'active': 1, 'name': 'Arg_C', 'p': 'P', 'specificity': 'R'},
                  '6': {'active': 0, 'name': 'Asp_N', 'p': '-', 'specificity': 'D'},
                  '7': {'active': 1, 'name': 'CNBr', 'p': '-', 'specificity': 'M'},
                  '8': {'active': 1, 'name': 'Glu_C', 'p': 'P', 'specificity': 'DE'},
                  '9': {'active': 1, 'name': 'PepsinA', 'p': 'P', 'specificity': 'FL'},
                  '10': {'active': 1, 'name': 'Chymotrypsin', 'p': 'P', 'specificity': 'FWYL'},
                  '11': {'active': 1, 'name': 'Glu_C/Tryp', 'p': 'P', 'specificity': 'DEKR'}}

defaultVarmods = [{'mass' : 15.994915,
                   'residues' : 'M',
                   'binary' : 0,
                   'max_mods_per_peptide' : 5,
                   'term_distance' : 0,
                   'N/C-term' : 0,
                   'required' : 0}]

def readVarmods():
    from multiplierz import load_mods
    modList = load_mods()
    
    varmods = {}
    for mod, site, mass in modList:
        varmods["%s(%s)" % (mod, site)] = str(mass)
    return varmods
        
            
#unimodDB = unimod.UnimodDatabase(os.path.join(myData, 'unimod.sqlite'))
from multiplierz.mass_biochem import unimod as unimodDB

# modLookup = unimodDB.get_pycomet_lookup()
# modLookup.update(readVarmods())
# modLookup.update({'iTRAQ8plex':{'add_Nterm_peptide':'304.20536', 'add_K_lysine':'304.20536'},
#                   'iTRAQ4plex':{'add_Nterm_peptide':'144.10206245', 'add_K_lysine':'144.10206245'},
#                   'Carbamidomethyl(C)': {'add_C_cysteine':'57.021464'},
#                   'Methylthio(C)': {'add_C_cysteine':'45.98772'}})
                            
toStandardPSMConversions = {'scan':'Query',
                            'num':'Peptide Rank',
                            'charge':'Charge',
                            'exp_neutral_mass':'Predicted mz',
                            'calc_neutral_mass':'Predicted mr',
                            'e-value':'Expectation Value',
                            'xcorr':'Cross-Correlation',
                            'delta_cn':'Delta CN',
                            'sp_score':'SP Score',
                            'ions_matched':'Ions Matched',
                            'ions_total':'Total Predicted Ions',
                            'plain_peptide':'Peptide Sequence',
                            'peptide_modifications':'Variable Modifications', # Old version.
                            'modified_peptide':'Variable Modifications', # New version.
                            'prev_aa':'Preceeding Residue',
                            'next_aa':'Following Residue',
                            'protein':'Accession Number',
                            'duplicate_protein_count':'Protein Count'}

       
#     varmod = {'mass' : bestType(mass),
#           'residues' : res,
#           'binary' : bestType(binary),
#           'max_mods_per_peptide' : bestType(maxMods),
#           'term_distance' : bestType(termDistance),
#           'N/C-term' : ncTerm,
#           'required' : req}
# self.varmods.append(varmod)
    
#         ist(modLookup.items())[0]
# ('Acetylation', {'add_H_histidine': 42.010565})
# list(modLookup.items())[200]
# ('deuterated methyl ester', {'add_K_lysine': 17.03448})
    

def render_to_fixmod(mod):
    if isinstance(mod, str):
        site, modmass = mod.split()
        site = site.strip(': ')
        modmass = modmass.strip(': ')
    elif isinstance(mod, tuple):
        site, modmass = mod
    else:
        raise IOError("Mod (%s) must be tuple or 'site: modification' string")
    
    return fixedModTypesForward[site], float(modmass)

def render_to_varmod(mod):
    if isinstance(mod, dict):
        return mod
    elif isinstance(mod, str) and os.path.exists(mod):
        mod_info = json.load(open(mod, 'r'))
        assert(set(mod_info.keys()) == {'mass', 'residues', 'binary', 'max_mods_per_peptide', 
                                        'term_distance', 'N/C-term', 'neutral_loss'})
        return 
    elif isinstance(mod, str):
        site, modmass = mod.split()
        site = site.strip(': ')
        modmass = modmass.strip(': ')
    elif isinstance(mod, tuple):
        site, modmass = mod
    else:
        raise IOError("Mod (%s) must be tuple or 'site: modification' string")
    modmass = float(modmass)

    # If given just a site+modmass pair, use default values for the rest.
    return {'mass':modmass,
            'residues':site,
            'binary':0,
            'max_mods_per_peptide':5,
            'term_distance':-1,
            'N/C-term':0,
            'neutral_loss':0.0}


def perform_comet_search(mgffile, database, fixed_mods = None, var_mods = None, **kwargs):
    assert os.path.exists(mgffile), "Input %s not found!" % mgffile
    assert os.path.exists(database), "Database %s not found!" % database
    
    try:
        performSearch, saveParameters = kwargs['runType']
    except KeyError:
        performSearch, saveParameters = True, False 
    if not (performSearch or saveParameters):
        print("Neither 'Perform Search' nor 'Save Parameter File' selected.  Aborting.")
        return
    
    searchObj = CometSearch()
    searchObj.database_name = database
    searchObj.allowed_missed_cleavage = missed_cleavages
    searchObj.peptide_mass_tolerance = prectol
    searchObj.peptide_mass_units = unitsDict[precunit]
    searchObj.fragment_bin_tol = fragbintol
    searchObj.fragment_bin_offset = fragbinoffset
    searchObj.num_threads = kwargs.get('NumThreads', 1)
    searchObj.spectrum_batch_size = kwargs.get('NumThreads', 0)
        
    if fixed_mods[:4] == '-FM=':
        fixed_mods = fixed_mods[4:]
    if var_mods[:4] == '-VM=':
        var_mods = var_mods[4:] 
    
    fixed_mod_list = [x.strip() for x in fixed_mods.split(',') if x]
    var_mod_list = [x.strip() for x in var_mods.split(',') if x]

    # It turns out that this function had just been left unfinished altogether.  Oops.
    searchObj.varmods = [render_to_varmod(x) for x in var_mod_list]
    searchObj.update([render_to_fixmod(x) for x in fixed_mod_list])


    
 
        
findmods = re.compile('\\[(\\d+.\\d+)\\]')
def convertVarmods(psm):
    peptide = psm['Peptide Sequence']
    cometvm = psm['Variable Modifications'][2:-2]
    mods = [(x.start(), x.group(1)) for x in findmods.finditer(cometvm)]
    assert all([cometvm[i-1].isalpha() for i, _ in mods])
    modstrs = ['%s%d: %s' % (cometvm[i-1], i, m) for i, m in mods]
    psm['Variable Modifications'] = '; '.join(modstrs)
    return psm
        
def format_text_report(reportfile, outputfile = None, mgffile = None, parameters = None,
                       most_rank = None, most_exp = None):
    """
    Renders a native Comet output .txt file into an mzReport-compatible and
    prettier .xlsx format.
    
    (Native .txt output is noncompatible mostly due to a space instead of 
    underscore in the 'modified peptide' column; hopefully that will be fixed
    soon.)
    """
    if most_rank:
        most_rank = int(most_rank)
    if most_exp:
        most_exp = float(most_exp)
    
    if mgffile:
        from multiplierz.mgf import parse_to_generator
        mgfgen = parse_to_generator(mgffile)
        queryToDesc = dict(enumerate(x['title'] for x in mgfgen), start = 1)
    else:
        queryToDesc = {}    
    
    columns = []
    rows = []
    report = open(reportfile, 'r')
    
    headeritems = next(report).split('\t')
    header = {'Program':headeritems[0],
              'Data':headeritems[1],
              'Search Run Time':headeritems[2],
              'Database':headeritems[3].strip()}
    columnline = next(report)
    
    # Fix for presumed bug; omit if this column title is changed in later Comet versions.
    columnline = columnline.replace('peptide\tmodifications', 'peptide_modifications')
    
    def tryNum(thing):
        try:
            return int(thing)
        except ValueError:
            try:
                return float(thing)
            except ValueError:
                return thing
    
    columns = [toStandardPSMConversions.get(x, x) for x in columnline.strip().split('\t')]
    for line in report:
        values = [tryNum(x.strip()) for x in line.split('\t')]
        row = dict(list(zip(columns, values)))
        row = convertVarmods(row)
        row['Spectrum Description'] = queryToDesc.get(row['Query'], 'Unknown')
        rows.append(row)
    
    report.close()
    if not outputfile:
        outputfile = '.'.join(reportfile.split('.')[:-1] + ['xlsx'])
    
    if outputfile.lower().endswith('xlsx') or outputfile.lower().endswith('xls'):
        headerwriter = writer(outputfile, columns = ['Program', 'Data',
                                                     'Search Run Time', 'Database'],
                                          sheet_name = 'Comet_Header')
        headerwriter.write(header)
        headerwriter.write(['', '', '', ''])
        if parameters:
            for setting, value in sorted(parameters.items()):
                headerwriter.write({'Program':setting, 'Data':value,
                                    'Search Run Time':'', 'Database':''})
        headerwriter.close()
    mainwriter = writer(outputfile, columns = ['Spectrum Description'] + columns, sheet_name = 'Data')
    for row in rows:
        if most_rank and row['Peptide Rank'] > most_rank:
            continue
        if most_exp and row['Expectation Value'] > most_exp:
            continue
        mainwriter.write(row)
    mainwriter.close()
    
    return outputfile
    
    

def format_xml_report(reportfile, outputfile, highest_rank = None, highest_exp = None):
    import xml.etree.ElementTree as ET
    
    if highest_rank:
        highest_rank = int(highest_rank)
    if highest_exp:
        highest_exp = int(highest_exp)
    
    root = ET.getroot()
    queries = [x for x in root[0]
               if x.tag == '{http://regis-web.systemsbiology.net/pepXML}spectrum_query']
    
     
    
    
    
    
    
    
    
def bestType(value):
    try:
        try:
            return int(value)
        except ValueError:
            return float(value)
    except ValueError:
        return str(value)
    
class CometSearch(dict):
    """
    Represents a comet parameters file, so that they can be manipulated via
    Python commands easily, and executes a Comet search using said parameters.
    
    If an existant parameters file is given, this is read in; if a non-existant
    parameters file name is given, defaults are used and saved to the given
    file name; if no file name is given, parameters are saved in a temp file
    only for use with Comet.
    
    .fields gives the available basic parameters; .enzymes contains the list
    of currently given enzymes; .varmods contains the current list of variable
    modifications.  Any of these three lists can be modified, and changes will
    be reflected in the resulting file.
    
    Any field can be accessed as a property of the object; for instance
    
    >settings = CometSearch()
    >settings.database_name
    'some/path/db.fasta'
    >settings.database_name = 'C:/real/path/to/fasta/database.fasta'
    """
    
    def __init__(self, parameterfile = None, database = None, save_parameter_file = False,
                 program_path = None):
        if program_path:
            assert os.path.exists(program_path), "Executable not found at %s" % program_path
            self.comet_path = program_path
        else:
            if not (cometPath and os.path.exists(cometPath)):
                raise RuntimeError("Comet executable not found at default location %s; "
                                   "update the multiplierz settings file to indicate "
                                   "your Comet installation." % cometPath)
            self.comet_path = cometPath
        
        self.file_name = parameterfile
        self.fasta_file = database

        self.fields = []
        self.enzymes = {}
        self.varmods = []
        self.fixmods = {}
        
        # Parameters for the run itself, used by run_comet_search().
        self.enzyme_selection = None
        
        if parameterfile:
            with open(parameterfile, 'r') as parameters:
                for line in parameters:
                    line = line.split('#')[0].strip()
                    
                    if not line:
                        continue
                    
                    try:
                        field, value = line.split('=')
                    except ValueError:
                        if line.strip() == '[COMET_ENZYME_INFO]':
                        #self.enzymes.append(line)                  
                            while True:
                                try:
                                    line = next(parameters)
                                    num, name, active, specificity, p = line.split()
                                except (IOError, ValueError, StopIteration):
                                    break
                                enzyme = {'name' : name,
                                          'active' : int(active),
                                          'specificity' : specificity,
                                          'p' : p
                                          }
                                self.enzymes[num.strip('.')] = enzyme
                            continue # Should be done; could be a return.
                        
                    if field[:12] == 'variable_mod' :
                        words = value.split()
                        if len(words) == 4:
                            mass, res, binary, maxMods = words
                            termDistance = -1
                            ncTerm = 0
                            req = 0
                        else:
                            mass, res, binary, maxMods, termDistance, ncTerm, req = words
                            
                        if float(mass) == 0.0:
                            continue
                        varmod = {'mass' : bestType(mass),
                                  'residues' : res,
                                  'binary' : bestType(binary),
                                  'max_mods_per_peptide' : bestType(maxMods),
                                  'term_distance' : bestType(termDistance),
                                  'N/C-term' : ncTerm,
                                  'required' : req}
                        self.varmods.append(varmod)

                        # Fixed mods are handled as generic modifications!
                    else:
                        self.fields.append(field.strip())
                        self[field.strip()] = bestType(value.strip())
                        

                        
        else:
            #self.update(defaultParameters)
            for field, value in list(defaultParameters.items()):
                self[field] = value
                self.fields.append(field)
            self.enzymes = defaultEnzymes
            self.varmods = defaultVarmods
        
    def write(self, file_name = None):
        assert file_name or self.file_name, "No output file specified!"
        
        assert len(self.varmods) <= 9, "Comet only allows up to 9 variable modifications."
        
        if (not file_name):
            file_name = self.file_name
        
        
        main_pars = []
        for field in self.fields:
            value = self[field]
            main_pars.append(('%s = %s\n' % (field, value)))
        main_pars.sort()
        with open(file_name, 'w') as parfile:
            for line in parameterPreamble:
                parfile.write(line)
                
            for line in main_pars:
                parfile.write(line)
            
            for num, mod in enumerate(self.varmods, start = 1):
                residues = mod['residues'].upper().replace('N-TERM', 'n').replace('C-TERM', 'c')
                residues = ''.join(x for x in residues if x.isalpha())
                parfile.write('variable_mod0%s = %.4f %s %s %s %s %s %s\n' % (num, mod['mass'], residues,
                                                                              mod['binary'], mod['max_mods_per_peptide'],
                                                                              mod['term_distance'], mod['N/C-term'],
                                                                              mod['required']))
                
            parfile.write('[COMET_ENZYME_INFO]\n')
            for num, enz in list(self.enzymes.items()):
                parfile.write('%s.\t%s\t%s\t%s\t%s\n' % (num, enz['name'], enz['active'],
                                                         enz['specificity'], enz['p']))
                
    
    def write_pretty(self, file_name=None):
        assert file_name or self.file_name, "No output file specified!"
        
        assert len(self.varmods) <= 9, "Comet only allows up to 9 variable modifications."
        
        if (not file_name):
            file_name = self.file_name
        
        # example_par = os.path.join(os.path.dirname(cometPath), 'comet.params.new')
        # call([cometPath, '-p'])
        assert(os.path.exists(example_par))

        seen_params = set()
        with open(example_par, 'r') as template:
            with open(file_name, 'w') as output:
                for line in template:
                    if '[COMET_ENZYME_INFO]' in line:
                        output.write('[COMET_ENZYME_INFO]\n')
                        for num, enz in list(self.enzymes.items()):
                            output.write('%s.\t%s\t%s\t%s\t%s\n' % (num, enz['name'], enz['active'],
                                                                    enz['specificity'], enz['p']))
                        break
                    if '#' in line:
                        text, comment = line.split('#', 1)
                        text = text.strip()
                        comment = '# ' + comment.strip()
                    else:
                        text = line.strip()
                        comment = ''
                    if not text:
                        output.write('%s\n' % comment)
                    else:
                        try:
                            key, default_val = text.split(' =')
                        except ValueError:
                            print(line)
                            print(text)
                            raise ValueError
                        
                        if key in self.keys():
                            output.write('%s = %s %s\n' % (key, self[key], comment))
                        elif 'variable_mod' in key:
                            modnum = int(key[-1])  # Only 9 mods in total supported, so
                            if modnum <= len(self.varmods):
                                mod = self.varmods[modnum-1]
                                residues = mod['residues']
                                residues = ''.join(x for x in residues if x.isalpha())
                                output.write('variable_mod0%s = %.4f %s %s %s %s %s %s\n' % (modnum, mod['mass'], residues,
                                                                                            mod['binary'], mod['max_mods_per_peptide'],
                                                                                            mod['term_distance'], mod['N/C-term'],
                                                                                            mod['required']))
                        else:
                            print("Using default value for missing parameter %s (%s)" % (key, default_val.strip()))
                            output.write('%s = %s %s\n' % (key, default_val.strip(), comment))
                        seen_params.add(key)

        unseen_params = seen_params - set(defaultParameters.keys())
        if unseen_params:
            print("Missing some (outdated?) parameters: %s" % (' '.join(unseen_params)))
        print("Wrote parameter file (%s)" % file_name)
        
    def run_search(self, data_file, outputfile = None, most_rank = None, most_exp = None, verbose = False):        
        self['digest_mass_range'] = '450.0 6000.0' # Remove this if this parameter gets added to the GUI!
        self['output_txtfile'] = '1'
        
        ext = data_file.split('.')[-1]
        if ext.lower() in ['raw', 'wiff', 'd']:
            from multiplierz.mgf import extract
            if verbose:
                print("Extracting to MGF...")
            data_file = extract(data_file)
            if verbose:
                print(("Extracted %s" % data_file))
                
        parfile = os.path.join(os.path.dirname(data_file), 'COMET.par.temp')
        print("Comet temp parameters location: %s" % parfile)
        # while os.path.exists(parfile): # Paranoia; avoid collisions.
        #     parfile += '.temp'
        if os.path.exists(parfile):
            print("Removing pre-existing %s" % parfile)
            os.remove(parfile)
        self.write_pretty(parfile)
        assert(os.path.exists(parfile))
        print(parfile)
        
        if not outputfile:
            outputfile = data_file + '.xlsx'
        
        try:
            expectedResultFile = data_file.rsplit('.', 1)[0] + '.txt'
            
            print('Initiating Comet search...')
            comet_command =[self.comet_path,
                            '-P' + parfile,
                            data_file]
            if self.fasta_file:
                comet_command.insert(2, '-D' + self.fasta_file)
            result = call(comet_command)
            print(('Comet search completed with return value %s' % result))    
            assert os.path.exists(expectedResultFile), "Comet failed to produce expected result file."
            
            return outputfile

            # if outputfile.split('.')[-1].lower() in ['xlsx', 'xls', 'mzd']:
            #     resultfile = format_report(expectedResultFile, outputfile, data_file,
            #                                parameters = dict(self),
            #                                most_rank = most_rank,
            #                                most_exp = most_exp)
            # else:
            #     resultfile = expectedResultFile
            
            # if outputfile:
            #     os.rename(resultfile, outputfile)
            # else:
            #     outputfile = resultfile
                
            # return outputfile
        
        finally:
            pass
            # os.remove(parfile)



def run_comet(mgf_file, fastafile, basic_config=True, 
              varmods = None, **kwargs):
    assert(mgf_file.endswith('.mgf'))
    pin_file = mgf_file.replace('.mgf', '.pin')
    txt_file = mgf_file.replace('.mgf', '.txt')
    if os.path.exists(pin_file) and os.path.exists(txt_file):
        print("%s and %s already exist" % (pin_file, txt_file))
        return pin_file, txt_file
    
    print("Running Comet on %s" % mgf_file)
    comet = CometSearch(database=fastafile)
    
    if basic_config:
        print("Using preset Comet parameters")
        comet['decoy_prefix'] = 'rev_'
        comet['fragment_bin_tol'] = '0.05'
        comet['num_output_lines'] = '100'
        comet['num_results'] = '5'
        
    comet.update(kwargs)
    assert(comet['decoy_prefix']=='rev_'), "Is the decoy prefix set correctly?"
    
    if varmods is not None:
        # Oxidized methionine is default; input varmods as empty list to remove.
        comet.varmods = varmods


    # Outputs hard-set because they're relied upon by the pipeline itself.
    comet['output_percolatorfile'] = '1'
    comet['output_txtfile'] = '1'
    comet.run_search(mgf_file)

    assert(os.path.exists(pin_file))
    assert(os.path.exists(txt_file))

    print("Checking %s" % txt_file)
    results = pd.read_csv(txt_file, sep='\t', skiprows=1)
    assert(len(results)>0)

    print("Done running Comet on %s" % mgf_file)
    return pin_file, txt_file


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Run Comet on an MGF file.')
    parser.add_argument('mgf_file', type=str, help='Input MGF file')
    parser.add_argument('fasta_file', type=str, help='FASTA file')
    parser.add_argument('--config')
    parser.add_argument('--varmods')
    args = parser.parse_args()

    
    comet_params = json.load(open(args.config))
    comet_varmods = json.load(open(args.varmods))
    run_comet(args.mgf_file, args.fasta_file, basic_config=True, varmods=comet_varmods, **comet_params)
