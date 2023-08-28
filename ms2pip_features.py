from psm_utils.io.peptide_record import PeptideRecordReader
from collections import defaultdict
from psm_utils.io import peptide_record
from psm_utils.io import write_file
from tqdm import tqdm
import pandas as pd
from argparser import args
from psm_utils.io import convert
from Data_parser import read_pin_file

CONFIG = {  'ms2rescore': 
            {   'tmp_path': '', 
                'spectrum_path': '', 
                'output_path': '',
                'psm_file': '',
                'psm_id_pattern': None, 
                'spectrum_id_pattern': ".*_(controllerType=0 controllerNumber=1 scan=[0-9]+)_.*", #to take the predictions of all rank PSMs
                'processes': 32,
                'num_cpu': 4}, 
            'ms2pip': {
                'model': 'HCD', 
                'frag_error': 0.02}}

def initilize_CONFIG(mgf_file : str, out_pin_file : str, psm_file: str):
    """
    Update config file according to args
    """
    global CONFIG
    CONFIG = {  'ms2rescore': 
                {   'tmp_path': '', 
                    'spectrum_path': mgf_file, 
                    'output_path': out_pin_file,
                    'psm_file': psm_file,
                    'psm_id_pattern': None, 
                    'spectrum_id_pattern': ".*_(controllerType=0 controllerNumber=1 scan=[0-9]+)_.*", 
                    'processes': 32,
                    'num_cpu': 32}, 
                'ms2pip': {
                    'model': 'HCD', 
                    'frag_error': 0.02}}
    return CONFIG

def get_psm_list(inputfile):
    """
    Read from psm_list from peprec
    """
    id_file = PeptideRecordReader(inputfile)
    id_file.filename = inputfile

    return id_file.read_file()

def Take_ms2pip_rescore_features(psm_list):
    """
    Extract MSPIP rescore features (as PSMS_list)
    """
    print("Ms2PIP-rescore features gathering")
    from ms2rescore.feature_generators import ms2pip
    
    n_duplicate = defaultdict(lambda: 1)
    number_duplicates_per_spec = 1

    indices_list = []
    for spec_id in psm_list["spectrum_id"]:
        indices_list.append(n_duplicate[spec_id])
        if n_duplicate[spec_id] > number_duplicates_per_spec:
            number_duplicates_per_spec = n_duplicate[spec_id]

        n_duplicate[spec_id] += 1

    for i in range(1,number_duplicates_per_spec+1):
        
        ms2pip.MS2PIPFeatureGenerator(CONFIG,processes=CONFIG["ms2rescore"]["processes"], spectrum_path = CONFIG["ms2rescore"]["spectrum_path"], spectrum_id_pattern = CONFIG["ms2rescore"]["spectrum_id_pattern"]).add_features(psm_list[[True if x == i else False for x in indices_list ]])      
            
def write_pin(psm_list):
    """
    write MS2PIP rescore features in pin file
    """
    write_file(
    psm_list,
    filename=CONFIG['ms2rescore']["output_path"],
    filetype="percolator",
    style="pin",
    feature_names=psm_list[0].rescoring_features.keys(),
        )

def unique(list1):
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return unique_list

def Take_ms2pip_features(psm_list, out_file):
    """
    Extract MSPIP features (DataFrame of intensities), will furthur extract intensity from output file
    """       
    print("update  CONFIG for MS2PIP feature-----")
    #"ptm": config_up._get_modification_config(psm_list),
    from ms2pip.ms2pipC import MS2PIP

    CONFIG["ms2pip"].update(
    {
        "ptm": [
             "Oxidation,15.994915,opt,M",
             "Carbamidomethyl,57.021464,opt,C",
        ],
        "sptm": [],
        "gptm": [],
    }
    ) 

    print("Extracting_ms2pip_features------")
    n_duplicate = defaultdict(lambda: 1)
    number_duplicates_per_spec = 1

    indices_list = []
    spec_id_unique_ = unique(psm_list["spectrum_id"])
    for spec_id in psm_list["spectrum_id"]:
        indices_list.append(n_duplicate[spec_id])
        if n_duplicate[spec_id] > number_duplicates_per_spec:
            number_duplicates_per_spec = n_duplicate[spec_id]

        n_duplicate[spec_id] += 1
    
    Dataframe_list = []
    x=0
    for i in range(1,number_duplicates_per_spec+1):
        psm_list_ = psm_list[[True if x == i else False for x in indices_list]]
        MSPIP_feature = MS2PIP(
                        peptide_record.to_dataframe(psm_list_),
                        spec_file=CONFIG["ms2rescore"]["spectrum_path"],
                        spectrum_id_pattern=CONFIG["ms2rescore"]["spectrum_id_pattern"],
                        params=CONFIG,
                        return_results=True,
                        num_cpu=CONFIG["ms2rescore"]["num_cpu"],
                    )
        pred_and_emp = MSPIP_feature.run()
        spec_id_unique = list(pred_and_emp["spec_id"].unique())
        spec_id_ls = []
        charge_ls = []
        ion_series_ls = []
        ion_series_mz_ls = []
        ion_series_pred_ls = []
        ion_series_targ_ls = []
        rank_series = []
        for spec_id in spec_id_unique:
            spec_id_ls.append(spec_id)
            rank_series.append(str(x))
            all_ions = pred_and_emp.loc[pred_and_emp["spec_id"]==spec_id]
            charge_ls.append(list(all_ions["charge"].unique()))
            ion_series_ls.append(list(all_ions["ion"].astype(str) + all_ions["ionnumber"].astype(str)))
            ion_series_mz_ls.append(list(all_ions["mz"]))
            ion_series_pred_ls.append(list(all_ions["prediction"]))
            ion_series_targ_ls.append(list(all_ions["target"]))

        columns = ["spec_id", "rank", "ions_series", "ions_charge", "ions_mz","ions_pred", "ions_targ"]
        data = {
            "spec_id":spec_id_ls, 
            "rank": rank_series,
            "ions_series":ion_series_ls,
            "ions_charge":charge_ls, 
            "ions_mz":ion_series_mz_ls,
            "ions_pred":ion_series_pred_ls,
             "ions_targ":ion_series_targ_ls
             }
        
        Dataframe_list.append(pd.DataFrame(data, columns=columns))
        x=x+1
    
    concat_df = pd.concat(Dataframe_list, join="inner")
    
    final_PSMs = []
    for spec in spec_id_unique_:
        rows_with_rank = concat_df.loc[concat_df["spec_id"]==str(spec)]
        #rows_with_rank_ = rows_with_rank.drop(["rank"], axis=1)
        final_PSMs.append(rows_with_rank)

    final_psms_df = pd.concat(final_PSMs, join="inner")
    file_name_ = out_file.split('/')
    file_name = file_name_[len(file_name_)-1]
    final_psms_df.to_csv(args.out + file_name+"_MSPIP.csv")
    print("MSPIP features written at: ",args.out + file_name+"_MSPIP.csv")    
    return args.out + file_name+"_MSPIP.csv"

def Take_MS2PIP_features():
    """
    Extract MSPIP features (DataFrame of intensities)
    Extract MSPIP rescore features (.pin file)
    """
    file_id = (args.id).split('.')
    file_id_ = file_id[0].split('/')
  
    out_pin_file = args.out + file_id_[len(file_id_)-1] +'.pin'
    if args.peprec is None:
        print("converting .idXML to .peprec format")
        peprec_file = args.out + file_id_[len(file_id_)-1] + '.peprec'
        convert(args.id, peprec_file)
        print(".peprec written at: ", peprec_file)
        args.peprec = peprec_file

    ms2pip_features_out = None
    if args.mgf is not None:
        CONFIG = initilize_CONFIG(args.mgf, out_pin_file , args.peprec)
        print("Initialized MS2PIP CONFIG----\n", CONFIG)
        psm_list = get_psm_list(CONFIG["ms2rescore"]["psm_file"])
        ms2pip_features_out = Take_ms2pip_features(psm_list, file_id[0])
    else:
        print("Error: unable_initialized please provide (.mgf) file")

    return ms2pip_features_out

def Take_MS2PIP_rescore_features():
    """
    Extract MSPIP features (DataFrame of intensities)
    Extract MSPIP rescore features (.pin file)
    """
    file_id = (args.id).split('.')
    file_id_ = file_id[0].split('/')
    out_pin_file = args.out + file_id_[len(file_id_)-1] +'.pin'
    if args.peprec is None:
        print("converting .idXML to .peprec format")
        peprec_file = args.out + file_id_[len(file_id_)-1] + '.peprec'
        convert(args.id, peprec_file)
        print(".peprec written at: ", peprec_file)
        args.peprec = peprec_file
        
    ms2pip_rescore_feat = None 
    if args.mgf is not None:
        CONFIG = initilize_CONFIG(args.mgf, out_pin_file , args.peprec)
        print("Initialized MS2PIP CONFIG----\n", CONFIG)
        psm_list = get_psm_list(CONFIG["ms2rescore"]["psm_file"])
        psm_list["rescoring_features"] =  [{} for _ in range(len(psm_list))] 
        Take_ms2pip_rescore_features(psm_list)
        write_pin(psm_list)
        ms2pip_rescore_feat = read_pin_file(out_pin_file)
        print("MSPIP rescore features written at: ", out_pin_file)
    else:
        print("Error: unable_initialized please provide .mgf file")

    return ms2pip_rescore_feat 
