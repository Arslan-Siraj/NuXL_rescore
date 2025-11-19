import pandas as pd
from pyopenms import *
from pkg_resources import get_distribution

from .argparser import args
from .plotting import plot_weights_perc, comparison_PSMs, plot_FDR_plot
from .FDR_calculation import FDR_filtering_perc, run_percolator, FDR_unique_PSMs
from .Data_parser import peptide_ids_to_dataframe, read_pin_file, read_fasta, annotate_features
from .entrapment import entrapment_calculations
from .RT_features import predict_from_DeepLC, calculate_RTfeatures
from .ms2pip_features import Take_MS2PIP_features, Take_MS2PIP_rescore_features

def process(id=None, calibration=None, unimod=None, feat_config=None,
            model_path=None, ms2pip=None, ms2pip_path=None,
            ms2pip_rescore=None, ms2pip_rescore_path=None,
            rt_model=None, entrap=None, actual_db=None, out=None):

    # If arguments are provided explicitly, override CLI parser
    if id is not None:
        args.id = id
    if calibration is not None:
        args.calibration = calibration
    if unimod is not None:
        args.unimod = unimod
    if feat_config is not None:
        args.feat_config = feat_config
    if model_path is not None:
        args.model_path = model_path
    if ms2pip is not None:
        args.ms2pip = ms2pip
    if ms2pip_path is not None:
        args.ms2pip_path = ms2pip_path
    if ms2pip_rescore is not None:
        args.ms2pip_rescore = ms2pip_rescore
    if ms2pip_rescore_path is not None:
        args.ms2pip_rescore_path = ms2pip_rescore_path
    if rt_model is not None:
        args.rt_model = rt_model
    if entrap is not None:
        args.entrap = entrap
    if actual_db is not None:
        args.actual_db = actual_db
    if out is not None:
        args.out = out

    print("==> idXML Loading")
    protein_ids = []
    peptide_ids = []
    IdXMLFile().load(args.id, protein_ids, peptide_ids) 

    RT_predictions_feat_df = None
    if args.rt_model is not None:
        print("==> RT columns extracting")
        RT_id_cols = peptide_ids_to_dataframe(peptide_ids)
        
        if args.rt_model == "DeepLC":
            print("-->>> selected RT model DeepLC") 
            calibration_data = pd.read_csv(args.calibration)
            RT_predictions = predict_from_DeepLC(RT_id_cols, calibration_data)
            RT_predictions_feat_df = calculate_RTfeatures(RT_predictions)
            print("Successfully extracted RT_features: ", RT_predictions_feat_df.shape)
            RT_predictions_feat_df.to_csv(args.out+"RT_features.csv")

    if RT_predictions_feat_df is None:
        print("Warning RT_predictions not extracted, use -rt_model DeepLC option")

    MS2PIP_feat_df = None 
    if args.ms2pip:
        if args.ms2pip_path is not None:
            MS2PIP_feat_df =  pd.read_csv(args.ms2pip_path)
        else:
            MS2PIP_path = Take_MS2PIP_features()
            MS2PIP_feat_df =  pd.read_csv(MS2PIP_path)
            
        print("Successfully extracted MS2PIP_Feature :", MS2PIP_feat_df.shape)
    else:
        print("Warning MS2PIP features (intensities) are not included, use -ms2pip True")

    MS2PIP_rescore_feat_df = None
    if args.ms2pip_rescore:
        if args.ms2pip_rescore_path is not None:
            MS2PIP_rescore_feat_df = read_pin_file(args.ms2pip_rescore_path)
        else:
            MS2PIP_rescore_feat_df = Take_MS2PIP_rescore_features() 

        print("Successfully extracted MS2PIP_rescore Features :", MS2PIP_rescore_feat_df.shape)     
    else:
        print("Warning MS2PIP rescore features are not included, use -ms2pip_rescore True")
    

    print("==> writing features in idXML file")
    prot_ids, pep_ids, extra_feat_names = annotate_features(protein_ids, peptide_ids, RT_predictions_feat_df, MS2PIP_feat_df, MS2PIP_rescore_feat_df)

    out_file_ = (args.id).split("/")
    Feat_idXML_out_path = args.out+"updated_"+out_file_[len(out_file_)-1]
    IdXMLFile().store(Feat_idXML_out_path, prot_ids, pep_ids)
    print("==>extra featured idXML stored at: ", Feat_idXML_out_path) 
    
    perc_result_file = run_percolator(args.id, args.perc_exec , args.perc_adapter)
    FDR_perc_file = FDR_filtering_perc(perc_result_file+'.idXML')
    
    print("==>Percolator and FDR calculations with extra features")
    Feat_perc_result_file = run_percolator(Feat_idXML_out_path, args.perc_exec , args.perc_adapter)
    plot_weights_perc(Feat_perc_result_file+'.weights', extra_feat_names)
    Feat_FDR_perc_file = FDR_filtering_perc(Feat_perc_result_file+'.idXML')

    if args.entrap == False:
        comparison_PSMs(Feat_perc_result_file+ '_0.0100_XLs.idXML', perc_result_file+ '_0.0100_XLs.idXML')
        XL_100_file = perc_result_file+'_1.0000_XLs.idXML'
        XL_100_feat_file = Feat_perc_result_file+'_1.0000_XLs.idXML'
        plot_FDR_plot(XL_100_file, XL_100_feat_file)
        
    if args.entrap:
        actual_prot = read_fasta(args.actual_db)
        un_100_XL_file = FDR_unique_PSMs(perc_result_file+'.idXML')
        un_100_XL_feat_file = FDR_unique_PSMs(Feat_perc_result_file+'.idXML')
        entrapment_calculations(un_100_XL_file, un_100_XL_feat_file, actual_prot)
        


def main():

    print("-----Configuation-----")
    for attr, value in vars(args).items():
        print(f"{attr}: {value}")

    if args.ms2pip and args.ms2pip_rescore : 
        print("Error! please select ms2rescore features or ms2pip intensity features or combine features like e-g RT+intensities or RT+ms2rescore")
    
    else :
        if args.ms2pip and args.ms2pip_path is None:
            ms2pip_curr_version = get_distribution("ms2pip").version
            ms2pip_desire_version = "3.11.0"
            print("ms2pip version: ", ms2pip_curr_version)
            if ms2pip_curr_version != ms2pip_desire_version :
                print("Error! ms2pip version ", ms2pip_desire_version , "required ", "For help, about dependencies see requirements.txt")
                print("Try pip install ms2pip==3.11.0")
            else: 
                process()
        
        elif args.ms2pip_rescore and args.ms2pip_rescore_path is None:
            try:
                ms2pip_curr_version = get_distribution("ms2pip").version
                ms2rescore_curr_version = get_distribution("ms2rescore").version
                ms2pip_desire_version = "4.0.0.dev1"
                print("ms2pip version: ", ms2pip_curr_version)
                print("ms2rescore version: ", ms2rescore_curr_version)
                if int(ms2pip_curr_version[0]) != int(4):
                    print("Error! ms2pip desire version 4.0.0.dev1, .., 4.0.0.dev5")
                    print("For help, about dependencies see requirements.txt")
                elif int(ms2rescore_curr_version[0]) < int(3):
                    print("Error! ms2rescore desire version 3.0.b4")
                    print("For help, about dependencies see requirements.txt")
                else:
                    process()
            
            except Exception as e:
                print("An error occurred:", e)
                print("For help, about dependencies see requirements.txt")
              
        else:
            process()

            
    