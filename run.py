import os
import re
import json
import pandas as pd
from pyopenms import *
from argparser import args
import matplotlib.pyplot as plt
from plotting import plot_weights_perc, plot_comparison_PSMs, comparison_PSMs, plot_FDR_plot
from FDR_calculation import FDR_filtering_perc, run_percolator, FDR_unique_PSMs
from Data_parser import peptide_ids_to_dataframe, read_pin_file, read_features_config, extract_intensities, read_fasta, annotate_features
from entrapment import entrapment_calculations

def process():
    print("-----Configuations-----")
    print(args)
    print("==> idXML Loading")
    protein_ids = []
    peptide_ids = []
    IdXMLFile().load(args.id, protein_ids, peptide_ids) 

    RT_predictions_feat_df = None
    if args.rt_model is not None:
        print("==> RT columns extracting")
        RT_id_cols = peptide_ids_to_dataframe(peptide_ids)
        
        if args.rt_model == "DeepLC":
            from RT_features import predict_from_DeepLC, calculate_RTfeatures
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
            from ms2pip_features import Take_MS2PIP_features
            MS2PIP_feat_df = Take_MS2PIP_features()
        print("Successfully extracted MS2PIP_Feature :", MS2PIP_feat_df.shape)
    else:
        print("Warning MS2PIP features (intensities) are not included, use -ms2pip True")

    MS2PIP_rescore_feat_df = None
    if args.ms2pip_rescore:
        if args.ms2pip_rescore_path is not None:
            MS2PIP_rescore_feat_df = read_pin_file(args.ms2pip_rescore_path)
        else:
            from ms2pip_features import Take_MS2PIP_rescore_features
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
    


if __name__ == "__main__":
    process()
    