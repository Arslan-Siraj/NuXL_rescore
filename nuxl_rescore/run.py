import pandas as pd
from pyopenms import *
from pkg_resources import get_distribution

from .argparser import build_parser
from .plotting import plot_weights_perc, comparison_PSMs, plot_FDR_plot
from .FDR_calculation import FDR_filtering_perc, run_percolator, FDR_unique_PSMs
from .Data_parser import peptide_ids_to_dataframe, read_pin_file, read_fasta, annotate_features
from .entrapment import entrapment_calculations
from .RT_features import predict_from_DeepLC, calculate_RTfeatures
from .ms2pip_features import Take_MS2PIP_features, Take_MS2PIP_rescore_features

def run_pipeline(_id=None, _calibration=None, _unimod=None, _feat_config=None,
    _model_path=None, _ms2pip=None, _ms2pip_path=None,
    _ms2pip_rescore=None, _ms2pip_rescore_path=None,
    _rt_model=None, _entrap=None, _actual_db=None, _out=None
):
    """
    explicit function arguments when called as Python API.
    """
    print("==> idXML Loading")
    protein_ids = []
    peptide_ids = []
    IdXMLFile().load(_id, protein_ids, peptide_ids)

    # -----------------------------
    # RT FEATURES
    # -----------------------------
    RT_predictions_feat_df = None

    if _rt_model is not None:
        print("==> RT columns extracting")

        RT_id_cols = peptide_ids_to_dataframe(peptide_ids)

        if _rt_model == "DeepLC":
            print("-->>> selected RT model: DeepLC")

            calibration_data = pd.read_csv(_calibration)

            RT_predictions = predict_from_DeepLC(RT_id_cols, calibration_data)
            RT_predictions_feat_df = calculate_RTfeatures(RT_predictions)

            print("Successfully extracted RT_features:", RT_predictions_feat_df.shape)
            RT_features_path = f"{_out}/RT_features.csv"
            RT_predictions_feat_df.to_csv(RT_features_path)

    if RT_predictions_feat_df is None:
        print("Warning: RT predictions not extracted. Use -rt_model DeepLC")

    # -----------------------------
    # MS2PIP INTENSITY FEATURES
    # -----------------------------
    MS2PIP_feat_df = None

    if _ms2pip:
        if _ms2pip_path is not None:
            MS2PIP_feat_df = pd.read_csv(_ms2pip_path)
        else:
            MS2PIP_path = Take_MS2PIP_features()
            MS2PIP_feat_df = pd.read_csv(MS2PIP_path)

        print("Successfully extracted MS2PIP_Features:", MS2PIP_feat_df.shape)

    else:
        print("Warning: MS2PIP intensity features disabled.")

    # -----------------------------
    # MS2RESCORE FEATURES
    # -----------------------------
    MS2PIP_rescore_feat_df = None

    if _ms2pip_rescore:
        if _ms2pip_rescore_path is not None:
            MS2PIP_rescore_feat_df = read_pin_file(_ms2pip_rescore_path)
        else:
            MS2PIP_rescore_feat_df = Take_MS2PIP_rescore_features()

        print("Successfully extracted MS2PIP_rescore Features:",
              MS2PIP_rescore_feat_df.shape)

    else:
        print("Warning: MS2PIP_rescore features disabled.")

    # -----------------------------
    # ANNOTATE FEATURES & STORE NEW idXML
    # -----------------------------
    print("==> writing features in idXML file")

    prot_ids, pep_ids, extra_feat_names = annotate_features(
        protein_ids, peptide_ids,
        RT_predictions_feat_df,
        MS2PIP_feat_df,
        MS2PIP_rescore_feat_df
    )

    out_file = _id.split("/")[-1]
    Feat_idXML_out_path = f"{_out}/updated_{out_file}"

    IdXMLFile().store(Feat_idXML_out_path, prot_ids, pep_ids)
    print("==> Updated idXML stored at:", Feat_idXML_out_path)

    # -----------------------------
    # PERCOLATOR
    # -----------------------------
    perc_result_file = run_percolator(_id, _perc_exec, _perc_adapter)
    FDR_perc_file = FDR_filtering_perc(perc_result_file + '.idXML')

    print("==> Percolator and FDR with extra features")
    Feat_perc_result_file = run_percolator(
        Feat_idXML_out_path, _perc_exec, _perc_adapter
    )

    plot_weights_perc(Feat_perc_result_file + '.weights', extra_feat_names)
    Feat_FDR_perc_file = FDR_filtering_perc(Feat_perc_result_file + '.idXML')

    # -----------------------------
    # ENTRAPMENT OR COMPARISON
    # -----------------------------
    if not _entrap:
        comparison_PSMs(
            Feat_perc_result_file + '_0.0100_XLs.idXML',
            perc_result_file + '_0.0100_XLs.idXML'
        )

        XL_all = perc_result_file + '_1.0000_XLs.idXML'
        XL_feat_all = Feat_perc_result_file + '_1.0000_XLs.idXML'

        plot_FDR_plot(XL_all, XL_feat_all)

    else:
        actual = read_fasta(_actual_db)
        unq = FDR_unique_PSMs(perc_result_file + '.idXML')
        unq_feat = FDR_unique_PSMs(Feat_perc_result_file + '.idXML')
        entrapment_calculations(unq, unq_feat, actual)
        


'''def main():

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
                run_pipeline()
        
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
                    run_pipeline()
            
            except Exception as e:
                print("An error occurred:", e)
                print("For help, about dependencies see requirements.txt")
              
        else:
            run_pipeline()

            
    '''