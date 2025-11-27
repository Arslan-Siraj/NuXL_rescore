import os
import argparse
from .run import run_pipeline

def build_parser():
    parser = argparse.ArgumentParser(
        description="NuXL rescoring pipeline"
    )

    parser.add_argument(
        '-id', 
        type=str,
        required=False,
        help="Input file (idXML format) …",
        metavar='id'
    )
    parser.add_argument(
        '-rt_model',
        type=str,
        required=False,
        default="DeepLC",
        metavar='rt_model'
    )
    parser.add_argument(
        '-calibration',
        type=str,
        required=False,
        default=os.getcwd()+"/calibration_data/RNA_All.csv",
        metavar='calibration'
    )
    parser.add_argument(
        '-model_path',
        type=str,
        required=False,
        default=os.getcwd()+"/RT_deeplc_model/generic_model/full_hc_Train_RNA_All",
        metavar='model_path'
    )
    parser.add_argument(
        '-unimod',
        type=str,
        required=False,
        default=os.getcwd()+"/unimod/unimod_to_formula.csv",
        metavar='unimod'
    )
    parser.add_argument(
        '-out',
        type=str,
        required=False,
        default=os.getcwd()+"/Output/",
        metavar='out'
    )
    parser.add_argument('-ms2pip', required=False, action='store_true')
    parser.add_argument('-ms2pip_rescore', required=False, action='store_true')
    parser.add_argument(
        '-perc_exec',
        type=str,
        required=False,
        default="/home/ubuntu/home/siraj/Percolator/percolator",
        metavar='perc_exec'
    )
    parser.add_argument(
        '-perc_adapter',
        type=str,
        required=False,
        default="/home/ubuntu/Openms_Test/OpenMS-build/bin/PercolatorAdapter",
        metavar='perc_adapter'
    )
    parser.add_argument('-feat_out', required=False, action='store_true')
    parser.add_argument('-ms2pip_path', type=str, required=False, default=None)
    parser.add_argument('-ms2pip_rescore_path', type=str, required=False, default=None)
    parser.add_argument('-mgf', type=str, required=False, default=None)
    parser.add_argument('-peprec', type=str, required=False, default=None)
    parser.add_argument(
        '-feat_config',
        type=str,
        required=False,
        default=os.getcwd()+"/features-config.json",
        metavar='feat_config'
    )
    parser.add_argument('-entrap', required=False, action='store_true')
    parser.add_argument(
        '-actual_db',
        type=str,
        required=False,
        default="None",
        metavar='actual_db'
    )

    return parser

def run_from_CLI(args):
      
    print("-----Configuation-----")
    for attr, value in vars(args).items():
        print(f"{attr}: {value}")
   
    run_pipeline(_id=args.id, _calibration=args.calibration, _unimod=args.unimod, _feat_config=args.feat_config, _feat_out=args.feat_out,
    _model_path=args.model_path, _ms2pip=args.ms2pip, _ms2pip_path=args.ms2pip_path,
    _ms2pip_rescore=args.ms2pip_rescore, _ms2pip_rescore_path=args.ms2pip_rescore_path,
    _rt_model=args.rt_model, _entrap=args.entrap, _actual_db=args.actual_db, _out=args.out, _peprec_file=args.peprec, _mgf_file=args.mgf
    )
