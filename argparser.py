import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    '-id', 
    type=str, 
    required=True,
    help="Input file (idXML format) path: where search for .preperc file if MS2PIP feature active",
    metavar='id'
    )

parser.add_argument(
    '-rt_model',
    type=str,
    required=False,
    default="DeepLC",
    help="if None no RT feature consider",
    metavar='rt_model'
)

parser.add_argument(
    '-calibration',
    type=str,
    required=False,
    default=os.getcwd()+"/calibration_data/RNA_All.csv",
    help="DeepLC calibration data path (.csv)",
    metavar='calibration'
)

parser.add_argument(
    '-model_path',
    type=str,
    required=False,
    default=os.getcwd()+"/RT_deeplc_model/generic_model/full_hc_Train_RNA_All",
    help="model path with name like full_hc_Train_RNA_All",
    metavar='model_path'
)

parser.add_argument(
    '-unimod',
    type=str,
    required=False,
    default=os.getcwd()+"/unimod/unimod_to_formula.csv",
    help="unimod/NuXL modification example /unimod/unimod_to_formula.csv",
    metavar='unimod'
)

parser.add_argument(
     '-out',
    type=str,
    required=False,
    default=os.getcwd()+"/Output/",
    help="output folder path ",
    metavar='out'
)

parser.add_argument(
     '-ms2pip',
    required=False,
    action='store_true',
    help="Extract ms2pip features {bool}",
)
  
parser.add_argument(
     '-ms2pip_rescore',
    required=False,
    action='store_true',
    help="Extract ms2pip_rescore features {bool}",
)
     
parser.add_argument(
     '-perc_exec',
    type=str,
    required=False,
    default="/home/ubuntu/home/siraj/Percolator/percolator",
    help="percolater executable path (Full path)",
    metavar='perc_exec'
)

parser.add_argument(
     '-perc_adapter',
    type=str,
    required=False,
    default="/home/ubuntu/Openms_Test/OpenMS-build/bin/PercolatorAdapter",
    help="percolater Adapter path (Full path)",
    metavar='perc_adapter'
)

parser.add_argument(
     '-feat_out',
    required=False,
    action='store_true',
    help="Write all extra feature output file at Out-folder (.csv)"
)

parser.add_argument(
     '-ms2pip_path',
    type= str,
    required=False,
    default= None,
    help="MS2PIP features file path (.csv file)",
    metavar='ms2pip_path'
)

parser.add_argument(
     '-ms2pip_rescore_path',
    type= str,
    required=False,
    default= None,
    help="MS2PIP rescore features file path (.pin file)",
    metavar='ms2pip_rescore_path'
)

parser.add_argument(
     '-mgf',
    type= str,
    required=False,
    default= None,
    help="Path for mgf file (.mgf file) need for MS2Rescore if results files not given",
    metavar='mgf'
)

parser.add_argument(
     '-peprec',
    type= str,
    required=False,
    default= None,
    help="Path for peprec (.peprec file) need for MS2Rescore if file not given generate from idXML",
    metavar='peprec'
)

parser.add_argument(
     '-feat_config',
    type= str,
    required=False,
    default=os.getcwd()+"/features-config.json",
    help="Path for feature config (.json file) need for features used for rescoring",
    metavar='feat_config'
)

parser.add_argument(
     '-entrap',
    required=False,
    action='store_true',
    help="Entrapment testing {bool}",
)

parser.add_argument(
     '-actual_db',
    type= str,
    required=False,
    default="None",
    help="Path for database (.fasta file) actual protein correspond to actual protocol",
    metavar='actual_db'
)

args = parser.parse_args()
