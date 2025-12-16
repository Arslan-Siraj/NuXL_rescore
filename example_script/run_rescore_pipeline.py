# !pip install rescore_pckg/NuXL_rescore
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

from pathlib import Path
from resource_manager import ensure_resources  

# Download data for run pipeline
import subprocess

# Download data for run pipeline
data_path = Path.cwd() / "data"

if not data_path.exists():
    data_path.mkdir(parents=True, exist_ok=True)

idxml_file = data_path / "Example_RNA_UV_XL.idXML"
if not idxml_file.exists():
    # Use subprocess to download via wget
    subprocess.run(
        [
            "wget",
            "-O",
            str(idxml_file),
            "https://raw.githubusercontent.com/Arslan-Siraj/nuxl-app/rescore/example-data/idXMLs/Example_RNA_UV_XL.idXML",
        ],
        check=True,
    )

mgf_file = data_path / "Example_RNA_UV_XL.mgf"
if not mgf_file.exists():
    # Use subprocess to download via wget
    subprocess.run(
        [
            "wget",
            "-O",
            str(mgf_file),
            "https://raw.githubusercontent.com/Arslan-Siraj/nuxl-app/rescore/example-data/mzML/Example_RNA_UV_XL.mgf",
        ],
        check=True,
    )


# Download and extract nuxl rescore resources if not already present
nuxl_resources_url = (
    "https://github.com/Arslan-Siraj/NuXL_rescore_resources/"
    "releases/download/0.0.1/nuxl_rescore_resource.zip"
)

nuxl_rescore_resource_dir = Path.cwd()  / "nuxl_rescore_resource"

if not nuxl_rescore_resource_dir.exists():
    nuxl_rescore_resource_dir.mkdir(parents=True, exist_ok=True)
    print("==> NuXL rescore resources not found. Downloading to:", nuxl_rescore_resource_dir)
    print("==> Directory created:", nuxl_rescore_resource_dir)
    ensure_resources(nuxl_rescore_resource_dir, nuxl_resources_url)
else:
    print("==> NuXL rescore resources found. Exist here:", nuxl_rescore_resource_dir)
    ensure_resources(nuxl_rescore_resource_dir, nuxl_resources_url)

# run as function
#from nuxl_rescore import run_pipeline
#run_pipeline(_id="/home/arslan/test_rescore/RNA_UV_Ecoli_S100_LB_bRfrac_9.idXML",
#    _calibration=nuxl_rescore_resource_dir / "nuxl_rescore_resource/calibration_data/RNA_All.csv",
#    _unimod=nuxl_rescore_resource_dir / "nuxl_rescore_resource/unimod/unimod_to_formula.csv",
#    _feat_config=nuxl_rescore_resource_dir / "nuxl_rescore_resource/features-config.json",
#    _rt_model="DeepLC",
#    _model_path=str(nuxl_rescore_resource_dir / "nuxl_rescore_resource/RT_deeplc_model/generic_model/full_hc_Train_RNA_All"),
#    _ms2pip=True,
#    _mgf_path = "/home/arslan/test_rescore/RNA_UV_Ecoli_S100_LB_bRfrac_9.mgf", 
#    _out = str(Path.cwd() / "rescore_out/"), # change directory if want
#    )

#how to add run with subprocess?
# or no need to import, run from CLI

import subprocess

run_from_cmd = [
    "nuxl_rescore", "run",
    "-id", str(idxml_file),
    "-calibration", str(nuxl_rescore_resource_dir / "nuxl_rescore_resource/calibration_data/RNA_All.csv"),
    "-unimod", str(nuxl_rescore_resource_dir / "nuxl_rescore_resource/unimod/unimod_to_formula.csv"),
    "-feat_config", str(nuxl_rescore_resource_dir / "nuxl_rescore_resource/features-config.json"),
    "-rt_model", "DeepLC",
    "-ms2pip",
    "-model_path", str(nuxl_rescore_resource_dir / "nuxl_rescore_resource/RT_deeplc_model/generic_model/full_hc_Train_RNA_All"),
    "-mgf", str(mgf_file),
    "-perc_exe", "percolator",
    "-perc_adapter", "PercolatorAdapter",
    "-out", Path.cwd() / "rescore_out/"
]

subprocess.run(run_from_cmd, check=True)
    
print("pipeline ran successfully")


