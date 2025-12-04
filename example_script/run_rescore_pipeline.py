# !pip install rescore_pckg/NuXL_rescore
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

from pathlib import Path
from resource_manager import ensure_resources  

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
    "-id", "/home/arslan/test_rescore/RNA_UV_Ecoli_S100_LB_bRfrac_9.idXML",
    "-calibration", str(nuxl_rescore_resource_dir / "nuxl_rescore_resource/calibration_data/RNA_All.csv"),
    "-unimod", str(nuxl_rescore_resource_dir / "nuxl_rescore_resource/unimod/unimod_to_formula.csv"),
    "-feat_config", str(nuxl_rescore_resource_dir / "nuxl_rescore_resource/features-config.json"),
    "-rt_model", "DeepLC",
    #"-ms2pip",
    "-model_path", str(nuxl_rescore_resource_dir / "nuxl_rescore_resource/RT_deeplc_model/generic_model/full_hc_Train_RNA_All"),
    "-mgf", "/home/arslan/test_rescore/RNA_UV_Ecoli_S100_LB_bRfrac_9.mgf",
    "-out", Path.cwd() / "rescore_out/",
]

subprocess.run(run_from_cmd, check=True)
    
print("pipeline ran successfully")


