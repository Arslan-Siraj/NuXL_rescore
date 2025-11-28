# !pip install rescore_pckg/NuXL_rescore

from pathlib import Path
from resource_manager import ensure_resources  
from nuxl_rescore import run_pipeline

# Download and extract nuxl rescore resources if not already present
nuxl_resources_url = (
    "https://github.com/Arslan-Siraj/NuXL_rescore_resources/"
    "releases/download/0.0.1/nuxl_rescore_resource.zip"
)

nuxl_rescore_resource_dir = Path.cwd()  / "nuxl_rescore_resource"

if not nuxl_rescore_resource_dir.exists():
    nuxl_rescore_resource_dir.mkdir(parents=True, exist_ok=True)
    print("NuXL rescore resources not found. Downloading to:", nuxl_rescore_resource_dir)
    print("Directory created:", nuxl_rescore_resource_dir)
    ensure_resources(nuxl_rescore_resource_dir, nuxl_resources_url)
else:
    print("NuXL rescore resources not found. Downloading to:", nuxl_rescore_resource_dir)
    ensure_resources(nuxl_rescore_resource_dir, nuxl_resources_url)


# make output directory
output_dir = Path.cwd() / "rescore_out/"
if not output_dir.exists():
    output_dir.mkdir(parents=True, exist_ok=True)
    print("Output directory created at:", output_dir)


run_pipeline(_id="/home/arslan/nuxl_app/nuxl-app/example-data/idXMLs/Example_perc_0.0100_XLs.idXML",
    _calibration=nuxl_rescore_resource_dir / "nuxl_rescore_resource/calibration_data/RNA_All.csv",
    _unimod=nuxl_rescore_resource_dir / "nuxl_rescore_resource/unimod/unimod_to_formula.csv",
    _feat_config=nuxl_rescore_resource_dir / "nuxl_rescore_resource/nuxl_rescore/features-config.json",
    _rt_model="DeepLC",
    _model_path=str(nuxl_rescore_resource_dir / "nuxl_rescore_resource/RT_deeplc_model/generic_model/full_hc_Train_RNA_All"),
    _ms2pip=True,
    _mgf_path = "/home/arslan/test_rescore/RNA_UV_Ecoli_S100_LB_bRfrac_9.mgf", 
    _out = output_dir,
    )

print("pipeline ran successfully")
