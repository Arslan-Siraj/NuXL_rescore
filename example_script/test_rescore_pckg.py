
from nuxl_rescore import run_pipeline

print("package loaded successfully")
run_pipeline(_id="/home/arslan/nuxl_app/nuxl-app/example-data/idXMLs/Example_perc_0.0100_XLs.idXML",
    _calibration="/home/arslan/rescore_pckg/NuXL_rescore/calibration_data/RNA_All.csv",
    _unimod="/home/arslan/rescore_pckg/NuXL_rescore/unimod/unimod_to_formula.csv",
    _feat_config="/home/arslan/rescore_pckg/NuXL_rescore/nuxl_rescore/features-config.json",
    _rt_model="DeepLC",
    _model_path="/home/arslan/rescore_pckg/NuXL_rescore/RT_deeplc_model/generic_model/full_hc_Train_RNA_All",
    _ms2pip=True,
    _mgf_path = "/home/arslan/test_rescore/RNA_UV_Ecoli_S100_LB_bRfrac_9.mgf", 
    _out = "/home/arslan/test_rescore/test_out/",
    )

print("pipeline ran successfully")
