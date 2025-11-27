
from nuxl_rescore import run_pipeline, readAndProcessIdXML


readed_df = readAndProcessIdXML("/home/arslan/test_rescore/RNA_UV_Ecoli_S100_LB_bRfrac_9.idXML", top=1)
print(readed_df.head())

print("package loaded successfully")
run_pipeline(id="/home/arslan/nuxl_app/nuxl-app/example-data/idXMLs/Example_perc_0.0100_XLs.idXML",
    calibration="/home/arslan/rescore_pckg/NuXL_rescore/calibration_data/RNA_All.csv",
    unimod="/home/arslan/rescore_pckg/NuXL_rescore/unimod/unimod_to_formula.csv",
    feat_config="/home/arslan/rescore_pckg/NuXL_rescore/nuxl_rescore/features-config.json",
    model_path="/home/arslan/rescore_pckg/NuXL_rescore/RT_deeplc_model/generic_model/full_hc_Train_RNA_All")

print("pipeline ran successfully")
