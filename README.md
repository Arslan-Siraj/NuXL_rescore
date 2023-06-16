# NuXL_rescore
This repository created for percolator base rescoring investigations of protein-nucleic acid crosslinking protocols with the addition of retention time prediction features from fine-tuned deepLC models and predictions of MSPIP base peak intensity features. <br /><br />
# Usage<br />
### predict retention time of RNA XL peptides and get .idXML percolator rescore output results<br />
run.py -id (.idxml) -model_path (model path from RT_deeplc_model) -calibration (calibration data path) -perc_exec (percolator executable) -perc_adapter (percolator_adapter) <br />

### addition of intensity features (if peak intensities already extracted)<br />
run.py -id (.idxml) -model_path (model path from RT_deeplc_model) -calibration (calibration data path) -perc_exec (percolator executable) -perc_adapter (percolator adapter) -ms2pip (store_true) -ms2pip_path (peak intensity file (.csv) path)<br />

### addition of intensity features (if peak intensities need to be extract)<br />
run.py -id (.idxml) -model_path (model path from RT_deeplc_model) -calibration (calibration data path) -perc_exec (percolator executable) -perc_adapter (percolator adapter) -ms2pip (store_true) -mgf (.mgf)<br />

### entrapment testing <br />
run.py -entrap (store_true)  -actual_db (actual database used during analysis)<br />

#### optional*<br />
-feat_config (path to .json feature configuration file)<br />
-out (path to output dir)<br />
-rt_model (path to new model, support of other model in future)<br />
-feat_out (extra feature file (.csv) will create in output directory)<br />
-peprec (for extract intensity features from file (.peprec), otherwise generate from .idXML and .mgf given file)<br />

#### help*<br />
run.py --help
