#!/bin/bash 
comb_path="/home/siraj/new_Data_20220922/"  

echo " ================>>>>>>>>>>>>>>>>>> RNA-UV (4SU) (ECOLI) In house <<<<<<<<<<<<<<<============================"

#RNA_4SU_SAX
OpenNuXL -NuXL:presets "RNA-UV (4SU)" -in ${comb_path}RNA_UV_4SU_ECOLI_365-4thioU_SAX_mzML_corrected.mzML -out RNA_UV_4SU_ECOLI_AChernev_070119_Ecoli_365_4thioU_SAX_mzML_corrected.idXML -database Ecoli_SwPr_canon_20220310_4400seq_MQ2030contaminants.fasta -NuXL:scoring "slow" -report:top_hits 4 -modifications:variable_max_per_peptide 2 -NuXL:length 2 -threads 30  -modifications:variable "Oxidation (M)" -peptide:max_size 30 -peptide:min_size 05 -peptide:missed_cleavages 2 -filter "autotune" "idfilter" "optimize"
MzTabExporter -in RNA_UV_4SU_ECOLI_AChernev_070119_Ecoli_365_4thioU_SAX_mzML_corrected_perc_0.0100_XLs.idXML -out RNA_UV_4SU_ECOLI_AChernev_070119_Ecoli_365_4thioU_SAX_mzML_corrected_perc_0.0100_XLs.mzTab
MzTabExporter -in RNA_UV_4SU_ECOLI_AChernev_070119_Ecoli_365_4thioU_SAX_mzML_corrected_perc_0.0100_peptides.idXML -out RNA_UV_4SU_ECOLI_AChernev_070119_Ecoli_365_4thioU_SAX_mzML_corrected_perc_0.0100_peptides_df.mzTab
FeatureFinderIdentification -in ${comb_path}RNA_UV_4SU_ECOLI_365-4thioU_SAX_mzML_corrected.mzML -id RNA_UV_4SU_ECOLI_365-4thioU_SAX_mzML_corrected_perc_0.0100_XLs.idXML -out RNA_UV_4SU_ECOLI_365-4thioU_SAX_mzML_corrected_perc_0.0100_XLs.featureXML
FeatureFinderIdentification -in ${comb_path}RNA_UV_4SU_ECOLI_365-4thioU_SAX_mzML_corrected.mzML -id RNA_UV_4SU_ECOLI_365-4thioU_SAX_mzML_corrected_perc_0.0100_peptides.idXML -out RNA_UV_4SU_ECOLI_365-4thioU_SAX_mzML_corrected_perc_0.0100_peptides.featureXML
MzTabExporter -in RNA_UV_4SU_ECOLI_365-4thioU_SAX_mzML_corrected_perc_0.0100_XLs.featureXML -out RNA_UV_4SU_ECOLI_365-4thioU_SAX_mzML_corrected_perc_0.0100_XLs_featureXML.mzTab
MzTabExporter -in RNA_UV_4SU_ECOLI_365-4thioU_SAX_mzML_corrected_perc_0.0100_peptides.featureXML -out RNA_UV_4SU_ECOLI_365-4thioU_SAX_mzML_corrected_perc_0.0100_peptides_featureXML.mzTab

#RNA_4SU_QiSi
OpenNuXL -NuXL:presets "RNA-UV (4SU)" -in ${comb_path}RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected.mzML -out RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected.idXML -database Ecoli_SwPr_canon_20220310_4400seq_MQ2030contaminants.fasta -NuXL:scoring "slow" -report:top_hits 4 -modifications:variable_max_per_peptide 2 -NuXL:length 2 -threads 30  -modifications:variable "Oxidation (M)" -peptide:max_size 30 -peptide:min_size 05 -peptide:missed_cleavages 2 -filter "autotune" "idfilter" "optimize"
MzTabExporter -in RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected_perc_0.0100_XLs.idXML -out RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected_perc_0.0100_XLs.mzTab
MzTabExporter -in RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected_perc_0.0100_peptides.idXML -out RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected_perc_0.0100_peptides_df.mzTab
FeatureFinderIdentification -in ${comb_path}RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected.mzML -id RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected_perc_0.0100_XLs.idXML -out RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected_perc_0.0100_XLs.featureXML
FeatureFinderIdentification -in ${comb_path}RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected.mzML -id RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected_perc_0.0100_peptides.idXML -out RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected_perc_0.0100_peptides.featureXML
MzTabExporter -in RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected_perc_0.0100_XLs.featureXML -out RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected_perc_0.0100_XLs_featureXML.mzTab
MzTabExporter -in RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected_perc_0.0100_peptides.featureXML -out RNA_UV_4SU_ECOLI_Ecoli_365-4thioU_QiSi-2_mzML_corrected_perc_0.0100_peptides_featureXML.mzTab

echo " ================>>>>>>>>>>>>>>>>>> pRBSID_HF (4SU) <<<<<<<<<<<<<<<============================"

pRBSID_file_path="/home/siraj/Protein_XL/OpenNuXL_Compare/pRBSID_HF_RNA/"

OpenNuXL -NuXL:presets 'RNA-UV Extended (4SU)' -in ${pRBSID_file_path}mRBS_dTamicon30K_4SU_HCD_half_1.raw -out mRBS_dTamicon30K_4SU_HCD_half_1.idXML -database ${RBSID_file_path}Swiss_Human_uniprot_05032022.fasta -NuXL:scoring "slow" -report:top_hits 4 -modifications:variable_max_per_peptide 2 -NuXL:length 1 -threads 30 -modifications:variable "Oxidation (M)" "Carbamidomethyl (C)" -NET_executable ${Net_executable} -ThermoRaw_executable ${ThermoRaw_executable} -precursor:mass_tolerance 10.0 -filter "idfilter" "autotune" 
MzTabExporter -in mRBS_dTamicon30K_4SU_HCD_half_1_perc_0.0100_XLs.idXML -out mRBS_dTamicon30K_4SU_HCD_half_1_perc_0.0100_XLs.mzTab
MzTabExporter -in mRBS_dTamicon30K_4SU_HCD_half_1_perc_0.0100_peptides.idXML -out mRBS_dTamicon30K_4SU_HCD_half_1_perc_0.0100_peptides_df.mzTab
FeatureFinderIdentification -in ${pRBSID_file_path}mRBS_dTamicon30K_4SU_HCD_half_1.raw.mzML -id mRBS_dTamicon30K_4SU_HCD_half_1_perc_0.0100_XLs.idXML -out mRBS_dTamicon30K_4SU_HCD_half_1_perc_0.0100_XLs.featureXML
FeatureFinderIdentification -in ${pRBSID_file_path}mRBS_dTamicon30K_4SU_HCD_half_1.raw.mzML -id mRBS_dTamicon30K_4SU_HCD_half_1_perc_0.0100_peptides.idXML -out mRBS_dTamicon30K_4SU_HCD_half_1_perc_0.0100_peptides.featureXML
MzTabExporter -in mRBS_dTamicon30K_4SU_HCD_half_1_perc_0.0100_XLs.featureXML -out mRBS_dTamicon30K_4SU_HCD_half_1_perc_0.0100_XLs_featureXML.mzTab
MzTabExporter -in mRBS_dTamicon30K_4SU_HCD_half_1_perc_0.0100_peptides.featureXML -out mRBS_dTamicon30K_4SU_HCD_half_1_perc_0.0100_peptides_featureXML.mzTab

OpenNuXL -NuXL:presets 'RNA-UV Extended (4SU)' -in ${pRBSID_file_path}mRBS_dTamicon30K_4SU_HCD_half_2.raw -out mRBS_dTamicon30K_4SU_HCD_half_2.idXML -database ${RBSID_file_path}Swiss_Human_uniprot_05032022.fasta -NuXL:scoring "slow" -report:top_hits 4 -modifications:variable_max_per_peptide 2 -NuXL:length 1 -threads 30 -modifications:variable "Oxidation (M)" "Carbamidomethyl (C)" -NET_executable ${Net_executable} -ThermoRaw_executable ${ThermoRaw_executable} -precursor:mass_tolerance 10.0 -filter "autotune" "idfilter"
MzTabExporter -in mRBS_dTamicon30K_4SU_HCD_half_2_perc_0.0100_XLs.idXML -out mRBS_dTamicon30K_4SU_HCD_half_2_perc_0.0100_XLs.mzTab
MzTabExporter -in mRBS_dTamicon30K_4SU_HCD_half_2_perc_0.0100_peptides.idXML -out mRBS_dTamicon30K_4SU_HCD_half_2_perc_0.0100_peptides_df.mzTab
FeatureFinderIdentification -in ${pRBSID_file_path}mRBS_dTamicon30K_4SU_HCD_half_2.raw.mzML -id mRBS_dTamicon30K_4SU_HCD_half_2_perc_0.0100_XLs.idXML -out mRBS_dTamicon30K_4SU_HCD_half_2_perc_0.0100_XLs.featureXML
FeatureFinderIdentification -in ${pRBSID_file_path}mRBS_dTamicon30K_4SU_HCD_half_2.raw.mzML -id mRBS_dTamicon30K_4SU_HCD_half_2_perc_0.0100_peptides.idXML -out mRBS_dTamicon30K_4SU_HCD_half_2_perc_0.0100_peptides.featureXML
MzTabExporter -in mRBS_dTamicon30K_4SU_HCD_half_2_perc_0.0100_XLs.featureXML -out mRBS_dTamicon30K_4SU_HCD_half_2_perc_0.0100_XLs_featureXML.mzTab
MzTabExporter -in mRBS_dTamicon30K_4SU_HCD_half_2_perc_0.0100_peptides.featureXML -out mRBS_dTamicon30K_4SU_HCD_half_2_perc_0.0100_peptides_featureXML.mzTab

######################################## iTRAPP (4SU) ########################################################

iTRAPP_file_path="/home/siraj/Protein_XL/OpenNuXL_Compare/iTRAPP/"
TRAPP_file_path="/home/siraj/Protein_XL/OpenNuXL_Compare/TRAPP_/"
ThermoRaw_executable="/home/sachsenb/OMS/THIRDPARTY/All/ThermoRawFileParser/ThermoRawFileParser.exe"
NET_executable="/home/sachsenb/mono/bin/bin/mono"

OpenNuXL -NuXL:presets 'RNA-UV Extended (4SU)' -in ${iTRAPP_file_path}iTRAPP_lambda.raw -out iTRAPP_lambda.idXML -database ${iTRAPP_file_path}Swiss_S.cerevisiae_06032022.fasta -NuXL:scoring "slow" -report:top_hits 4 -NuXL:length 4 -threads 30 -modifications:variable "Oxidation (M)" -modifications:fixed 'Carbamidomethyl (C)' -NET_executable ${Net_executable} -ThermoRaw_executable ${ThermoRaw_executable} -filter "autotune"
MzTabExporter -in iTRAPP_lambda_perc_0.0100_XLs.idXML -out iTRAPP_lambda_perc_0.0100_XLs.mzTab
MzTabExporter -in iTRAPP_lambda_perc_0.0100_peptides.idXML -out iTRAPP_lambda_perc_0.0100_peptides_df.mzTab
FeatureFinderIdentification -in ${iTRAPP_file_path}iTRAPP_lambda.raw.mzML -id iTRAPP_lambda_perc_0.0100_XLs.idXML -out iTRAPP_lambda_perc_0.0100_XLs.featureXML
FeatureFinderIdentification -in ${iTRAPP_file_path}iTRAPP_lambda.raw.mzML -id iTRAPP_lambda_perc_0.0100_peptides.idXML -out iTRAPP_lambda_perc_0.0100_peptides.featureXML
MzTabExporter -in iTRAPP_lambda_perc_0.0100_XLs.featureXML -out iTRAPP_lambda_perc_0.0100_XLs_featureXML.mzTab
MzTabExporter -in iTRAPP_lambda_perc_0.0100_peptides.featureXML -out iTRAPP_lambda_perc_0.0100_peptides_featureXML.mzTab

#################################### iTRAPP:iTRAPP_no_lambda  ########################################################

OpenNuXL -NuXL:presets 'RNA-UV Extended (4SU)' -in ${iTRAPP_file_path}iTRAPP_no_lambda.raw -out iTRAPP_no_lambda.idXML -database ${iTRAPP_file_path}Swiss_S.cerevisiae_06032022.fasta -NuXL:scoring "slow" -report:top_hits 4 -NuXL:length 4 -threads 30 -modifications:variable "Oxidation (M)" -modifications:fixed 'Carbamidomethyl (C)' -NET_executable ${Net_executable} -ThermoRaw_executable ${ThermoRaw_executable} -filter "autotune"
MzTabExporter -in iTRAPP_no_lambda_perc_0.0100_XLs.idXML -out iTRAPP_no_lambda_perc_0.0100_XLs.mzTab
MzTabExporter -in iTRAPP_no_lambda_perc_0.0100_peptides.idXML -out iTRAPP_no_lambda_perc_0.0100_peptides_df.mzTab
FeatureFinderIdentification -in ${iTRAPP_file_path}iTRAPP_no_lambda.raw.mzML -id iTRAPP_no_lambda_perc_0.0100_XLs.idXML -out iTRAPP_no_lambda_perc_0.0100_XLs.featureXML
FeatureFinderIdentification -in ${iTRAPP_file_path}iTRAPP_no_lambda.raw.mzML -id iTRAPP_no_lambda_perc_0.0100_peptides.idXML -out iTRAPP_no_lambda_perc_0.0100_peptides.featureXML
MzTabExporter -in iTRAPP_no_lambda_perc_0.0100_XLs.featureXML -out iTRAPP_no_lambda_perc_0.0100_XLs_featureXML.mzTab
MzTabExporter -in iTRAPP_no_lambda_perc_0.0100_peptides.featureXML -out iTRAPP_no_lambda_perc_0.0100_peptides_featureXML.mzTab
