#!/bin/bash 

echo " ================>>>>>>>>>>>>>>>>>> RNA_NM <<<<<<<<<<<<<<<============================"

python correlation_analysis.py -in_peps "/home/ubuntu/mspip_all_new/NM_S30_7_MSPIP_test/updated_MRaabe_LW_091221_171221_Expl2_XL_Ecoli_NM_S30_bRP_rep1_7_perc_0.1000_peptides.idXML" \
        -in_XLs "/home/ubuntu/mspip_all_new/NM_S30_7_MSPIP_test/updated_MRaabe_LW_091221_171221_Expl2_XL_Ecoli_NM_S30_bRP_rep1_7_perc_0.0100_XLs.idXML" \
        -protocol "NM_S30_7"

python correlation_analysis.py -in_peps "/home/ubuntu/mspip_all_new/NM_S30_8_MSPIP_test/updated_MRaabe_LW_091221_171221_Expl2_XL_Ecoli_NM_S30_bRP_rep1_8_perc_0.1000_peptides.idXML" \
        -in_XLs "/home/ubuntu/mspip_all_new/NM_S30_8_MSPIP_test/updated_MRaabe_LW_091221_171221_Expl2_XL_Ecoli_NM_S30_bRP_rep1_8_perc_0.0100_XLs.idXML" \
        -protocol "NM_S30_8"

python correlation_analysis.py -in_peps "/home/ubuntu/mspip_all_new/NM_S30_9_MSPIP_test/updated_MRaabe_LW_091221_171221_Expl2_XL_Ecoli_NM_S30_bRP_rep1_9_perc_0.1000_peptides.idXML" \
        -in_XLs "/home/ubuntu/mspip_all_new/NM_S30_9_MSPIP_test/updated_MRaabe_LW_091221_171221_Expl2_XL_Ecoli_NM_S30_bRP_rep1_9_perc_0.0100_XLs.idXML" \
        -protocol "NM_S30_9"

python correlation_analysis.py -in_peps "/home/ubuntu/mspip_all_new/NM_S100_6_MSPIP_test/updated_MRaabe_LW_091221_171221_Expl2_XL_Ecoli_NM_S100_bRP_rep1_6_perc_0.1000_peptides.idXML" \
        -in_XLs "/home/ubuntu/mspip_all_new/NM_S100_6_MSPIP_test/updated_MRaabe_LW_091221_171221_Expl2_XL_Ecoli_NM_S100_bRP_rep1_6_perc_0.0100_XLs.idXML" \
        -protocol "NM_S100_6"

echo " ================>>>>>>>>>>>>>>>>>> RNA_UV <<<<<<<<<<<<<<<============================"

python correlation_analysis.py -in_peps "/home/ubuntu/mspip_all_new/UV_S30_9_MSPIP_test/updated_M_Raabe_A_Wulf_220421_270421_Expl3_Ecoli_XL_UV_S30_LB_bRPfrac_9_perc_0.1000_peptides.idXML" \
        -in_XLs "/home/ubuntu/mspip_all_new/UV_S30_9_MSPIP_test/updated_M_Raabe_A_Wulf_220421_270421_Expl3_Ecoli_XL_UV_S30_LB_bRPfrac_9_perc_0.0100_XLs.idXML" \
        -protocol "UV_S30_9"

python correlation_analysis.py -in_peps "/home/ubuntu/mspip_all_new/UV_S100_9_MSPIP_test/updated_M_Raabe_A_Wulf_220421_260421_Expl3_Ecoli_XL_UV_S100_LB_bRPfrac_9_perc_0.1000_peptides.idXML" \
        -in_XLs "/home/ubuntu/mspip_all_new/UV_S100_9_MSPIP_test/updated_M_Raabe_A_Wulf_220421_260421_Expl3_Ecoli_XL_UV_S100_LB_bRPfrac_9_perc_0.0100_XLs.idXML" \
        -protocol "UV_S100_9"

python correlation_analysis.py -in_peps "/home/ubuntu/mspip_all_new/UV_S100_10_MSPIP_test/updated_M_Raabe_A_Wulf_220421_260421_Expl3_Ecoli_XL_UV_S100_LB_bRPfrac_10_perc_0.1000_peptides.idXML" \
        -in_XLs "/home/ubuntu/mspip_all_new/UV_S100_10_MSPIP_test/updated_M_Raabe_A_Wulf_220421_260421_Expl3_Ecoli_XL_UV_S100_LB_bRPfrac_10_perc_0.0100_XLs.idXML" \
        -protocol "UV_S100_10"

echo " ================>>>>>>>>>>>>>>>>>> RNA_DEB <<<<<<<<<<<<<<<============================"

python correlation_analysis.py -in_peps "/home/ubuntu/mspip_all_new/DEB_S30_8_MSPIP_test/updated_M_Raabe_A_Wulf_220421_290421_Expl3_Ecoli_XL_DEB_S30_LB_bRPfrac_8_perc_0.1000_peptides.idXML" \
        -in_XLs "/home/ubuntu/mspip_all_new/DEB_S30_8_MSPIP_test/updated_M_Raabe_A_Wulf_220421_290421_Expl3_Ecoli_XL_DEB_S30_LB_bRPfrac_8_perc_0.0100_XLs.idXML" \
        -protocol "DEB_S30_8"

python correlation_analysis.py -in_peps "/home/ubuntu/mspip_all_new/DEB_S30_9_MSPIP_test/updated_M_Raabe_A_Wulf_220421_290421_Expl3_Ecoli_XL_DEB_S30_LB_bRPfrac_9_perc_0.1000_peptides.idXML" \
        -in_XLs "/home/ubuntu/mspip_all_new/DEB_S30_9_MSPIP_test/updated_M_Raabe_A_Wulf_220421_290421_Expl3_Ecoli_XL_DEB_S30_LB_bRPfrac_9_perc_0.0100_XLs.idXML" \
        -protocol "DEB_S30_9"


echo " ================>>>>>>>>>>>>>>>>>> RNA_4SU <<<<<<<<<<<<<<<============================"

python correlation_analysis.py -in_peps "/home/ubuntu/mspip_all_new/_4SU_SAX_MSPIP_test/updated_RNA_UV_4SU_ECOLI_AChernev_070119_Ecoli_365_4thioU_SAX_mzML_corrected_perc_0.1000_peptides.idXML" \
        -in_XLs "/home/ubuntu/mspip_all_new/_4SU_SAX_MSPIP_test/updated_RNA_UV_4SU_ECOLI_AChernev_070119_Ecoli_365_4thioU_SAX_mzML_corrected_perc_0.0100_XLs.idXML" \
        -protocol "_4SU_SAX_"

python correlation_analysis.py -in_peps "/home/ubuntu/mspip_all_new/_4SU_Qisi_MSPIP_test/updated_RNA_UV_4SU_ECOLI_AChernev_030119_Ecoli_365-4thioU_QiSi-2_mzML_corrected_perc_0.1000_peptides.idXML" \
        -in_XLs "/home/ubuntu/mspip_all_new/_4SU_Qisi_MSPIP_test/updated_RNA_UV_4SU_ECOLI_AChernev_030119_Ecoli_365-4thioU_QiSi-2_mzML_corrected_perc_0.0100_XLs.idXML" \
        -protocol "_4SU_Qisi_"