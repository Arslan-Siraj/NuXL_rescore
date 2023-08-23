from idXML2df import readAndProcessIdXML
import pandas as pd
from scipy.stats import pearsonr
import numpy as np

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
     '-in_peps',
    type=str,
    metavar='in_peps'
)

parser.add_argument(
     '-in_XLs',
    type=str,
    metavar='in_XLs'
)

parser.add_argument(
     '-protocol',
    type=str,
    metavar='protocol'
)

args = parser.parse_args()

def convert_string_to_float_list(number_strings):
    result = []
    
    for num_str in number_strings:
        num_str = num_str.strip('[]')  # Remove brackets
        nums = num_str.split(', ')      # Split into individual number strings
        float_nums = [float(num) for num in nums]  # Convert to float
        
        result.append(float_nums)
    
    return result


def custom_function(row):
    ions_b_targ_xls_list = convert_string_to_float_list([row["ions_b_targ_xls"]])
    ions_b_targ_pep_list = convert_string_to_float_list([row["ions_b_targ_pep"]])
    first_three_same = all(value1 == value2 for value1, value2 in zip(ions_b_targ_xls_list[0][:3], ions_b_targ_pep_list[0][:3]))
    return not first_three_same

def custom_function_same(row):
    ions_b_targ_xls_list = convert_string_to_float_list([row["ions_b_targ_xls"]])
    ions_b_targ_pep_list = convert_string_to_float_list([row["ions_b_targ_pep"]])
    ions_y_targ_xls_list = convert_string_to_float_list([row["ions_y_targ_xls"]])
    ions_y_targ_pep_list = convert_string_to_float_list([row["ions_y_targ_pep"]])

    first_three_any_constant = (
        all(value == ions_b_targ_xls_list[0][0] for value in ions_b_targ_xls_list[0][:3]) or
        all(value == ions_b_targ_pep_list[0][0] for value in ions_b_targ_pep_list[0][:3]) or
        all(value == ions_y_targ_xls_list[0][0] for value in ions_y_targ_xls_list[0][:3]) or
        all(value == ions_y_targ_pep_list[0][0] for value in ions_y_targ_pep_list[0][:3])
    )
    
    return not first_three_any_constant

def custom_function_y(row):
    ions_b_targ_xls_list = convert_string_to_float_list([row["ions_y_targ_xls"]])
    ions_b_targ_pep_list = convert_string_to_float_list([row["ions_y_targ_pep"]])
    first_three_same = all(value1 == value2 for value1, value2 in zip(ions_b_targ_xls_list[0][:3], ions_b_targ_pep_list[0][:3]))
    return not first_three_same
    
#transform log to anti log 
def logantitransform(ls):
    ls_log_out = []
    for list in ls:
        list_temp = []
        for i in list_temp:
            list_temp.append((2 ** float(i)) - 0.001)
        ls_log_out.append(list_temp)
    ls_log_out

def first_three_elements_constant(lst):
    if len(lst) >= 3 and lst[0] == lst[1] == lst[2]:
        return True
    else:
        return False

def all_elements_constant(lst):
    if len(lst) < 2:
        return True
    
    first_element = lst[0]
    
    for element in lst:
        if element != first_element:
            return False
    
    return True

def correlation_analysis(pept_file_, XL_file_, protocol):
    peptides_df = readAndProcessIdXML(pept_file_)
    xls_df = readAndProcessIdXML(XL_file_)

    peptides_df['pep_mod'] = peptides_df["Peptide"] + '|' + peptides_df["NuXL:NA"]
    xls_df['pep_mod'] = xls_df["Peptide"] + '|' + xls_df["NuXL:NA"]
    print("Total_XLs: ", xls_df.shape)

    pep_ob_df = peptides_df[["Peptide", 'pep_mod', "ions_b_targ", "ions_y_targ"]]
    xls_ob_df = xls_df[["Peptide",  'pep_mod', "ions_b_targ", "ions_y_targ"]]

    #Get all same peptides
    merged_df_ = pd.merge(xls_ob_df, pep_ob_df, on="Peptide", suffixes=('_xls', '_pep'))#, how='left')
    print("AFter merge: ", merged_df_.shape)

    print("after merged: ", merged_df_.shape)

    condition_result = merged_df_.apply(custom_function, axis=1)
    merged_df_ = merged_df_[condition_result]
    print("after removing b same ions: ", merged_df_.shape)

    condition_result_y = merged_df_.apply(custom_function_y, axis=1)
    merged_df_ = merged_df_[condition_result_y]
    print("after removing same y ions: ", merged_df_.shape)

    '''condition_result_s = merged_df_.apply(custom_function_same, axis=1)
    merged_df_ = merged_df_[condition_result_s]
    print("after removing constant ions: ", merged_df_.shape)'''

    merged_df = merged_df_.drop_duplicates(subset=["pep_mod_xls"], keep='first')
    #print("after just one pep_mod_xls: ", merged_df.shape)

    peptides = merged_df["Peptide"]
    ions_b_targ_xls_ = convert_string_to_float_list(list(merged_df["ions_b_targ_xls"]))
    ions_y_targ_xls_ = convert_string_to_float_list(list(merged_df["ions_y_targ_xls"]))
    ions_b_targ_pep_ = convert_string_to_float_list(list(merged_df["ions_b_targ_pep"]))
    ions_y_targ_pep_ = convert_string_to_float_list(list(merged_df["ions_y_targ_pep"]))
    
    #intensity = (2 ** log2_intensity) - 0.001
    #ions_b_targ_xls_tran = [[(2 ** log2_intensity) - 0.001 for log2_intensity in sublist] for sublist in ions_b_targ_xls_]
    #merged_df["ions_b_targ_xls_tran"] = ions_b_targ_xls_tran

    #ions_b_targ_xls_log = logantitransform(ions_b_targ_xls_)
    #ions_y_targ_xls_log = logantitransform(ions_y_targ_xls_)
    #ions_b_targ_pep_log = logantitransform(ions_b_targ_pep_)
    #ions_y_targ_pep_log = logantitransform(ions_y_targ_pep_)

    #merged_df["ions_b_targ_xls_log"] = ions_b_targ_xls_log
    #merged_df["ions_y_targ_xls_log"] = ions_y_targ_xls_log
    #merged_df["ions_b_targ_pep_log"] = ions_b_targ_pep_log
    #merged_df["ions_y_targ_pep_log"] = ions_y_targ_pep_log

    ################## All Ions Correlation ####################


    correlation_values = []
    for i in range(len(ions_y_targ_xls_)):
        if all_elements_constant(ions_b_targ_xls_[i]) or all_elements_constant(ions_b_targ_pep_[i]):
            correlation_values.append(0.0)
        else:
            correlation, _ = pearsonr(np.array(ions_b_targ_xls_[i]), np.array(ions_b_targ_pep_[i]))
            correlation_values.append(correlation)
        
    merged_df["b_correlation(All)"] = correlation_values
    merged_df['protocol'] = protocol
    correlation_values = []
    for i in range(len(ions_y_targ_xls_)):
        if all_elements_constant(ions_y_targ_xls_[i]) or all_elements_constant(ions_y_targ_pep_[i]):
            correlation_values.append(0.0)
        else:
            correlation, _ = pearsonr(np.array(ions_y_targ_xls_[i]), np.array(ions_y_targ_pep_[i]))
            correlation_values.append(correlation)
        
    merged_df["y_correlation(All)"] = correlation_values

    merged_df["max_b_y(All)"] = merged_df[["b_correlation(All)", "y_correlation(All)"]].max(axis=1)
    print("##################### "+ protocol + "######################### \n")
    print("max (b_corr, y_corr) (All) :", merged_df["max_b_y(All)"].mean(), "+- ", merged_df["max_b_y(All)"].std())
    print("b_corr", merged_df["b_correlation(All)"].mean(), "+- ", merged_df["b_correlation(All)"].std())
    print("y_corr", merged_df["y_correlation(All)"].mean(), "+- ", merged_df["y_correlation(All)"].std()) 

    file_ = open(protocol + "_intensities_correlation.txt", "w") ##file to write the report
    file_.write("-----max_b_y(All) analysis---\n")

    file_.write("max (b_corr, y_corr) (All) : " + str(merged_df["max_b_y(All)"].mean())+ "+- " + str(merged_df["max_b_y(All)"].std()) + '\n')
    file_.write("b_corr: " + str(merged_df["b_correlation(All)"].mean()) +  "+- " + str(merged_df["b_correlation(All)"].std()) + '\n')
    file_.write("y_corr: " + str(merged_df["y_correlation(All)"].mean()) +  "+- " + str(merged_df["y_correlation(All)"].std()) + '\n')

    ########################## first three ions correlation ########################
    correlation_values = []
    for i in range(len(ions_y_targ_xls_)):
        x = ions_b_targ_xls_[i]
        y = ions_b_targ_pep_[i]

        if first_three_elements_constant(ions_b_targ_xls_[i]) or first_three_elements_constant(ions_b_targ_pep_[i]):
            correlation_values.append(0.0)
        else:
            correlation, _ = pearsonr(np.array(x[0:3]), np.array(y[0:3]))
            correlation_values.append(correlation)
        
    merged_df["b_correlation(3-prefix/suffix)"] = correlation_values

    correlation_values = []
    for i in range(len(ions_y_targ_xls_)):
        x = ions_y_targ_xls_[i]
        y = ions_y_targ_pep_[i]
        if first_three_elements_constant(ions_y_targ_xls_[i]) or first_three_elements_constant(ions_y_targ_pep_[i]):
            correlation_values.append(0.0)
        else:
            correlation, _ = pearsonr(np.array(x[0:3]), np.array(y[0:3]))
            correlation_values.append(correlation)
        
    merged_df["y_correlation(3-prefix/suffix)"] = correlation_values

    merged_df["max_b_y(3-prefix/suffix)"] = merged_df[["b_correlation(3-prefix/suffix)", "y_correlation(3-prefix/suffix)"]].max(axis=1)

    print("max (b_corr, y_corr) (3-prefix/suffix) :", merged_df["max_b_y(3-prefix/suffix)"].mean(), "+- ", merged_df["max_b_y(3-prefix/suffix)"].std())
    print("b_corr(3-prefix/suffix)", merged_df["b_correlation(3-prefix/suffix)"].mean(), "+- ", merged_df["b_correlation(3-prefix/suffix)"].std())
    print("y_corr(3-prefix/suffix)", merged_df["y_correlation(3-prefix/suffix)"].mean(), "+- ", merged_df["y_correlation(3-prefix/suffix)"].std()) 

    file_.write("-----max_b_y(3-prefix/suffix) analysis---\n")

    file_.write("max (b_corr, y_corr) (3-prefix/suffix) : " + str(merged_df["max_b_y(3-prefix/suffix)"].mean())+ "+- " + str(merged_df["max_b_y(3-prefix/suffix)"].std()) + '\n')
    file_.write("b_corr: " + str(merged_df["b_correlation(3-prefix/suffix)"].mean()) +  "+- " + str(merged_df["b_correlation(3-prefix/suffix)"].std()) + '\n')
    file_.write("y_corr: " + str(merged_df["y_correlation(3-prefix/suffix)"].mean()) +  "+- " + str(merged_df["y_correlation(3-prefix/suffix)"].std()) + '\n')

    file_.close()
    merged_df = merged_df.drop(columns=["ions_b_targ_xls", "ions_y_targ_xls", "ions_b_targ_pep", "ions_y_targ_pep"])
    merged_df.to_csv(protocol+'.csv')

if __name__ == "__main__":
    correlation_analysis(args.in_peps, args.in_XLs, args.protocol)