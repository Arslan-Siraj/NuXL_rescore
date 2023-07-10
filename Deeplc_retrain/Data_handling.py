#####################################################
###### Modification feature generation ##############
import pandas as pd
import numpy as np
    
def find_indices_of_substring(full_string, sub_string):
    return [index for index in range(len(full_string)) if full_string.startswith(sub_string, index)]


def alphadeep_col(Data_df, mod_seq_col = 'opt_global_cv_MS:1000889_peptidoform_sequence', NA_mod_col='opt_global_NuXL:NA', NA = 'RNA', mod1_ = 'Oxidation@M', mod2_ = 'Carbamidomethyl@C'):
    sequence = list(Data_df[mod_seq_col])
    mod_sites_ls = []
    seq_mod_ls = []

    mod1 = mod1_[:mod1_.find('@')]
    mod2 = mod2_[:mod2_.find('@')]

    #Carbamidomethyl@C
    for seq in sequence:
        mod1_indices = mod2_indices = []
        mod1_indices = find_indices_of_substring(seq, '(' + mod1 + ')')
        mod2_indices = find_indices_of_substring(seq, '(' + mod2 + ')')
        
        seq_mod = ''
        mod_sites = ''
        if (len(mod1_indices)>0 or len(mod2_indices)>0):
            if(len(mod1_indices)>0 and len(mod2_indices)==0):
                    check_ox = 0
                    for i in mod1_indices:
                        if (check_ox==0):
                            seq_mod = mod1_ 
                            mod_sites = str(i)
                            check_ox = check_ox + 1
                        else:
                            #print(i, check_ox, (check_ox)*(len(mod1)+2))
                            x = int(i) - int(check_ox)*(len(mod1)+2)
                            mod_sites = mod_sites + ';' + str(x)
                            seq_mod = seq_mod + ';' + mod1_
                            check_ox = check_ox + 1 
                            
            if(len(mod1_indices)==0 and len(mod2_indices)>0):
                    check_car = 0
                    for i in mod2_indices:
                        if (check_car==0):
                            seq_mod = mod2_
                            mod_sites = str(i)
                            check_car = check_car + 1
                        else:
                            #print(i, check_car, (check_car)*(len(mod2)+2))
                            x = int(i) - int(check_car)*(len(mod2)+2)
                            mod_sites = mod_sites + ';' + str(x)
                            seq_mod = seq_mod + ';' + mod2_
                            check_car = check_car + 1
                    #print(seq_mod, mod_sites)

            check_mod1 = 0
            check_mod2 = 0
            seq_mod_ls_ = []
            mod_sites_ls_ = []
            if(len(mod1_indices)>0 and len(mod2_indices)>0):
                #print(seq, len(seq),'Oxi: ', mod1_indices, 'Carbo: ', mod2_indices)

                mod1_len = int(len(mod1)+2)
                mod2_len = int(len(mod2)+2)
                while((len(mod1_indices)!=0) or (len(mod2_indices)!=0)):

                    #print(check_mod1,check_mod2, mod1_len*(check_mod1), mod2_len*(check_mod2))
                    maximium = len(mod1_indices)+len(mod2_indices)

                    for i in range(maximium):
                        if(len(mod1_indices)!=0 and len(mod2_indices)!=0):
                            mod1_min = min(mod1_indices)
                            mod2_min = min(mod2_indices)
                            
                            #print('mod1_min: ', mod1_min, 'mod2_min: ', mod2_min)
                            if(mod1_min < mod2_min):
                                #print('case1: ')
                                mod1_indices.pop(0)
                                mod_sites_ls_.append((mod1_min - mod1_len*(check_mod1) - mod2_len*(check_mod2)))
                                check_mod1 = check_mod1 + 1
                                seq_mod_ls_.append(mod1_)
                            else:
                                #print('case2: ')
                                mod2_indices.pop(0)
                                mod_sites_ls_.append((mod2_min - mod1_len*(check_mod1) - mod2_len*(check_mod2)))
                                check_mod2 = check_mod2 + 1
                                seq_mod_ls_.append(mod2_)
                    

                    if(len(mod1_indices)!=0 and len(mod2_indices)==0):
                        #print('case3: ')
                        for i in mod1_indices:
                            mod_sites_ls_.append((int(i) - mod1_len*(check_mod1) - mod2_len*(check_mod2)))
                            check_mod1 = check_mod1 + 1
                            seq_mod_ls_.append(mod1_)
                    
                    if(len(mod1_indices)==0 and len(mod2_indices)!=0):
                        #print('case4: ')
                        for i in mod2_indices:
                            mod_sites_ls_.append((int(i) - mod1_len*(check_mod1) - mod2_len*(check_mod2)))
                            check_mod2 = check_mod2 + 1
                            seq_mod_ls_.append(mod2_)

                    mod1_indices.clear()
                    mod2_indices.clear()

                #print(mod_sites_ls, seq_mod_ls)

                for i in range(len(mod_sites_ls_)):
                        if(i==0):
                            mod_sites = str(mod_sites_ls_[i])
                            seq_mod = str(seq_mod_ls_[i])
                        else:
                            mod_sites = mod_sites + ';' + str(mod_sites_ls_[i])
                            seq_mod = seq_mod + ';' + str(seq_mod_ls_[i])

        '''if(len(mod_sites)!=0):
            print(seq, mod_sites, seq_mod)'''
        mod_sites_ls.append(mod_sites)            
        seq_mod_ls.append(seq_mod)
    
    Data_df.rename(columns={NA_mod_col:'NA_modification'}, inplace=True)

    NA_ = list(Data_df['NA_modification'])
    NA_mod = []
    for i in NA_:
        if (i!= 'none'):
            NA_mod.append(str(i)+'@' + NA + '_XL')
        else:
            NA_mod.append('')

    Data_df['mods'] = seq_mod_ls
    Data_df['mod_sites'] = mod_sites_ls
    Data_df['NA_mod'] = NA_mod
    Data_df['nAA'] = Data_df['sequence'].apply(len)
    
    return Data_df

def Dataprepare_XL(Data_df, mod_seq_col = 'opt_global_cv_MS:1000889_peptidoform_sequence', NA_mod_col='opt_global_NuXL:NA', NA = 'RNA', mod1_ = 'Oxidation@M', mod2_ = 'Carbamidomethyl@C'):
    sequence = list(Data_df[mod_seq_col])
    mod_sites_ls = []
    seq_mod_ls = []

    mod1 = mod1_[:mod1_.find('@')]
    mod2 = mod2_[:mod2_.find('@')]

    #Carbamidomethyl@C
    for seq in sequence:
        mod1_indices = mod2_indices = []
        mod1_indices = find_indices_of_substring(seq, '(' + mod1 + ')')
        mod2_indices = find_indices_of_substring(seq, '(' + mod2 + ')')
        
        seq_mod = ''
        mod_sites = ''
        if (len(mod1_indices)>0 or len(mod2_indices)>0):
            if(len(mod1_indices)>0 and len(mod2_indices)==0):
                    check_ox = 0
                    for i in mod1_indices:
                        if (check_ox==0):
                            seq_mod = mod1_ 
                            mod_sites = str(i)
                            check_ox = check_ox + 1
                        else:
                            #print(i, check_ox, (check_ox)*(len(mod1)+2))
                            x = int(i) - int(check_ox)*(len(mod1)+2)
                            mod_sites = mod_sites + ';' + str(x)
                            seq_mod = seq_mod + ';' + mod1_
                            check_ox = check_ox + 1 
                            
            if(len(mod1_indices)==0 and len(mod2_indices)>0):
                    check_car = 0
                    for i in mod2_indices:
                        if (check_car==0):
                            seq_mod = mod2_
                            mod_sites = str(i)
                            check_car = check_car + 1
                        else:
                            #print(i, check_car, (check_car)*(len(mod2)+2))
                            x = int(i) - int(check_car)*(len(mod2)+2)
                            mod_sites = mod_sites + ';' + str(x)
                            seq_mod = seq_mod + ';' + mod2_
                            check_car = check_car + 1
                    #print(seq_mod, mod_sites)

            check_mod1 = 0
            check_mod2 = 0
            seq_mod_ls_ = []
            mod_sites_ls_ = []
            if(len(mod1_indices)>0 and len(mod2_indices)>0):
                #print(seq, len(seq),'Oxi: ', mod1_indices, 'Carbo: ', mod2_indices)

                mod1_len = int(len(mod1)+2)
                mod2_len = int(len(mod2)+2)
                while((len(mod1_indices)!=0) or (len(mod2_indices)!=0)):

                    #print(check_mod1,check_mod2, mod1_len*(check_mod1), mod2_len*(check_mod2))
                    maximium = len(mod1_indices)+len(mod2_indices)

                    for i in range(maximium):
                        if(len(mod1_indices)!=0 and len(mod2_indices)!=0):
                            mod1_min = min(mod1_indices)
                            mod2_min = min(mod2_indices)
                            
                            #print('mod1_min: ', mod1_min, 'mod2_min: ', mod2_min)
                            if(mod1_min < mod2_min):
                                #print('case1: ')
                                mod1_indices.pop(0)
                                mod_sites_ls_.append((mod1_min - mod1_len*(check_mod1) - mod2_len*(check_mod2)))
                                check_mod1 = check_mod1 + 1
                                seq_mod_ls_.append(mod1_)
                            else:
                                #print('case2: ')
                                mod2_indices.pop(0)
                                mod_sites_ls_.append((mod2_min - mod1_len*(check_mod1) - mod2_len*(check_mod2)))
                                check_mod2 = check_mod2 + 1
                                seq_mod_ls_.append(mod2_)
                    

                    if(len(mod1_indices)!=0 and len(mod2_indices)==0):
                        #print('case3: ')
                        for i in mod1_indices:
                            mod_sites_ls_.append((int(i) - mod1_len*(check_mod1) - mod2_len*(check_mod2)))
                            check_mod1 = check_mod1 + 1
                            seq_mod_ls_.append(mod1_)
                    
                    if(len(mod1_indices)==0 and len(mod2_indices)!=0):
                        #print('case4: ')
                        for i in mod2_indices:
                            mod_sites_ls_.append((int(i) - mod1_len*(check_mod1) - mod2_len*(check_mod2)))
                            check_mod2 = check_mod2 + 1
                            seq_mod_ls_.append(mod2_)

                    mod1_indices.clear()
                    mod2_indices.clear()

                #print(mod_sites_ls, seq_mod_ls)

                for i in range(len(mod_sites_ls_)):
                        if(i==0):
                            mod_sites = str(mod_sites_ls_[i])
                            seq_mod = str(seq_mod_ls_[i])
                        else:
                            mod_sites = mod_sites + ';' + str(mod_sites_ls_[i])
                            seq_mod = seq_mod + ';' + str(seq_mod_ls_[i])

        '''if(len(mod_sites)!=0):
            print(seq, mod_sites, seq_mod)'''
        mod_sites_ls.append(mod_sites)            
        seq_mod_ls.append(seq_mod)
    
    Data_df.rename(columns={NA_mod_col:'NA_modification'}, inplace=True)

    NA_ = list(Data_df['NA_modification'])
    NA_mod = []
    for i in NA_:
        if (i!= 'none'):
            NA_mod.append(str(i)+'@' + NA + '_XL')
        else:
            NA_mod.append('')

    modifications = []
    #mod_ = ""
    for i in range(len(seq_mod_ls)):
        mod_f = ""
        x = seq_mod_ls[i].split(';')
        y = mod_sites_ls[i].split(';')
        #print(x,len(x),y,len(y))
        for j in range(len(x)): 
            if(j==0):
                if (x[j] == ''):
                    break
                    print("jf", x[j])
                else:
                    if x[j] == mod1_ :  
                        mod_f = y[j]+ "|" + mod1
                    else: 
                        mod_f = y[j] + "|" + mod2
            elif(j>0):
                if x[j] == mod1_ :  
                    mod_f = mod_f +"|" + y[j]+ "|" + mod1
                else: 
                    mod_f = mod_f +"|" + y[j] + "|" + mod2

        #print("mod_f: ", mod_f)
        modifications.append(mod_f)

    sequence_simple = list(Data_df["seq"])
    #Data_df['nAA'] = Data_df['sequence'].apply(len)
    
    print("XL_file...")
    loc_position=list(Data_df["localization_position"])
    loc_NA_Mod = []
    for j in range(len(loc_position)):
        if(loc_position[j]!=-1):
                NA_mod_f = str(loc_position[j]+1)+ "|" + str(NA_mod[j])
                loc_NA_Mod.append(NA_mod_f)
        else:
            pos_ = int(len(sequence_simple[j])/2)
            NA_mod_f = str(pos_)+ "|" + str(NA_mod[j])
            #print(sequence_simple[j], pos_ , loc_position[j])
            loc_NA_Mod.append(NA_mod_f)

    Data_df['NA_mod_loc'] = loc_NA_Mod
    #NA_modifications=list(Data_df['opt_global_NuXL:NA']) #can deeplc column and make to RNA/DNA
    #for i in range(len(NA_modifications)):
    full_modifications=[]
    Alpha_mod_ls = []
    Alpha_mod_sites_ls = []
    for i in range(len(modifications)):
        mod__ = []
        fixed_mod_col = modifications[i].split('|')
        NA_modifications = loc_NA_Mod[i].split('|')
        check_NA = True
        check_only_NA = True
        #print(i, "TOgather: ", fixed_mod_col, NA_modifications)
        for y in range(0,len(fixed_mod_col)-1,2):
            check_only_NA = False
            if (int(fixed_mod_col[y]) == int(NA_modifications[0])):
                #print("   case_1: ==", fixed_mod_col[y], NA_modifications[0])
                mod__.append(fixed_mod_col[y])
                mod__.append(fixed_mod_col[y+1])
                if check_NA:
                    mod__.append(int(NA_modifications[0])+1)
                    mod__.append(NA_modifications[1])
                    check_NA = False

            elif(int(fixed_mod_col[y])<int(NA_modifications[0])):
                #print("   case_2: ==", len(fixed_mod_col), fixed_mod_col[y], NA_modifications[0])
                mod__.append(fixed_mod_col[y])
                mod__.append(fixed_mod_col[y+1])
                if (check_NA):
                    if(len(fixed_mod_col)==y+2):
                            mod__.append(int(NA_modifications[0]))
                            mod__.append(NA_modifications[1])
                            #print("case_2                          Final_arrangment.....",mod__)
                            check_NA = False
                    else:
                        if (fixed_mod_col[y+2] != NA_modifications[0]):
                            mod__.append(int(NA_modifications[0]))
                            mod__.append(NA_modifications[1])
                            #print("case_2                          Final_arrangment.....",mod__)
                            check_NA = False
            
            elif(int(fixed_mod_col[y]) > int(NA_modifications[0])):
                #print("   case_3: ==", "fixed_mod_col[y]", fixed_mod_col[y], "NA_modifications[0] ", NA_modifications[0])
                if check_NA:
                        mod__.append(int(NA_modifications[0]))
                        mod__.append(NA_modifications[1])
                        #print("case_2                          Final_arrangment.....",mod__)
                        check_NA = False
                mod__.append(fixed_mod_col[y])
                mod__.append(fixed_mod_col[y+1])
            #print("                          Final_arrangment.....",mod__)

        if (check_only_NA):
            mod__.append(int(NA_modifications[0]))
            mod__.append(NA_modifications[1])

        
        ####Can handle here the alpha peptdeep column
        #Creating_DeeplC modifications
        mod__str = ""
        for y in range(0,len(mod__)-1,2):
            #print(y)
            if(y==0):
                mod__str = str(mod__[y])+ '|' + str(mod__[y+1])
                #y=y+2
            else:
                mod__str = mod__str + '|' + str(mod__[y])+ '|' + str(mod__[y+1])
                #y=y+2
            
        #print(mod__str, mod__)
        full_modifications.append(mod__str)

        #Creating Alphapeptdeepcol
        Alpha_mod = ""
        Alpha_mod_sites = ""
        for y in range(0,len(mod__)-1,2):
            if(y==0):
                Alpha_mod_sites = str(mod__[y])
                if mod__[y+1] == mod1:
                    Alpha_mod = str(mod1_)
                elif mod__[y+1] == mod2:
                    Alpha_mod = str(mod2_)
                else:
                    Alpha_mod = str(mod__[y+1])
            else:
                Alpha_mod_sites = Alpha_mod_sites + ';' + str(mod__[y])

                if mod__[y+1] == mod1:
                    Alpha_mod = Alpha_mod +';' +str(mod1_)
                elif mod__[y+1] == mod2:
                    Alpha_mod = Alpha_mod +';' + str(mod2_)
                else:
                    Alpha_mod = Alpha_mod +';' + str(mod__[y+1])

        Alpha_mod_ls.append(Alpha_mod)
        Alpha_mod_sites_ls.append(Alpha_mod_sites)

    #AlphapepdeepCOl
    Data_df['mods'] = Alpha_mod_ls
    Data_df['mod_sites'] = Alpha_mod_sites_ls
    #Data_df['NA_mod'] = NA_mod
    Data_df['nAA'] = Data_df['seq'].apply(len)
    #Data_df['modifications'] = modifications

    #DeepLC_column
    Data_df['modifications'] = full_modifications

    return Data_df

def Dataprepare_nonXL(Data_df, mod_seq_col = 'opt_global_cv_MS:1000889_peptidoform_sequence', NA_mod_col='opt_global_NuXL:NA', NA = 'RNA', mod1_ = 'Oxidation@M', mod2_ = 'Carbamidomethyl@C'):
    sequence = list(Data_df[mod_seq_col])
    mod_sites_ls = []
    seq_mod_ls = []

    mod1 = mod1_[:mod1_.find('@')]
    mod2 = mod2_[:mod2_.find('@')]

    #Carbamidomethyl@C
    for seq in sequence:
        mod1_indices = mod2_indices = []
        mod1_indices = find_indices_of_substring(seq, '(' + mod1 + ')')
        mod2_indices = find_indices_of_substring(seq, '(' + mod2 + ')')
        
        seq_mod = ''
        mod_sites = ''
        if (len(mod1_indices)>0 or len(mod2_indices)>0):
            if(len(mod1_indices)>0 and len(mod2_indices)==0):
                    check_ox = 0
                    for i in mod1_indices:
                        if (check_ox==0):
                            seq_mod = mod1_ 
                            mod_sites = str(i)
                            check_ox = check_ox + 1
                        else:
                            #print(i, check_ox, (check_ox)*(len(mod1)+2))
                            x = int(i) - int(check_ox)*(len(mod1)+2)
                            mod_sites = mod_sites + ';' + str(x)
                            seq_mod = seq_mod + ';' + mod1_
                            check_ox = check_ox + 1 
                            
            if(len(mod1_indices)==0 and len(mod2_indices)>0):
                    check_car = 0
                    for i in mod2_indices:
                        if (check_car==0):
                            seq_mod = mod2_
                            mod_sites = str(i)
                            check_car = check_car + 1
                        else:
                            #print(i, check_car, (check_car)*(len(mod2)+2))
                            x = int(i) - int(check_car)*(len(mod2)+2)
                            mod_sites = mod_sites + ';' + str(x)
                            seq_mod = seq_mod + ';' + mod2_
                            check_car = check_car + 1
                    #print(seq_mod, mod_sites)

            check_mod1 = 0
            check_mod2 = 0
            seq_mod_ls_ = []
            mod_sites_ls_ = []
            if(len(mod1_indices)>0 and len(mod2_indices)>0):
                #print(seq, len(seq),'Oxi: ', mod1_indices, 'Carbo: ', mod2_indices)

                mod1_len = int(len(mod1)+2)
                mod2_len = int(len(mod2)+2)
                while((len(mod1_indices)!=0) or (len(mod2_indices)!=0)):

                    #print(check_mod1,check_mod2, mod1_len*(check_mod1), mod2_len*(check_mod2))
                    maximium = len(mod1_indices)+len(mod2_indices)

                    for i in range(maximium):
                        if(len(mod1_indices)!=0 and len(mod2_indices)!=0):
                            mod1_min = min(mod1_indices)
                            mod2_min = min(mod2_indices)
                            
                            #print('mod1_min: ', mod1_min, 'mod2_min: ', mod2_min)
                            if(mod1_min < mod2_min):
                                #print('case1: ')
                                mod1_indices.pop(0)
                                mod_sites_ls_.append((mod1_min - mod1_len*(check_mod1) - mod2_len*(check_mod2)))
                                check_mod1 = check_mod1 + 1
                                seq_mod_ls_.append(mod1_)
                            else:
                                #print('case2: ')
                                mod2_indices.pop(0)
                                mod_sites_ls_.append((mod2_min - mod1_len*(check_mod1) - mod2_len*(check_mod2)))
                                check_mod2 = check_mod2 + 1
                                seq_mod_ls_.append(mod2_)
                    

                    if(len(mod1_indices)!=0 and len(mod2_indices)==0):
                        #print('case3: ')
                        for i in mod1_indices:
                            mod_sites_ls_.append((int(i) - mod1_len*(check_mod1) - mod2_len*(check_mod2)))
                            check_mod1 = check_mod1 + 1
                            seq_mod_ls_.append(mod1_)
                    
                    if(len(mod1_indices)==0 and len(mod2_indices)!=0):
                        #print('case4: ')
                        for i in mod2_indices:
                            mod_sites_ls_.append((int(i) - mod1_len*(check_mod1) - mod2_len*(check_mod2)))
                            check_mod2 = check_mod2 + 1
                            seq_mod_ls_.append(mod2_)

                    mod1_indices.clear()
                    mod2_indices.clear()

                #print(mod_sites_ls, seq_mod_ls)

                for i in range(len(mod_sites_ls_)):
                        if(i==0):
                            mod_sites = str(mod_sites_ls_[i])
                            seq_mod = str(seq_mod_ls_[i])
                        else:
                            mod_sites = mod_sites + ';' + str(mod_sites_ls_[i])
                            seq_mod = seq_mod + ';' + str(seq_mod_ls_[i])

        '''if(len(mod_sites)!=0):
            print(seq, mod_sites, seq_mod)'''
        mod_sites_ls.append(mod_sites)            
        seq_mod_ls.append(seq_mod)
    
    Data_df.rename(columns={NA_mod_col:'NA_modification'}, inplace=True)

    NA_ = list(Data_df['NA_modification'])
    NA_mod = []
    for i in NA_:
        if (i!= 'none'):
            NA_mod.append(str(i)+'@' + NA + '_XL')
        else:
            NA_mod.append('')

    modifications = []
    #mod_ = ""
    for i in range(len(seq_mod_ls)):
        mod_f = ""
        x = seq_mod_ls[i].split(';')
        y = mod_sites_ls[i].split(';')
        #print(x,len(x),y,len(y))
        for j in range(len(x)): 
            if(j==0):
                if (x[j] == ''):
                    break
                    print("jf", x[j])
                else:
                    if x[j] == mod1_ :  
                        mod_f = y[j]+ "|" + mod1
                    else: 
                        mod_f = y[j] + "|" + mod2
            elif(j>0):
                if x[j] == mod1_ :  
                    mod_f = mod_f +"|" + y[j]+ "|" + mod1
                else: 
                    mod_f = mod_f +"|" + y[j] + "|" + mod2

        #print("mod_f: ", mod_f)
        modifications.append(mod_f)

    sequence_simple = list(Data_df["seq"])
    #Data_df['nAA'] = Data_df['sequence'].apply(len)
    
    fake_mod_NA = [ -1 for _ in range(len(sequence_simple))]
    Data_df['NA_mod_loc'] = fake_mod_NA

    
    #AlphapepdeepCOl
    Data_df['mods'] = seq_mod_ls
    Data_df['mod_sites'] = mod_sites_ls
    Data_df['nAA'] = Data_df['seq'].apply(len)

    #DeepLC_column
    Data_df['modifications'] = modifications
    return Data_df

