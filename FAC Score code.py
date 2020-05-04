"""
FAC SCORE: This code takes a list of peak areas from a set of FAC 
samples run at the same time and calculates FAC scores for each ion within
each set of FAC replicates in order to identify heterologously expressed
metabolites associated with the FAC under study.

Date started: 4/16/2020
Date completed: 5/4/2020
Coded by: Lindsay Caesar
"""

"""
-------------------------Importing libraries and files------------------------

This module imports libraries needed for the code, as well as the .csv file 
containing peak areas for FAC scoring.

NOTE: make sure .csv file also contains a "group" category where biological
replicates belong to same group. Group and mz values should be columns, with
FAC names and peak areas as rows.

-------------------------------------------------------------------------------
"""

import pandas as pd
import numpy as np

#read .csv file
df = pd.read_csv ('file:///C:/Users/lkc6464/Desktop/test_FAC_forpythondev.csv')
#copy file name here. uncomment the following to view first 10 rows and make 
#sure your file imported correctly. If you see too many columns and "NaN" under
#the columns, you can delete the columns, save the file, and reimport

#df.head(10)


"""
------------------Peak area averages for each ion by group--------------------

This module calculates the average peak area of each ion per group of 
biological replicates, as well as the average peak area of all other samples.

------------------------------------------------------------------------------
"""

# groupby to calculate mean of groups and means of samples outside of group.
group = df.groupby(['group'])

# calculate mean peak area within each FAC
def calc_FAC_mean(group):
    for x in group:
        FAC_mean=df.groupby(['group']).mean()
    return FAC_mean

# calculates mean of all other samples outside of each FAC
def calc_reference_mean(group):
    for x in group:
       total_sum=df.sum(axis=0)
       group_sum=df.groupby(['group']).sum()
       not_group_sum= total_sum-group_sum
       total_samples = df.count(axis=0)
       group_samples = df.groupby(['group']).count()
       not_group_samples = total_samples-group_samples
       reference_mean = not_group_sum/not_group_samples
       reference_mean = reference_mean.drop(columns="group")
    return reference_mean

"""
-------------------Ratio of ion abundance in FAC versus others-----------------

This module takes the peak areas calculated in the previous modules and
calculates a ratio of each ion's abundance in the FAC versus other samples. If
the ion is only detected in the FAC and not in the other samples, it gets a 
fixed ratio value of 99999.

-------------------------------------------------------------------------------
"""
#assign variables
reference_mean=calc_reference_mean(group)
FAC_mean = calc_FAC_mean(group)


#calculates ratio of FAC ion abundance/rest of samples. If only in FAC, ratio
#set to 99999.
def calc_ratio(group):
    for x in group:
        ratio = FAC_mean/reference_mean
        ratio_filtered = ratio.replace([np.inf, -np.inf], 99999)
    return ratio_filtered

"""
----------------------Log transformations--------------------------------------

This module does a log10 transformation of previous calculations, including the
ion abundance average within groups and the ratio of the FAC abundance to all
other samples. This module also calculates the average and standard deviation
of the log transformed averages of all ions detected within a group. Log
transformed values of 0 are excluded from total group average and st deviation.

-------------------------------------------------------------------------------
"""

#assign variables
ratio = calc_ratio(group)


#calculates log of FAC_mean
def calc_log_FAC_mean(group):
    for x in group:
        log_FAC_mean = np.log10(FAC_mean)
        log_FAC_mean_filtered = log_FAC_mean.replace([np.inf, -np.inf], np.nan)
    return log_FAC_mean_filtered

#calculates log of ratio
def calc_log_ratio(group):
    for x in group:
        log_ratio = np.log10(ratio)
        log_ratio_filtered = log_ratio.replace([np.inf, -np.inf], np.nan)
    return log_ratio_filtered




"""
---------------------------------Z-score---------------------------------------

This module calculates the Z-score value, using the log transformed 
average of each ion within a group, the total average of all log transformed 
ions within a group, and the stdev of all log transformed ions within a group.

-------------------------------------------------------------------------------
"""
#assign more dataframes
log_FAC_mean = calc_log_FAC_mean(group)

#calculate Z score
def calc_Z_score(log_FAC_mean):
    Zscore = (log_FAC_mean - np.nanmean(log_FAC_mean))/np.nanstd(log_FAC_mean)
    return Zscore


"""
--------------------------FAC Scoring------------------------------------------

This module calculates the FAC scores by group. Non-normalized FAC scores, 
normalized FAC scores with a standard subtraction of 3.4 (based on original
publication), and normalized FAC scores using data specific to this dataset are 
included. User should decide which normalized score to use.

-------------------------------------------------------------------------------
"""

#define variables
Zscore = calc_Z_score(log_FAC_mean)
log_ratio = calc_log_ratio(group)
     
def are_both_negative(group):
    for x in group:
        set_to_truefalse = (Zscore < 0) & (log_ratio < 0)
    return set_to_truefalse
    
#truefalse dataframe
truefalse = are_both_negative(group)

#converts truefalse dataframe to binary matrix
def set_to_binary(group):
    for x in group:
        set_to_binary = truefalse.replace([True, False], [0,1])
    return set_to_binary

#variable for binary dataframe
binarymultiplier = set_to_binary(group)

#calculates FAC score--multiplies by binary multiplier value to set values to 0
def FAC_score(group):
    for x in group:
        FAC_score = Zscore * log_ratio * binarymultiplier
    return FAC_score

#non-normalized FAC score dataframe
FAC_score_df = FAC_score(group)

#adds max value column to FAC_score_df dataframe
FAC_score_df['max_value'] = FAC_score_df.max(axis=1)

#variable to hold dST max
dST_max = FAC_score_df.at['dST','max_value']

#subtracts dST max from FAC scores. RECOMMENDED normalization approach
def normalized_FAC_score_dST(FAC_score_df):
    for x in group:
        normalized_FAC_score = FAC_score_df - dST_max
    return normalized_FAC_score

#subtracts 3.4 from FAC scores.
def normalized_FAC_score_3pt4(FAC_score_df):
    for x in group:
        normalized_FAC_score_3pt4=FAC_score_df - 3.4
    return normalized_FAC_score_3pt4


"""
--------------------------Export to .csv---------------------------------------

This module exports the normalized FAC scores to a .csv file in the same path
where the peak list is stored. User should specify which normalized score
to export based on choice in previous section (or produce multiple .csv files
with different types of normalization for comparison)

-------------------------------------------------------------------------------
"""
#export normalized FAC scores to file, where columns + rows represent ions and
# groups, but instead of peak areas, we have FAC scores.

#this is where user will decide which FAC score file to export
FAC_score_for_export = normalized_FAC_score_dST(FAC_score_df)

FAC_score_for_export.to_csv('CladoniaFACs_python_test_2.csv')