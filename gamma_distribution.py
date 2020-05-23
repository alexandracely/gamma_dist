# -*- coding: utf-8 -*-
"""
Created on Fri May 22 10:09:02 2020

@author: Dmitry Molokhov (molokhov@outlook.com)
"""
# READ IT FIRST:
# Gamma distribution fitting assumes C10 as starting component.
# The algoritm is based on Curtis Whitson's work derived primarily
# from Phase Behavior SPE Monograph (Volume 20) by Whitson and Brule.
# Additional information source are numerous papers on Gamma distribution
# which can be found at Whitson's website https://whitson.com/publications/
# Particularly useful was the 2019 paper by Bilal Younus, Curtis Whitson et al
# "Field-Wide Equation of State Model Development" that is also available from
# the download section of their website.
# The script take a .csv file with the following columns with names (in square brakets):
# 1st column [SCN]: SCN identifiers (e.g., C10, C11, C12 etc.);
# 2nd column [mfi_lab]: mole fraction of component as per full composition;
# 3rd column [wfi_lab]: weight fraction of component as per full composition.
# Average sample molecular weight is entered in the main section of the code.

import pandas as pd
import numpy as np
import scipy.optimize as optim
import scipy.special as sps
import scipy.stats as stats
pd.set_option('display.max_columns', 500)

# Function to prepare the data for distribution fitting.
# Within the function we backcalculate component MWs used by that lab (though they can be entered
# directly); normalise molecular weights and finally get our initial estimate of component molecular
# weight upper bounds.
def prepare_input(df, mw):
    
    # Back-calculating MWs
    df['MWi'] = df.apply(lambda x : x['wfi_lab']*mw/x['mfi_lab'], axis = 1)
    
    # Calculating initial values of upper bounds of individual MWn slices.
    # Upper bounds are considered as midway between SCN MW values.
    # C36+ upper bound is set to an arbitrary number (10000).
    df['ubound_init'] = df['MWi']+(df['MWi'].shift(-1)-df['MWi'])/2
    df['ubound_init'] = df['ubound_init'].fillna(100000)
    
    # Generating regresion variables for component molecular weight bounds:
    df['ubound'] = 'm'+df['SCN']
    # df['alpha'] = 'alpha'
    # df['ave_mC10plus'] = 'ave_mC10plus'
    
    # Adding top row to represent a lower boundary of C10 (or upper C9 boundary)
    df_top = pd.DataFrame(pd.DataFrame([['C9']+[np.nan] * (len(df.columns)-1)], columns=df.columns))
    df = df_top.append(df, ignore_index=True)
    
    # Calculating rescaled C10+ weight fractions.
    df['wni_lab'] = df['wfi_lab']/df['wfi_lab'].sum()
    
    # Adding a regression variable for lower C10 molecular weight boundary:
    df.at[0,'ubound'] = 'ita'
    df.at[27,'ubound'] = df.at[27,'ubound_init']
    
    return df

def gamma_distribution(reg_vals, reg_vars, df):
    
    # Ensuring consistency of input data
    assert len(reg_vals) == len(reg_vars)
    
    # Creating a dictionary of regression values indexed with variable names:
    lookup = dict()
    for i in range(len(reg_vals)):
        lookup[reg_vars[i]] = reg_vals[i]
    
    # Updating the dataframe with set regression values (replacing variables with values):
    df = df.replace(lookup)
    
    # Below equation references are from the SPE Phase Behavior monograph.
    beta = (lookup['ave_mC10plus']-lookup['ita'])/lookup['alpha'] # Equation 5.14
    df['y'] = (df['ubound']-lookup['ita'])/beta # Equation 5.22
    # Below is the equation 5.20 from the monograph but with correction of the typo.
    # Correct form can be derived from the Equation 5.13.
    # Another typo is in the Equation 5.15 of the original SPE monograph. Correct
    # form of the equation can be found in 1990 Whiton's paper "Application of the Gamma
    # Distribution Model to MW and Boiling Point Data For Petroleum Fractions", Equation 21
    df['Q'] = (np.exp(-df['y'])*(df['y']**lookup['alpha'])/
               sps.gamma(lookup['alpha']))
    df['P0'] = stats.gamma.cdf((df['ubound']-lookup['ita']), 
                               a=lookup['alpha'], scale=beta) # Equation 5.18
    df['P1'] = df['P0']-(df['Q']/lookup['alpha']) # Equation 5.19
    df['Mi'] = (lookup['ita']+lookup['alpha']*beta*
                ((df['P1']-df['P1'].shift())/(df['P0']-df['P0'].shift()))) # Equation 5.17
    df['Wi'] = df['Mi'] * (df['P0']-df['P0'].shift()) # Weight
    df['Wni'] = df['Wi']/df['Wi'].sum(skipna = True) # Normalised weight fraction
    # Finally calculating RMSE between lab and calculated data. Converting it to 
    # percentage as it is a bigger number and better for the solver
    rmse = 100*((df.loc[df.index[0:-1], 'Wni']-df.loc[df.index[0:-1], 'wni_lab'])**2).mean()**.5

    return rmse

if __name__ == "__main__":
    
    # Reading .csv file with SCN identifiers, mole and weight fractions as three input columns
    # Column names are assumed to be in the top row and  they are SCN, mfi_lab, wfi_lab
    # Note, fractions are not normalised.
    comp_input = pd.read_csv(r'.\DATA\gamma_dist_input.csv', header = 0, index_col = False)
    
    # Whole sample MW.
    sample_mw = 161.60
    # C10+ molecular weight if available in lab report. Otherwise a reasonable estimate.
    ave_MC10plus = 225.0
    # Initial value for C10 lower bound (C9 upper bound) or ita as per Whitson's monograph.
    # Assumed somewhere in between C9 amd C10 molecular weight. Can be calculated with the
    # Equation 5.15 from the monograph (correcting the typo) or just by subtracting 14
    # from the upper bound.
    ita = 131.0
    
    # Preparing input for the regression:
    df = prepare_input(comp_input, sample_mw)
    
    # Forming an array of regression variables (excluding MC36 upper bound):
    reg_variables = np.concatenate((df.loc[df.index[0:-1], 'ubound'].unique(), np.array(['alpha', 'ave_mC10plus'])))
    
    # Forming an array of initial values corresponding to these variables (excluding MC36 upper bound):
    init_vals = np.concatenate((df.loc[df.index[0:-1], 'ubound_init'], np.array([1.0, ave_MC10plus])))
    init_vals[np.isnan(init_vals)] = ita
    
    res = optim.minimize(gamma_distribution, args=(reg_variables, df), x0=init_vals, 
                            method = 'Nelder-Mead', options={'maxiter':10000})
    res_df = pd.DataFrame({'Variables': reg_variables, 'Values':res.x})

    print(res_df)
    print('RMSE:', res.fun)
    # test.to_csv(r'.\DATA\out_test.csv')