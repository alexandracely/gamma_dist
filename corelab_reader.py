# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 19:26:05 2020

@author: Dmitry Molokhov (molokhov@outlook.com)
"""

from openpyxl import load_workbook
import pandas as pd
import numpy as np
import re
import scipy.special as sps
import scipy.stats as stats
import scipy.optimize as optim
import time

# Class to store sample data and relevant fluid properties.
class FlashExperimentData:
    
    def __init__(self, lqd=None, gas=None, res=None, ave_lqd_mw=None):
        self.liquid = lqd
        self.gas = gas
        self.reservoir = res
        self.av_lqd_mw = ave_lqd_mw
        # Estimation of average MW of C7 & C10 plus fractions
        # can be implemented later.
        self._ave_C10_mw = None
        self._ave_C7_mw = None
        self._c10_heavy_end_lqd = None
        self._c7_heavy_end_lqd = None
        self.gamma_input = None
        self.gamma_output = None
    
    def assign_composition(self):
        pass
    
    # The following property decorators are to define the heavy 
    # end of the liquid fraction.
    # Don't like this way of splitting the heavy end but I reckon
    # I exhausted my knowledge of Pandas.
    @property
    def c10_heavy_end_lqd(self):
        i = self.liquid[self.liquid['scn'] == 'C10'].index[0]
        self._c10_heavy_end_lqd = self.liquid.iloc[i:]
        return self._c10_heavy_end_lqd

    @property
    def ave_C10_mw(self):
        self._ave_C10_mw = 220
        return self._ave_C10_mw
    
    @property
    def c7_heavy_end_lqd(self):
        i = self.liquid[self.liquid['scn'] == 'C7'].index[0]
        self._c7_heavy_end_lqd = self.liquid.iloc[i:]
        return self._c7_heavy_end_lqd
    
    def prepare_input(self, id, n=10):
        lqd=self.c10_heavy_end_lqd
        mw=self.av_lqd_mw
        self.gamma_input = lqd.copy()
        self.gamma_input['MWi_lab'] = self.gamma_input.apply(lambda x : x['lqd_wp']*mw/x['lqd_mp'], axis = 1)
        # Calculating initial values of upper bounds of individual MWn slices.
        # Upper bounds are considered as midway between SCN MW values.
        # C36+ upper bound is set to an arbitrary number (10000).
        self.gamma_input['ubound_init'] = self.gamma_input['MWi_lab']+(self.gamma_input['MWi_lab'].shift(-1)-self.gamma_input['MWi_lab'])/2
        self.gamma_input['ubound_init'] = self.gamma_input['ubound_init'].fillna(100000)
        
        # Generating regresion variables for component molecular weight bounds:
        self.gamma_input['ubound'] = 'm'+self.gamma_input['scn']
        
        # Adding top row to represent a lower boundary of C10 (or upper C9 boundary)
        top_index = 'C'+str(n-1)
        df_top = pd.DataFrame(pd.DataFrame([[top_index]+[np.nan] * (len(self.gamma_input.columns)-1)], columns=self.gamma_input.columns))
        self.gamma_input = df_top.append(self.gamma_input, ignore_index=True)
        
        # Calculating rescaled C10+ weight fractions.
        self.gamma_input['wni_lab'] = self.gamma_input['lqd_wp']/self.gamma_input['lqd_wp'].sum()
        
        # Adding a regression variable for lower C10 molecular weight boundary:
        self.gamma_input.at[0,'ubound'] = 'ita'
        self.gamma_input.iloc[-1, self.gamma_input.columns.get_loc('ubound')] = self.gamma_input.iloc[-1, self.gamma_input.columns.get_loc('ubound_init')]
        self.gamma_input['sample_id'] = id.replace('.', '_')
        self.gamma_input['heavy_end_mw'] = id.replace('.', '_')+'_heavy_mw'


# Class to store multiple sample data.
class FlashExpDataCollection(FlashExperimentData, dict):
    
    def __init__(self, worksheet_names):
        dict.__init__(self, ((wn, FlashExperimentData.__init__(self)) for wn in worksheet_names))
        
    @property
    def sample_names(self):
        return list(self.keys())
    
    def add_sample(self, sample_name, sample):
        self[sample_name] = sample
        
    def _prepare_regression(self):
        for key, item in self.items():
            item.prepare_input(key)
            
    def gamma_distribution(self, reg_vals, reg_vars, rmse_switch = False):

        # Ensuring consistency of input data
        assert len(reg_vals) == len(reg_vars)
        
        # Creating a dictionary of regression values indexed with variable names:
        lookup = dict()
        for i in range(len(reg_vals)):
            lookup[reg_vars[i]] = reg_vals[i]
        
        # Updating the dataframe with set regression values (replacing variables with values):
        # Note, this time it is not a single dataframe but a class containing a number of dataframes.
        error_array = pd.Series(dtype='float64')
        for key, item in self.items():
            df = item.gamma_input
            df = df.replace(lookup)
            ind = key.replace('.', '_')+'_heavy_mw'
            # Below equation references are from the SPE Phase Behavior monograph.
            beta = (lookup[ind]-lookup['ita'])/lookup['alpha'] # Equation 5.14
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
            
            temp = (df.loc[df.index[1:-1], 'Wni']-df.loc[df.index[1:-1], 'wni_lab'])**2
            error_array = pd.concat([error_array, temp], ignore_index=True)
            
            if rmse_switch:
                df['Zni'] = df['Wni']/df['Mi']*df['Wi'].sum(skipna = True)
                item.gamma_output = df

        # rmse = 100*error_array.mean()**0.5
    
        return 100*error_array.mean()**0.5
        
    def gamma_distribution_fit(self, n=10, alpha=1):
        # Gamma distribution fit from C7 is not implemented yet.
        if n == 7:
            raise NotImplementedError
        self._prepare_regression()
        # reg_variables = np.concatenate(( 
        #                                 self[sample_names[0]].gamma_input.loc[self[sample_names[0]].gamma_input.index[0:-1], 'ubound'].unique())#,
                                        # self[sample_names[0]].gamma_input.loc[self._fit_df.index[0:-1], 'heavy_end_mw'].unique()))
        # reg_variables = np.delete(reg_variables, np.where(reg_variables == 100000))
        df = self[self.sample_names[0]].gamma_input
        reg_variables = np.concatenate((np.array(['alpha']), df.loc[df.index[0:-1], 'ubound'].unique(),
                                        np.array([sample.replace('.', '_')+'_heavy_mw' for sample in self.sample_names])))
        init_vals = df.loc[df.index[0:-1], 'ubound_init']
        init_vals.iloc[0] = init_vals.iloc[1]-14
        init_vals = pd.Series(alpha).append(init_vals, ignore_index=True)
        for key, item in self.items():
            init_vals = init_vals.append(pd.Series(item.ave_C10_mw), ignore_index=True)
            
        ub = init_vals+init_vals*0.02
        ub[17:] = init_vals[17:]+init_vals[17:]*0.05
        lb = init_vals-init_vals*0.02
        lb[17:] = init_vals[17:]-init_vals[17:]*0.05
        lb[0] = -np.inf
        ub[0] = np.inf
        
        # self.gamma_distribution(init_vals, reg_variables)
        res = optim.minimize(self.gamma_distribution, args=(reg_variables), x0=init_vals,
                             method = 'SLSQP', bounds=optim.Bounds(lb, ub), options={'maxiter':10000})
        res_df = pd.DataFrame({'Variables': reg_variables, 'Values':res.x})
        print('RMSE: ', res.fun)

        self.gamma_distribution(res.x, reg_variables, rmse_switch = True)
        
        for key, item in self.items():
            item.gamma_output.to_csv(r'.\DATA\out_new.csv', mode='a') # [['SCN', 'Mi', 'Wni', 'Zni']]
        
        # lookup = dict()
        # for i in range(len(init_vals)):
        #     lookup[reg_variables[i]] = init_vals[i]
        
        # for key, item in self.items():
        #     item.gamma_input = item.gamma_input.replace(lookup)
        #     item.gamma_input.to_csv('out.csv', mode='a')# = item.gamma_input.replace(lookup)
        
        # print(len(reg_variables))
        # print(len(init_vals))
        # print(len(reg_variables) == len(init_vals))
        # assert len(reg_variables) == len(init_vals)


# Class that contains functionality to parse a Core Labs Excel report
# and pass the data to FlashExperimentData class instance.
# openpyxl only works with Excel 2010+ file format.
class CoreLabsXLSXLoader:
    
    # Depending whether a specific worksheet is input a class instance has
    # a functionality either to read the whole flash data found in the
    # workbook or just that individual worksheet.
    def __init__(self, report_path, worksheet=None):
        self.file_path = report_path
        self.worksheet = worksheet
    
    # Method to parse an individual worksheet.
    def read_flash_data(self, worksheet):
        # worksheet = ws
        # Reading the composition table as the whole.
        df = pd.read_excel(self.file_path, worksheet, skiprows = 11, nrows = 52, 
                        usecols = 'B:I', header = None,
                        names = ['scn', 'cl_name', 'lqd_mp', 'lqd_wp', 'gas_mp',
                                 'gas_wp', 'res_mp', 'res_wp'], na_values = ' ')
        # Filling in empty cells in carbon group columns.
        df['scn'] = df['scn'].ffill()
        df.set_index(['scn', 'cl_name'], inplace=True)
        # df.sort_index(inplace=True)
        liq = df.iloc[:, 0:2].reset_index()
        gas = df.iloc[:, 2:4].reset_index()
        res = df.iloc[:, 4:].reset_index()
        return liq, gas, res
    
    def read(self):
        wb = load_workbook(self.file_path)
        if self.worksheet:
            flash_data_list = self.worksheet
        else:
            flash_data_list = [worksheet for worksheet in wb.sheetnames if
                               re.search('C\.\d+', worksheet)]
        samples = FlashExpDataCollection(flash_data_list)
        for worksheet in flash_data_list:
            sheet = wb[worksheet]
            desc = sheet['B8'].value+sheet['B9'].value
            depth, sample_num, cylinder = self.__parser(desc)
            # print('Depth: ', depth, 'Sample number: ', sample_num, 'Cylinder: ', cylinder)
            liq, gas, res = self.read_flash_data(worksheet)
            # Typically, flashed liquid average mole weight in CL reports is in
            # the cell O39:
            lqd_av_mw = sheet['O39'].value
            samples.add_sample(worksheet, FlashExperimentData(liq, gas, res,
                                                              lqd_av_mw))
        return samples
    
    # Parcer function is supposed to extract useful sample descriptors from
    # the text strings above the composition table. In most Core LAbs reports
    # these strings are in the cells B8 and B9. Presumably, they contain
    # sample depth, cylinder and sample numbers. At this stage I am only
    # extracting numerical part of the these numbers. Depth is assumed to be
    # MD depth in metres.
    def __parser(self, descript_string):
        try:
            depth = float(re.search('Depth\D+(\d+\.?\d*)', descript_string).group(1))
        except:
            depth = None
        try:
            sample_num = re.search('Sample N\D+(\d+\.?\d*)', descript_string).group(1)
        except:
            sample_num = 'N/A'
        try:
            cylinder = re.search('(Cylinder|Chamber)\D+(\d+)', descript_string).group(2)
        except:
            cylinder = 'N/A'
        return depth, sample_num, cylinder

if __name__ == "__main__":
    
    start_time = time.time()
    input_file = '.\DATA\PS1.xlsx'
    cl_report = CoreLabsXLSXLoader(input_file) #, worksheet=['C.1', 'C.4']
    sample_collection = cl_report.read()
    # print(sample_collection['C.1'].c10_heavy_end_lqd)
    sample_collection.gamma_distribution_fit()
    
    print("--- Execution time %s seconds ---" % (time.time() - start_time))
    
    # print(sample_collection[0].c10_heavy_end_lqd)
    # cl_report.read_flash_data()