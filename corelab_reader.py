# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 19:26:05 2020

@author: Dmitry Molokhov (molokhov@outlook.com)
"""

from openpyxl import load_workbook
import pandas as pd
import re


# Class to store sample data and relevant fluid properties.
class FlashExperimentData:
    
    def __init__(self, lqd=None, gas=None, res=None, ave_lqd_mw=None):
        self.liquid = lqd
        self.gas = gas
        self.reservoir = res
        self.av_lqd_mw = ave_lqd_mw 
        self._c10_heavy_end_lqd = None
        self._c7_heavy_end_lqd = None
    
    def assign_composition(self):
        pass
    
    # Don't like this way of splitting the heavy end but I reckon
    # I exhausted my knowledge of Pandas.
    @property
    def c10_heavy_end_lqd(self):
        i = self.liquid[self.liquid['scn'] == 'C10'].index[0]
        self._c10_heavy_end_lqd = self.liquid.iloc[i:]
        return self._c10_heavy_end_lqd
    
    @property
    def c7_heavy_end_lqd(self):
        i = self.liquid[self.liquid['scn'] == 'C7'].index[0]
        self._c7_heavy_end_lqd = self.liquid.iloc[i:]
        return self._c7_heavy_end_lqd
    
    # The method to define the heavy end of the liquid fraction.
    # By default assumes that the fraction starts at C10.
    # def split_heavy_end_lqd(self, start_SCN = 'C10'):
    #     heavy_end = self.liquid.loc[start_SCN:]
    #     return heavy_end


class FlashExpDataCollection(FlashExperimentData, dict):
    
    def __init__(self, *worksheet_names):
        dict.__init__(self, ((wn, FlashExperimentData.__init__()) for wn in worksheet_names))
        
    @property
    def sample_names(self):
        return self.keys()
    
    def add_sample(self, sample_name, sample):
        self[sample_name].append(sample)

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
            flash_data_list = [self.worksheet]
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
            samples.add_sample(worksheet, FlashExperimentData(liq, gas, res, lqd_av_mw))
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
    
    path = '..\..\PVT_Reports\PS1.xlsx'
    cl_report = CoreLabsXLSXLoader(path, worksheet='C.1')
    sample_collection = cl_report.read()
    print(sample_collection.sample_names)
    
    # print(sample_collection[0].c10_heavy_end_lqd)
    # cl_report.read_flash_data()