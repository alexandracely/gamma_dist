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
    
    def __init__(self, lqd, gas, res, ave_lqd_mw):
        self.liquid = lqd
        self.gas = gas
        self.res = res
        self.av_lqd_mw = ave_lqd_mw
    
    def assign_composition(self):
        pass
    
    def split_heavy_end(self, start_SCN):
        pass

# Class that contains functionality to parse a Core Labs Excel report
# and pass the data to FlashExperimentData class instance.
# openpyxl only works with Excel 2010+ file format.
class CoreLabsXLSXLoader:
    
    # Depending whether a specific worksheet is input a class instance has
    # a functionality either to read the whole flash data found in the
    # workbook or just that individual worksheet.
    def __init__(self, report_path, ws=None):
        self.file_path = report_path
        self.worksheet = ws
    
    # Method to parse an individual worksheet.
    def read_flash_data(self, ws):
        worksheet = ws
        # Reading the composition table as the whole.
        df = pd.read_excel(self.file_path, worksheet, skiprows = 11, nrows = 52, 
                        usecols = 'B:I', header = None,
                        names = ['scn', 'cl_name', 'lqd_mp', 'lqd_wp', 'gas_mp',
                                 'gas_wp', 'res_mp', 'res_wp'], na_values = ' ')
        # Filling in empty cells in carbon group columns.
        df['scn'] = df['scn'].ffill()
        df.set_index(['scn', 'cl_name'], inplace=True)
        liq = df.iloc[:, 0:2]
        gas = df.iloc[:, 2:4]
        res = df.iloc[:, 4:]
        return liq, gas, res
    
    def read(self):
        wb = load_workbook(self.file_path)
        if self.worksheet is None:
            flash_data_list = [worksheet for worksheet in wb.sheetnames if
                               re.search('C\.\d+', worksheet)]
        else:
            flash_data_list = [self.worksheet]
        for worksheet in flash_data_list:
            sheet = wb[worksheet]
            desc = sheet['B8'].value+sheet['B9'].value
            depth, sample_num, cylinder = self.__parser(desc)
            print('Depth: ', depth, 'Sample number: ', sample_num, 'Cylinder: ', cylinder)
            liq, gas, res = self.read_flash_data(worksheet)
            # Typically, flashed liquid average mole weight in CL reports is in
            # the cell O39:
            lqd_av_mw = sheet['B8'].value+sheet['O39'].value
            
        print(liq)
    
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
    cl_report = CoreLabsXLSXLoader(path, ws='C.1')
    cl_report.read()
    # cl_report.read_flash_data()