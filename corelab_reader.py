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
    
    def __init__(self, composition, ave_MW):
        pass
    
    def assign_composition(self):
        pass
    
    def split_heavy_end(self, start_SCN):
        pass

# Class that contains functionality to parse a Core Labs Excel report
# and pass the data to FlashExperimentData class instance.
# openpyxl only works with Excel 2010+ file format.
class CoreLabsXLSXLoader:
    
    def __init__(self, report_path):
        self.file_path = report_path
    
    def read_flash_data(self):
        worksheet = 'C.1'
        # Reading the composition table as the whole.
        df = pd.read_excel(self.file_path, worksheet, skiprows = 11, nrows = 52, 
                        usecols = 'B:I', header = None,
                        names = ['scn', 'cl_name', 'lqd_mp', 'lqd_wp', 'gas_mp',
                                 'gas_wp', 'res_mp', 'res_wp'], na_values = ' ')
        # Filling in empty cells in carbon group columns.
        df['scn'] = df['scn'].ffill()
        print(df)
    
    def read(self):
        wb = load_workbook(self.file_path)
        flash_data_list = [worksheet for worksheet in wb.sheetnames if
                           re.search('C\.\d+', worksheet)]
        print(flash_data_list)
        # pass
    
if __name__ == "__main__":
    
    path = '..\..\PVT_Reports\PS1.xlsx'
    cl_report = CoreLabsXLSXLoader(path)
    # cl_report.read()
    cl_report.read_flash_data()

# ws = wb['C.1']

# cell_range = ws['B12':'C63']

# print(cell_range)