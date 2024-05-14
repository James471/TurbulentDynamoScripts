E_MAG_COLUMN_INDEX = 11
E_KIN_COLUMN_INDEX = 9
V_RMS_COLUMN_INDEX = 13
TIME_COLUMN_INDEX = 0
Cs_RMS_COLUMN_INDEX = 14


SOLVER_DICT = {"8wave": "Split-Roe", "bouchut-split": "Split-Bouchut", "Roe": "USM-Roe", 
              "HLLD": "USM-HLLD", "HLLC": "USM-HLLC", "bk-usm": "USM-BK"}

COLOR_DICT = {"8wave": "#377eb8", "HLLD": "#984ea3", "HLLC": "#4daf4a", 
              "Roe": "#a65628", "bouchut-split": "#f781bf", "bk-usm": "#ff7f00"}

ORDER_DICT = {"8wave": 0, "bouchut-split": 1, "Roe": 2, "HLLD": 3, "HLLC": 4, "bk-usm": 5}

FACT_DICT = {"bk-usm": 1.0, "Roe": 0.1, "HLLC": 0.01, "HLLD": 0.001, "bouchut-split": 0.0001, "8wave": 0.00001}