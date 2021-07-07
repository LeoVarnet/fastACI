function [ACI,results,cfg_ACI] = Script4_getACI_calculate(cfg_ACI, y, y_correct, X, U)

ACI = [];
results = [];
cfg_ACI = [];

error('The script %s is outdated, the new name for this file is fastACI_getACI (since 7 July 2021).',mfilename);