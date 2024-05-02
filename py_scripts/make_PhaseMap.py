from itertools import chain, combinations
import scipy.io
import os
import matplotlib.pyplot as plt

import numpy as np

import csv
import glob
import random

import re

from ECGtools import *


#place for UQ files

output_dir_pm = "/Volumes/DATAS/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_reentry_hires_PhaseMaps/"
output_dir_df = "/Volumes/DATAS/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_reentry_hires_DF/"
output_dir_ht = "/Volumes/DATAS/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_reentry_hires_HTi/"
output_dir_pm = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_reentry_hires_PhaseMaps_sampled/"
output_dir_df = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_reentry_hires_DF_sampled/"
output_dir_ht = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_reentry_hires_HTi_sampled/"

#input_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_hires/"
input_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_reentry_hires_surface/"

#input_dir = "/Volumes/DATAS/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_reentry_hires_surface/"

input_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_reentry_hires_sampled/"

def get_input_data(input_file, varname=[]):
    
    tmp = scipy.io.loadmat(os.path.join(input_file))
    
    if len(varname)==0:
        vars=[k for k in tmp.keys() if not k[0]=='_']
        if len(vars)==0:
            raise ValueError("No data found in file",os.path.join(input_file))
        elif len(vars)>1:
            print("multiple variables found in file, using: ",vars[0])
        varname=vars[0]
    
    return np.transpose(tmp[varname])


def main():
    
  files = os.listdir(input_dir)
  files.sort()

  print(len(files))

  sf = 1000
#  num_runs = 1
  cnt = 0
  for f in files:
    if ".mat" in f and not f[0] == ".":
      print(f)
      TMPs = get_input_data(os.path.join(input_dir,f))
      print(TMPs.shape)
      
      HT, PM, cycles, mx_cycles = phaseAnalysis(TMPs.T, sf)
      
#      figs = plot_PM_check(TMPs.T, HT = HT, phase = PM, cycles = cycles, mx_cycles = mx_cycles, sample_rate = sf)
#      figs = plot_PM_check(TMPs, sampling_rate = sf)
      
      print(PM.shape)
      print(mx_cycles)
      cycles_pad = np.zeros((len(cycles), mx_cycles))
      for k in range(len(cycles)):
#                print(cycles[k][0])
#                print(type(cycles[k][0]))
#                print(len(cycles[k][0]))
        cycles_pad[k,0:len(cycles[k][0])] = cycles[k][0]
          
      scipy.io.savemat(os.path.join(output_dir_pm, f), dict( PM = PM, cycles = cycles_pad))
      
      scipy.io.savemat(os.path.join(output_dir_ht,f), dict(HTi = HT.imag))
                 
      
      DF_out = dominantFrequency(TMPs, sf)
      
      print(DF_out[0].shape)
      scipy.io.savemat(os.path.join(output_dir_df,f),dict(DF = DF_out[0]))
      
#      cnt+=1
      
#      if cnt>=num_runs:
#        break
          
    
#    print(TMPs[0].shape)
#    LATs = compute_LAT(TMPs)
    
  return

if __name__ == "__main__":
    main()
