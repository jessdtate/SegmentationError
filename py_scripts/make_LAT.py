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
#output_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_hires_LAT/"
output_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_reentry_hires_LAT_sampled/"
#output_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_reentry_hires_LAT/"
#input_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_hires/"
#input_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_hires_new/"
input_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_reentry_hires_sampled/"
#input_dir = "/Volumes/DATAS/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_reentry_hires_surface/"

#def compute_LAT(TMPs,dt = 1.0 ):
#
#    diff = (TMPs[:,2:]-TMPs[:,0:-2])/2
#    LATs = np.argmax(diff,axis=1)*dt
#    return LATs
#
#def compute_LAT_from_TMP(TMP,dt = 1):
#
#
#    dv_dt = (TMPs[1:]-TMPs[0:-1])/dt
#    t_max = np.argmax(dv_dt)
#    LAT = t_max*dt
##    plt.plot(TMP)
##    plt.plot(dv_dt*np.max(dv_dt)/np.max(TMP))
##    plt.show()
#    input()
#    return LAT

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

  for f in files:
    if ".mat" in f and not f[0] == ".":
#  f = "Solution_model_params_172_apex_stim_point_sampled.mat"
      print(f)
      TMPs = get_input_data(os.path.join(input_dir,f)).T
      print(TMPs.shape)
      LATs, dvdt, mx_p = findLAT(TMPs, sample_rate = 1000, mode = "max",
                    n_LAT = 0, window = 50, ignore_first = 0)
                    
      print("LATs done")
      mx_lat = 0
      print(type(LATs))
      print(len(LATs))
      print(len(LATs[0]))
      print(type(LATs[0]))

      LAT_pad = np.zeros((len(LATs), mx_p))

      for k in range(len(LATs)):
        LAT_pad[k,0:len(LATs[k])] = np.array(LATs[k])
        



      #            LATs = np.array(LATs)
      #            print(LATs.shape)
      scipy.io.savemat(os.path.join(output_dir,f),dict(LATs = LAT_pad))
        



  #    print(TMPs[0].shape)
  #    LATs = compute_LAT(TMPs)

  return

if __name__ == "__main__":
    main()
