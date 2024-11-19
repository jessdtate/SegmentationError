from itertools import chain, combinations
import scipy.io
import os
import argparse
import socket

import numpy as np
import math
#import scipy.signal as sig
#from scipy.fft import fft, fftfreq, ifft, fftshift
#from tabulate import tabulate

from UncertainSCI.distributions import BetaDistribution, NormalDistribution
from UncertainSCI.indexing import TotalDegreeSet
from UncertainSCI.pce import PolynomialChaosExpansion

import matplotlib.pyplot as plt

import csv
import glob
import random

import re
import sys

#hostname = socket.gethostname()
#print(hostname)
#
#if "metis" in hostname or "local" in hostname:
##if True:
#    sys.path.append("/Users/jess/CIBC/FP/UQ/py_scripts")
#    output_dir_root = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/"
#    shapedata_dir_root = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/"
#
#
#    torso_pots_dir = "/Users/jess/CIBC/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/Interventions/BSPM_pf/"
#
#
#    scirun_call = "/Users/jess/software/SCIRun_testing2/bin_515/SCIRun/SCIRun_test"
#    #scirun_call = "/Applications/SCIRun.app/Contents/MacOS/SCIRun"
#    seg3D_call = "/Applications/Seg3D2.app/Contents/MacOS/Seg3D2"
#    cleaver_call = "/Users/jess/software/cleaver2/build/bin/cleaver-cli"
#    
#elif "cibc" in hostname or "shell" in hostname:
#
#    sys.path.append("/home/sci/jess/UQ/py_scripts")
#    output_dir_root = "/home/sci/segmentation_error/Dalhousie_seg/UQ_data/Forward/"
#    shapedata_dir_root = "/home/sci/segmentation_error/Dalhousie_seg/shape_data/"
#
#
#    torso_pots_dir = "/usr/sci/projects/Edgar/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/Interventions/BSPM_pf/"


#print(sys.path)

from ECGtools import *
from EPUQtools import *
from run_heart_shape import *

#from run_torso_postition_for_dnn import run_pce_data, load_pce_data, run_pce_metrics, load_pce_files, make_summary_stats,plot_quantile_bands
#
#from check_DNN_results import load_geom, calc_correlation, calc_rmse, calc_rms, calc_rms_e

#calc_rms, calc_corr, calc_rmse, calc_rms_e, cal

#place for UQ files
exec(open("initValues.py").read())

def build_fwd_parser(parser):
  parser.add_argument('--resultset', required=False,
                      help='pacing profile to run.  {[pseudoECG], NGpseudoECG, inriaLAT, ECGsimLAT, ECGsimBSP, inriaTMP, ECG, all}',
                      default = "pseudoECG")
  return parser



def main():
    
  parser = build_parser()
  parser = build_fwd_parser(parser)
  options = parser.parse_args()
  
  options.pipeline = "Forward"
  
  print(type(options))
  print(options)
  
#  return
  #    global resultset, pacing
    
  pacing = options.pacing
  resultset = options.resultset
  distrib = options.distrib

  sample_params = default_sample_params

  sample_params["distribution"] = distrib

  pipePaths = set_dirs(options)

  sample_params = get_params_from_scores(sample_params, pipePaths)

  print(sample_params["std"])

  # TODO: make optional rerun
  compute_surface_variance(sample_params["std"], pipePaths)


  domain = np.vstack((sample_params["std"], -sample_params["std"]))

  print("domain")
  print(domain)



  if len(domain) == 0:
    domain = default_domain
    
    sample_params["covariance"] = default_cov
    sample_params["mean"] = default_mean
    
  sample_params["domain"] = domain



  if resultset == "all":
    options.compute_samples=False
    options.run_model=False
    options.make_FIDs=False
    options.hardchecks=True
    
    run_all_results(resultset_lists[resultset], sample_params, options, **pipePaths)
    
  else:
    run_all_results( [ resultset ], sample_params, options, **pipePaths)

  plt.show()

  return
    
if __name__ == "__main__":
  main()
