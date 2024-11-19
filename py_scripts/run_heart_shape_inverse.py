from itertools import chain, combinations
import scipy.io
import os
import argparse

import numpy as np
import math
import scipy.signal as sig
from tabulate import tabulate

from UncertainSCI.distributions import BetaDistribution
from UncertainSCI.indexing import TotalDegreeSet
from UncertainSCI.pce import PolynomialChaosExpansion

import matplotlib.pyplot as plt

import csv
import glob
import random

import re
import sys


from ECGtools import *
from EPUQtools import *
from run_heart_shape import *

exec(open("initValues.py").read())

def build_inv_parser(parser):
  parser.add_argument('--forwardset', required=False,
                      help='forward data to use.  {[pseudoECG], NGpseudoECG, ECGsimBSP, ECG, recorded}',
                      default = "pseudoECG")
  parser.add_argument('--forwardtype', required=False,
                      help='forward data to use.  {[mean_torso], individual, mean_heart_geom}',
                      default = "mean_torso")
  return parser
  
  
def run_pipeline(sample_params, options, **kwargs):

  pipePaths = set_dirs(options)
  
  if options.compute_samples:
    if os.path.exists(samples_file):
      print('checking file')
      recognized = False
      while not recognized:
        print("Sample file exists.  Overwrite? [y/N]")
        str_input = str(input())
        if len(str_input)==0 or str_input.lower()=='n' or str_input.lower()=='no':
          print('exiting')
          return False
        elif str_input.lower()=='y' or str_input.lower()=='yes':
          recognized=True
    
    print("running samples")
            
    make_samples(samples_file, pce)
    print("Samples compute.  Make sure to run the model")
#        return
    samples = pce.samples
    options.run_model = True
  else:
    tmp = scipy.io.loadmat(samples_file)
    samples = tmp["samples"]
      
  if options.run_model:
    files = make_shape_points(samples)
    ind_func_files = make_indicator_functions(files)
    print("make cleaver meshes")
    
    return
      
      
      
  pce = set_distribution(sample_params)
  
  fid = fid_map[options.pacing]
  
  if not fid and options.hardchecks:
    raise ValueError(" FID missing and cannot be run with hardchecks on or in batch mode ")
  
  if not fid and not "LAT" in options.resultset:
    options.make_FIDs=True
  
  
  sol_input = load_solutions(resultset, pacing, options)
  
  if sol_input:
    model_output = sol_input[0]
    model_size = sol_input[1]
    ofiles = sol_input[2]
  else:
    print("cannot load "+pacing+" data.  skipping ")
    return
  
  if options.make_FIDs:
    choose_FIDS(model_output, model_size)
    return


  print("loaded samples")
  print(len(model_output), model_output[0].shape)
  print(samples.shape)
  
  if np.shape(samples)[0]>len(model_output) and not options.hardchecks:
    print("Paring down samples to match solutions")
    samples = correct_samples(samples, pce, ofiles)
    if len(samples)==0:
      return
    print(samples)
    print("new sample size",samples.shape)
  elif np.shape(samples)[0]<len(model_output) and not options.hardchecks:
    raise ValueError("more solutions than samples. TODO: implement pare down for this case")
#    elif  (not np.shape(samples)[0]==len(model_output)) and options.hardchecks:
#        raise ValueError("sample and model mismatch that cannot be resolved in batch mode or with hardchecks")

#    samples = sample_workaround(samples, domain)

  print(output_dir)
  UQ_file = os.path.join(output_dir, resultset+"_"+pacing+"_UQ_values.mat")
  print(UQ_file)
  RMS_UQ_file = os.path.join(output_dir, resultset+"_"+pacing+"_RMS_UQ_values.mat")
  
  Q = 4  # Number of quantile bands to plot

  # only two dimensions currently
  N = np.prod(model_size)
  print(N)
  print(model_size)
  do_RMS=not np.any(N==model_size)
  
  if not os.path.exists(UQ_file):
    if options.hardchecks:
      print(" cannot run pce implicitly with hardchecks on, or in batch mode.  Skipping "+pacing)
      return
    options.run_pce = True
      
  
  if options.run_pce or not os.path.exists(UQ_file):
      
    pce.set_samples(samples)
    
    print(np.vstack(model_output).shape)
    
    run_pce_data(pce, np.vstack(model_output),model_size, UQ_file, Q)
    print("PCE model built")
      
  if do_RMS:
    print("running RMS")
    RMS = model_RMS(model_output, model_size)
    
    
    print(RMS.shape)
#        print(RMS)
    if not os.path.exists(RMS_UQ_file):
      if options.hardchecks:
        print(" cannot run pce implicitly with hardchecks on, or in batch mode.  Skipping "+pacing)
        return
        
      options.run_pce = True
        
        
    if options.run_pce:
      print("running RMS PCE")
      pce_rms = set_distribution(sample_params)
      pce_rms.set_samples(samples)
      run_pce_data(pce_rms, RMS, (1,model_size[0]), RMS_UQ_file, Q)
      
      
  pot_UQ = load_pce_data(UQ_file)
  
  print(" pce stats loaded from disk")

  summary_stats = make_summary_stats(pot_UQ, len(model_output))
  print("pot stats")
  print(tabulate(summary_stats, headers = "firstrow"))
      
  if do_RMS:
    rms_UQ = load_pce_data(RMS_UQ_file)
    
    dt = 0.5

    lead_num = 404

    print(model_size)
    quantiles = np.reshape(pot_UQ["quantiles"], (pot_UQ["quantiles"].shape[0], model_size[0], model_size[1]))

    time = np.linspace(0,model_size[1]*dt,model_size[1])
    RMS_med = calc_rms(pot_UQ["median"],axis=0)
    RMS_quant = calc_rms(quantiles,axis=1)
    
    print("median time")
    print(pot_UQ["median"].shape)
    
    if not model_size == pot_UQ["median"].shape:
        pot_UQ["median"] = np.reshape(pot_UQ["median"], model_size)
    print(pot_UQ["median"].shape)
    print(time.shape)
    plot_ecg_uq(time, pot_UQ["median"], quantiles, model_output, model_size,0, title =  pacing)

#    if options.plot_channels>=0:
#        plot_all_channels(model_output, LATs, dys, pot_size, options.plot_channels)
     
#    if options.plot_solutions:
#        channels = (stdev_LAT>3*np.mean(stdev_LAT)).nonzero()[0]
#        plot_all_solutions(model_output, LATs, dys*50, pot_size, channels)
        
    print("rms_UQ")
    print(RMS.shape)
    print(rms_UQ["median"].shape)
    plot_quantile_bands(time, rms_UQ["median"], rms_UQ["quantiles"], Q, title =  pacing, ylabel = "RMS (mV)", xlabel = "time (ms)")
      
      
  return
    


def main():
    
  parser = build_parser()
  parser = build_inv_parser(parser)
  options = parser.parse_args()
  
  options.pipeline = "Inverse"
  
  print(type(options))
  print(options)

#    global resultset, pacing
      
  pacing = options.pacing
  resultset = options.resultset
  
  sample_params = default_sample_params

  sample_params["distribution"] = distrib
  
  set_dirs(options.data_dir, resultset, pacing)
  
  domain = get_ranges_from_scores(sample_params)
  
  # TODO: make optional rerun
  compute_surface_variance(domain)
  
  print("domain")
  print(domain)
  
  
  
  if len(domain) == 0:
    domain = default_domain
      
  sample_params["domain"] = domain
  

  
  if pacing == "all" or pacing == "new" or pacing == "old":
    options.compute_samples=False
    options.run_model=False
    options.make_FIDs=False
    options.hardchecks=True
    
    for p in pacing_lists[pacing]:
      options.pacing = p
      run_pipeline(resultset, p, sample_params, options)
      
  else:
    run_pipeline(resultset, pacing, sample_params, options)

  plt.show()
  
  return
    

if __name__ == "__main__":
  main()
