from itertools import chain, combinations
import scipy.io
import os

import numpy as np

from base64 import b64encode
import random
import argparse
import math
import scipy.signal as sig

from tabulate import tabulate
import re
import sys
import csv
import glob

import matplotlib.pyplot as plt


from UncertainSCI.distributions import BetaDistribution, NormalDistribution
from UncertainSCI.indexing import TotalDegreeSet
from UncertainSCI.pce import PolynomialChaosExpansion

def loadMatlabMatrix(filename, varname_def=[]):

  tmp = scipy.io.loadmat(filename)
  if not varname_def:
    vars=[k for k in tmp.keys() if not k[0]=='_']
    if not vars:
      print("No data found in file",os.path.join(results_dir,f))
      return False
    elif len(vars)>1:
      print("multiple variables found in file, using: ",vars[0])
    varname=vars[0]
  else:
    varname = varname_def
    
  return tmp[varname]
  
def check_dimension(geom):
  #geom must be dict with "pts" and "fac" fields

  if not geom["pts"].shape[0]==3:
    geom["pts"] = geom["pts"].T
    
    if not  geom["pts"].shape[0]==3:
      raise ValueError("points must be 3D")

  if not geom["fac"].shape[0]==3:
    geom["fac"] = geom["fac"].T
    
    if not  geom["fac"].shape[0]==3:
      raise ValueError("points must be 3D")

  return geom

def load_geom(filename):
    
#    print(filename)
    
  tmp=scipy.io.loadmat(filename)

  geom = {}

  #    print(tmp.keys())

  no_pts = False
  no_fac = False

  if "pts" in tmp.keys():
    geom["pts"] = tmp["pts"]
  elif "cor" in tmp.keys():
    geom["pts"] = tmp["cor"]
  else:
    no_pts = True
    
  if "fac" in tmp.keys():
    geom["fac"]=tmp["fac"]
  elif "face" in tmp.keys():
    geom["fac"]=tmp["face"]
  else:
    no_fac = True
    
  if no_fac and no_pts:
  #        print(type(tmp))
    kys = [ k for k in tmp.keys() if not  k[0]=="_"]
    field_found = False
    it = 0
    
    # couldn't get this to work.
    # TODO: generalize reader a bit
  #        print(k)
  #        while no_fac or no_pts and it<len(kys):
  #            print(type(tmp[kys[it]]))
  #            if type(tmp[kys[it]])==type(np.array([])):
  #                tmp_new = tmp[kys[it]](0,0)
  #            if "pts" in tmp_new.keys():
  #                geom["pts"] = tmp_new["pts"]
  #                no_pts = False
  #
  #            if "fac" in tmp_new.keys():
  #                geom["fac"] = tmp_new["fac"]
  #                no_fac = False
  #
  #            no_fac = False
  #            no_pts = False
  #            it +=1

    geom["pts"] = tmp[kys[0]]["pts"][0,0]
    geom["fac"] = tmp[kys[0]]["fac"][0,0]
    
    no_fac = False
    no_pts = False
    
  #    print(geom.keys())
  if no_fac or no_pts:
    raise ValueError(filename+" is not a dict with a pts and fac field")
    
  geom= check_dimension(geom)

  return geom


def calc_correlation(r_pots, gt_pots, axis = 1):

  axis = check_inputs(r_pots, gt_pots, axis)
  cc = np.sum(r_pots*gt_pots,axis = axis)/(np.linalg.norm(r_pots,axis = axis)*np.linalg.norm(gt_pots, axis = axis ))

  return cc
    

def calc_rmse(r_pots, gt_pots, axis = 1):

  axis = check_inputs(r_pots, gt_pots, axis)
  diff = gt_pots-r_pots
  rmse = calc_rms(diff,axis = axis)/np.linalg.norm(gt_pots, axis=axis)
  
  return rmse
  
def calc_rms_e(r_pots, gt_pots, axis=1):
    
  axis = check_inputs(r_pots, gt_pots, axis)
  diff = gt_pots-r_pots
  rms_e = calc_rms(diff,axis = axis)
  
  return rms_e
    
def calc_rms(pots, axis = 0 ):
        
  rms = np.sqrt(np.mean(pots*pots/pots.shape[axis],axis = axis))
  
  return rms

def calc_correlation_time(r_pots, gt_pots):
  return calc_correlation(r_pots, gt_pots, 0)
def calc_correlation_space(r_pots, gt_pots):
  return calc_correlation(r_pots, gt_pots, 1)
def calc_correlation_1D(r_pots, gt_pots):
  return calc_correlation(r_pots.flatten(), gt_pots.flatten())
def calc_rmse_time(r_pots, gt_pots):
  return calc_rmse(r_pots, gt_pots, 0)
def calc_rmse_space(r_pots, gt_pots):
  return calc_rmse(r_pots, gt_pots, 1)
def calc_rmse_1D(r_pots, gt_pots):
  return calc_rmse(r_pots.flatten(), gt_pots.flatten())
def calc_rms_e_time(r_pots, gt_pots):
  return calc_rms_e(r_pots, gt_pots, 0)
def calc_rms_e_space(r_pots, gt_pots):
  return calc_rms_e(r_pots, gt_pots, 1)
def calc_rms_e_1D(r_pots, gt_pots):
  return calc_rms_e(r_pots.flatten(), gt_pots.flatten())

def run_pce_metrics(indices, dist, samples, metrics, met_UQ_files,Q=4):

  for key, value in metrics.items():

    pce = PolynomialChaosExpansion(indices, dist)
    pce.set_samples(samples)
    run_pce_data(pce, value ,(value.size,1), met_UQ_files[key], Q)
    
  return
  
  
def run_pce_data(pce, model_output, model_size, sample_params, output_dir, UQ_file, Q = 4):

  print("running PCE")
  #    print(model_output.shape)
  pce.build(model_output=model_output)

  mean = pce.mean()
  stdev = pce.stdev()

  # Power set of [0, 1, ..., dimension-1]
  variable_interactions = list(chain.from_iterable(combinations(range(sample_params["dimension"]), r) for r in range(1, sample_params["dimension"]+1)))

  # "Total sensitivity" is a non-partitive relative sensitivity measure per parameter.
  total_sensitivity = pce.total_sensitivity()

  # "Global sensitivity" is a partitive relative sensitivity measure per set of parameters.
  global_sensitivity = pce.global_sensitivity(variable_interactions)

  #    global_sensitivity, variable_interactions = pce.global_sensitivity(interaction_orders=interaction_orders)

  labels = [' '.join([pce.plabels[v] for v in varlist]) for varlist in variable_interactions]

  dq = 0.5/(Q+1)
  q_lower = np.arange(dq, 0.5-1e-7, dq)[::-1]
  q_upper = np.arange(0.5 + dq, 1.0-1e-7, dq)
  quantile_levels = np.append(np.concatenate((q_lower, q_upper)), 0.5)

  quantiles = pce.quantile(quantile_levels, M=int(2e3))

  median = pce.quantile(0.5, M=int(1e3))

  med_peak = np.argmax(median)
  mean_peak = np.argmax(mean)
  qrs_peak = np.mean([med_peak, mean_peak])

  print("peak= ")
  print(qrs_peak)

  print("model_size")
  print(np.prod(model_size))
  print(model_size)
  print(np.any(np.prod(model_size)==model_size))
  if not np.any(np.prod(model_size)==model_size):
    print("resizing pce data")
    mean = mean.reshape(model_size)
    stdev = stdev.reshape(model_size)
    median = median.reshape(model_size)

  print(mean.shape)
  if not os.path.exists(output_dir+"UQ_data"):
    os.makedirs(output_dir+"UQ_data")
  print("saving :", UQ_file)
  scipy.io.savemat(UQ_file, dict(mean = mean.T, stdev = stdev.T, tot_sensitivity = total_sensitivity.T, glob_sensitivity = global_sensitivity.T, quantiles = quantiles.T, median = median.T, labels = labels, QRS_peak = np.array(qrs_peak)))

  return
  
def load_pce_data(filename):
  tmp = scipy.io.loadmat(filename)
  
#    print(tmp.keys())
  print(filename)
  UQ_stats = {}
  for key, value  in tmp.items():
      if key[:2] =="__":  continue
      print(key)
      UQ_stats[key] = value.T
      
  return UQ_stats
  
def load_pce_files(filenames_dict):
  UQ_stats_dict ={}
  
  for key, value in filenames_dict.items():
      UQ_stats_dict[key] = load_pce_data(value)
      
  return UQ_stats_dict

def check_metrics_files(metrics_files):
  if not type(metrics_files)==type({}):
    raise ValueError("input must be dict with metric name and file paths")
    
  check = [os.path.exists(v) for v in metrics_files.values()]

  return np.array(check).all()
        

    
def compute_metrics_set(r_pots, gt_pots, metrics_list):

  diff = gt_pots-r_pots

  switcher = {
            "correlation_t": calc_correlation_time,
            "correlation_s": calc_correlation_space,
            "correlation" : calc_correlation_1D,
            "rmse_t": calc_rmse_time,
            "rmse_s": calc_rmse,
            "rmse": calc_rmse_1D,
            "rms_t": calc_rms_e_time,
            "rms_s": calc_rms_e_space,
            "rms": calc_rms_e_1D,
  }
  metrics = {}

  for met in metrics_list:
    func = switcher.get(met, lambda: "Invalid method")
    metrics[met] = func(r_pots, gt_pots)

  return metrics


def save_metrics_files(metrics_files, metrics):
  if not len(metrics_files) == len(metrics):
    raise ValueError("metrics file list doesn't match metrics list")

  for key, value in metrics.items():
    scipy.io.savemat(metrics_files[key], {key : value})
    
  return True
    
def load_metrics_files(metrics_files):
  if not type(metrics_files)==type({}):
    raise ValueError("input must be dict with metric name and file paths")

  metrics = {}
  for key, value in metrics.items():
    tmp = scipy.io.savemat(value)
    metrics[key] = tmp[key]
    
  return metrics
  
def append_metrics(metrics, met_add):
#    print(type(met_add))
  metrics_new={}
  if len(metrics) ==0:
    metrics_new = met_add
  elif len(metrics) == len(met_add):
    for key, value in met_add.items():
        metrics_new[key] = np.vstack((metrics[key], value))
  else:
    raise ValueError("metrics dicts don't align")

  return metrics_new


def make_summary_stats_from_metrics(UQ_stats_metrics, num_samples, title=""):
    
  for key, value in UQ_stats_metrics.items():
  #        print(value.keys())
    summary_stats = make_summary_stats(value, num_samples)
    print(title+" "+key)
    print(tabulate(summary_stats, headers = "firstrow"))
    
  return  summary_stats
    
def make_summary_stats(UQ_stats, num_samples):
  std_error = UQ_stats["stdev"]/np.sqrt(num_samples)
  std_mean = UQ_stats["stdev"]/UQ_stats["mean"]
  summary_stats = [[ "stat","shape", "min", "max", "mean", "std"],
                  [ "mean", UQ_stats["mean"].shape, np.min(UQ_stats["mean"]), np.max(UQ_stats["mean"]), np.mean(UQ_stats["mean"]), np.std(UQ_stats["mean"])],
                  [ "median", UQ_stats["median"].shape, np.min(UQ_stats["median"]), np.max(UQ_stats["median"]), np.mean(UQ_stats["median"]), np.std(UQ_stats["median"])],
                  [ "stdev", UQ_stats["stdev"].shape, np.min(UQ_stats["stdev"]), np.max(UQ_stats["stdev"]), np.mean(UQ_stats["stdev"]), np.std(UQ_stats["stdev"])],
                  [ "std error", std_error.shape, np.min(std_error), np.max(std_error), np.mean(std_error), np.std(std_error)],
                  [ "stdev/mean", std_mean.shape, np.min(np.abs(std_mean)), np.max(np.abs(std_mean)), np.mean(np.abs(std_mean)), np.std(np.abs(std_mean))]
                  ]
  return summary_stats
  
def make_metrics_filenames(res_dir, fileroot, metrics_list):

  met_file_names = {}
  fileending = ".mat"

  for met in metrics_list:
    met_file_names[met] = os.path.join(res_dir, fileroot+met+fileending)
    
  return met_file_names
  
def make_UQ_filenames(output_dir, sample_params, options, tag = ""):

  fileending = ".mat"
  
  if len(tag)==0:
    tag = "UQ_values"

  UQ_filename = os.path.join(output_dir, "UQ_data", options.resultset+"_"+options.pacing+"_"+sample_params["distribution"]+"_"+tag+".mat")
    
  return UQ_filename
  
def check_Model_outputs(samples, model_output, pce, ofiles, options):

  if np.shape(samples)[0]>len(model_output):
  #        if options.hardchecks:
  #            raise ValueError("can't pair down samples with hardchecks on (or batch mode)")
    print("Paring down samples to match solutions")
    samples = correct_samples(samples, pce, ofiles)
    if len(samples)==0:
      return False
    print("new sample size",samples.shape)
    
  elif np.shape(samples)[0]<len(model_output):
    raise ValueError("more solutions than samples. TODO: implement pare down for this case")
    if options.hardchecks:
      return False
      
  return True


    
    
def model_RMS(model_output, pot_size, axis = 0):
  RMS = np.zeros((len(model_output),pot_size[axis-1]))
  for k in range(len(model_output)):
    m = np.resize(model_output[k], pot_size)
    RMS[k,:] = np.sqrt(np.sum(m*m/pot_size[axis],axis=axis))
  return RMS
  
      
    
def model_RMS(model_output, pot_size, axis = 0):
  RMS = np.zeros((len(model_output),pot_size[axis-1]))
  for k in range(len(model_output)):
    m = np.resize(model_output[k], pot_size)
    RMS[k,:] = np.sqrt(np.sum(m*m/pot_size[axis],axis=axis))
  return RMS
  
  
 
def plot_quantile_bands(time, median, quantiles, Q, filename, **kwargs):
      
  fig = plt.figure(figsize=(4,2))
  plt.plot(time, median.flatten(), 'b', label='RMS median')
  #    plt.spines[['right', 'top']].set_visible(False)

  print(quantiles.shape)
  band_mass = 1/(2*(Q+1))
  for ind in range(Q):
    alpha = (Q-ind) * 1/Q - (1/(2*Q))
    if ind == 0:
      plt.fill_between(time,quantiles[ind,:], quantiles[Q+ind, :],
                       interpolate=True, facecolor='red', alpha=alpha,
                       label='{0:1.2f} probability mass (each band)'.format(band_mass))
    else:
      plt.fill_between(time, quantiles[ind, :], quantiles[Q+ind, :], interpolate=True, facecolor='red', alpha=alpha)
  ax = fig.gca()
  ax.set(**kwargs)

  #    plt.subplots_adjust(left = 0.1, bottom=0.1, right=0.9, top=0.9)
  plt.tight_layout()

  if len(filename)>0:
    plt.savefig(filename)
  return
  
def save_stats(UQ_file, summary_stats):
  print("saving stats table")
  print(UQ_file[:-4]+"_table.csv")
  with open(UQ_file[:-4]+"_table.csv", 'w') as fid:
    write = csv.writer(fid)
    write.writerows(summary_stats)
  
  
if __name__ == "__main__":
  print("this is a library and I haven't made any unit tests yet")
