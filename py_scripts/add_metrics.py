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

output_dir_root = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/"
default_datadir = "vent_MD_cheat/11samples/"

results_dir = os.path.join(output_dir_root, default_datadir)

tag="NGpseudoECG"
tag="ECGsimBSP"


files = os.listdir(results_dir)
files.sort()

print(len(files))

for f in files:
    if tag in f and "RMS_UQ_values.mat" in f:
        print(f)
        RMS_data = scipy.io.loadmat(os.path.join(results_dir, f))
        
        f_data = f.replace("_RMS_", "_")
        
        UQ_data = scipy.io.loadmat(os.path.join(results_dir, f_data))
        

        med_peak = np.argmax(RMS_data["median"])
        print(med_peak)
        mean_peak = np.argmax(RMS_data["mean"])
        print(mean_peak)
        qrs_peak = int(np.mean([med_peak, mean_peak]))
        
        print(qrs_peak)
        
        RMS_data["QRS_peak"] = qrs_peak
        
        print(RMS_data.keys())
        
        UQ_data["QRS_peak"] = qrs_peak
        
        print(UQ_data.keys())
        
        f_txt = f_data.replace(".mat", "_visslice_index.txt")
        np.savetxt(os.path.join(results_dir, f_txt), np.array([qrs_peak]))
        
        
#        f_txt = f_data.replace(".mat", "_visslice_index.mat")
#        scipy.io.savemat(os.path.join(results_dir, f_txt), {"QRS_peak" : np.array([qrs_peak])})


        scipy.io.savemat(os.path.join(results_dir, f), RMS_data)
        scipy.io.savemat(os.path.join(results_dir, f_data), UQ_data)
