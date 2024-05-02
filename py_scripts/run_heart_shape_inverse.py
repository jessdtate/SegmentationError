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

sys.path.append("/Users/jess/CIBC/FP/UQ/py_scripts")

print(sys.path)

from run_torso_postition_for_dnn import run_pce_data, load_pce_data, run_pce_metrics, load_pce_files, make_summary_stats,plot_quantile_bands

from check_DNN_results import load_geom, calc_correlation, calc_rmse, calc_rms, calc_rms_e

#calc_rms, calc_corr, calc_rmse, calc_rms_e, cal

#place for UQ files
output_dir_root = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Inverse/"
shapedata_dir_root = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/"
default_datadir = "vent_MD_cheat/11samples/"
Forward_dir_root = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/"


torso_pots_dir = "/Users/jess/CIBC/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/Interventions/BSPM_pf/"


scirun_call = "/Users/jess/software/SCIRun_testing2/bin_515/SCIRun/SCIRun_test"
#scirun_call = "/Applications/SCIRun.app/Contents/MacOS/SCIRun"




fid_map = { "sinus" : {"sf": 0.5, "qon" : 0, "qoff" : 37, "stoff" : 100, "toff": 200},
            "apex" : {"sf": 0.5, "qon" : 0, "qoff" : 54, "stoff" : 110, "toff": 220},
            "LV" : {"sf": 0.5, "qon" : 0, "qoff" : 65, "stoff" : 114, "toff": 230},
            "RV" : {"sf": 0.5, "qon" : 0, "qoff" : 66, "stoff" : 116, "toff": 230},
            "septal" : {"sf": 0.5, "qon" : 0, "qoff" : 50, "stoff" : 108, "toff": 220},
            "RVV" : {"sf": 0.5, "qon" : 0, "qoff" : 60, "stoff" : 110, "toff": 225},
            "LVV" : {"sf": 0.5, "qon" : 0, "qoff" : 60, "stoff" : 115, "toff": 220},
            "Rvot" : {"sf": 0.5, "qon" : 0, "qoff" : 60, "stoff" : 105, "toff": 220},
            "Lvot" : {"sf": 0.5, "qon" : 0, "qoff" : 60, "stoff" : 105, "toff": 220},
            "RVB" : {"sf": 0.5, "qon" : 2, "qoff" : 65, "stoff" : 105, "toff": 225}
}

pacing_lists = {"all" : [ k for  k in fid_map.keys()],
                "old" : ["sinus", "apex", "LV", "RV", "septal"],
                "new" : ["RVV", "LVV", "Rvot", "Lvot", "RVB"]
}


#domain = np.array([[-165, 165], [-112, 112], [-85, 85], [-55, 55], [-45,45]]).T
# 1.5 sigma
default_domain = np.array([[-101, 101], [-55, 55], [-40, 40], [-35, 35], [-25,25]]).T
# this is uniform, but probably should be guassian from ShapeWorks
sample_params = { "dimension": 5, "alpha": 1, "beta": 1}
# TODO: should probably be running these as gaussian



def build_parser():
    parser = argparse.ArgumentParser()

    # This will be implemented as rollout broadens
    parser.add_argument('--data_dir', required=False,
                        help='directory where the shape data is located, relative to the shape_data dir',
                        default = default_datadir)
    parser.add_argument('--pacing', required=False,
                        help='pacing profile to run.  {[sinus], apex, RV, LV, septal, RVV, LVV, Rvot, Lvot, RVB, all, old, or new}',
                        default = "sinus")
    parser.add_argument('--forwardset', required=False,
                        help='forward data to use.  {[pseudoECG], NGpseudoECG, ECGsimBSP, ECG, recorded}',
                        default = "pseudoECG")
    parser.add_argument('--forwardtype', required=False,
                        help='forward data to use.  {[mean_torso], individual, mean_heart_geom}',
                        default = "mean_torso")
    parser.add_argument('--run_model',
                        help='run the model with the samples.  Otherwise load from file.',
                        required=False, action = "store_true")
    parser.add_argument('--run_pce',
                        help='run the build PCE from samples and model solutions.  Otherwise load from file.',
                        required=False, action = "store_true")
    parser.add_argument('--run_model_values',
                        help='run the secondary model values for the model, LATs and RMS curves.  Otherwise load from file.',
                        required=False, action = "store_true")
    parser.add_argument('--plot_solutions',
                        help='plot the signals, dv/dt, and activation times for solutions.',
                        required=False, action = "store_true")
    parser.add_argument('--plot_channels', type=int,
                        dest='plot_channels', help='plot all channels a given solution.  The signals, dv/dt, and activation times are plotted.',
                        metavar='', required=False, default = -1)
    parser.add_argument('--make_FIDs',
                        help='an option to manually choose fiducials for the ecg',
                        required=False, action = "store_true")
    parser.add_argument('--hardchecks',
                        help='stop pipeline rather than run parameters implicitly',
                        required=False, action = "store_true")
    return parser
    
def set_dirs(data_dir,resultset,pacing):

    global shape_data_dir, shape_data_outdir, output_dir, tmp_dir, SR_output_file, samples_file, torso_pots_fname, torso_pots, solutions_dir, results_dir
    
    if pacing == "all" or pacing == "old":
        pacing = "sinus"
    elif pacing == "new":
        pacing = "RVV"
    
    output_dir = os.path.join(output_dir_root,data_dir)
    solutions_dir = os.path.join(output_dir,"solutions")
    
    samples_file = os.path.join(output_dir, "UQ_heart_forward_samples.mat")
        
        
    shape_data_dir = os.path.join(shapedata_dir_root,data_dir)
    shape_data_outdir = os.path.join(output_dir,"shape_models/")
    
    pacing_map = {"sinus" : "6105829a_01.120avg-pf.mat", "LV" : "61058472_03.120avg-pf.mat", "apex" : "6105b06d_10.120avg-pf.mat", "RV": "6105b134_16.120avg-pf.mat", "septal": "6105d439_43.120avg-pf.mat", "RVV" : "", "LVV" : "", "Rvot" : "", "Lvot" : "", "RVB" : "" }
    
    resultdir_map = {"pseudoECG": "Solutions_inria_iso/inria_iso_PseudoECG",
          "NGpseudoECG": "Solutions_inria_iso/inria_iso_NG_PseudoECG",
          "inriaLAT" : "Solutions_inria_iso/inria_iso_hires_LAT_sampled" ,
          "ECGsimLAT" : "Solutions_ECGsim/activation_times_sampled" ,
          "ECGsimBSP" : "Solutions_ECGsim/BSP_signals",
          "ECG" : "Solutions_ECGsim/triECG" ,
          "inriaTMP" : "Solutions_inria_iso/inria_iso_hires_sampled"
}
    results_dir = os.path.join(solutions_dir,resultdir_map[resultset])
    
#    torso_pots_fname = "6105d439_43.120avg-pf.mat"
    #torso_pots_fname = "6105829a_01.120avg-pf.mat" # sinus
#    torso_pots_fname = "61058472_03.120avg-pf.mat" # paced
    #torso_pots_fname = "6105b06d_10.120avg-pf.mat" # pace (phrenic)
#    torso_pots_fname = "6105b134_16.120avg-pf.mat" # RV paced

    torso_pots_fname = pacing_map[pacing]
    
    if not torso_pots_fname:
        print(pacing+" not mapped to recorded ECGs")
        return False
        

    torso_pots = os.path.join(torso_pots_dir,torso_pots_fname)
    
    if not os.path.exists(shape_data_outdir):
        os.makedirs(shape_data_outdir)
    
    
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    tmp_dir = os.path.join(output_dir,"tmp")
    
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
        
    SR_output_file = os.path.join(tmp_dir, torso_pots_fname[:-4]+"_UQ_heart_forward_SR_solutions.txt")
    
    
    
    
    
    # get matrix size
    tp = scipy.io.loadmat(torso_pots)
    tp_sz = tp['ts']['potvals'][0][0].shape
#    N = 1024*tp_sz[1]
#
#    pot_size = (1024, tp_sz[1])
    
    return True
    
def get_ranges_from_scores(sample_params, factor = 1.0):
    
    csv_files = glob.glob(shape_data_dir+"/*.csv" )
    
    if len(csv_files)==0:
        print("save scores from shapeworks studio as csv")
        return np.array([])
    elif len(csv_files)>1:
        times = [os.path.getmtime(f) for f in csv_files]
        ind = times.index(max(times))
        csv_files = csv_files[ind]
    else:
        csv_files=csv_files[0]
    
    with open(csv_files, newline='') as csvfile:
        data = list(csv.reader(csvfile))
        
#    print(data)
    
    data_np = np.array(data[1:])
    scores = np.zeros(data_np.shape)
    for k in range(len(data_np)):
        for l in range(len(data_np[0])):
            scores[k,l] = float(data_np[k,l])
#    print(scores)
    std_scores = np.std(scores,axis = 0)
#    print(std_scores)
    
    domain = []
    for k in range(sample_params['dimension']):
        domain.append([-std_scores[k+1]*factor, std_scores[k+1]*factor ])
    
    return np.array(domain).T
    
    
def set_distribution(sample_params):
    dimension = sample_params['dimension']
    alpha = sample_params['alpha']
    beta = sample_params['beta']
    domain = sample_params['domain']
    dist = BetaDistribution(alpha=alpha, beta=beta, dim=dimension, domain=domain)
    # # Expressivity setup
    order = 5
    indices = TotalDegreeSet(dim=dimension, order=order)
    print("number of indices = ",indices.get_indices().shape)
    pce = PolynomialChaosExpansion(indices, dist)
    return pce
    
   
    
def make_runscript():
        # convert samples into transformation matrices, save to disk, load in SCIRun, run model, save solutions to disk, load back to script.
        
    network_file = "/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/run_UQ_model_all.srn5"
    # experiment files.  It is here for now to make some things easier
    
    trans_file = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/staple_submission/geom/trans_to_seg.mat"
    man_trans_file = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/epi_centered/9samples/manual_trans.mat"
    heart_geom = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/epi_centered/9samples/epi_9_0-10.mat"
    torso_geom = "/Users/jess/CIBC/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/Meshes/Dalhousie_torso_morphed.mat"
    torso_elec_geom = "/Users/jess/CIBC/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/Meshes/Dalhousie_torso_subset.mat"
    script_file_tmp = os.path.join(tmp_dir, "UQ_heart_tmp.py")
    
    s_file=open(script_file_tmp,'w+')
    s_file.write("scirun_load_network('"+network_file+"')\n")
    s_file.write("scirun_set_module_state('ReadField:0','Filename','"+heart_geom+"')\n")
    s_file.write("scirun_set_module_state('ReadField:1','Filename','"+torso_geom+"')\n")
    s_file.write("scirun_set_module_state('ReadField:2','Filename','"+torso_elec_geom+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:0','Filename','"+man_trans_file+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:6','Filename','"+trans_file+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:7','Filename','"+torso_pots+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:8','Filename','"+samples_file+"')\n")
    s_file.write("scirun_set_module_state('ResizeMatrix:0','NoOfRows',"+str(N)+")\n")
    s_file.write("scirun_set_module_state('WriteMatrix:1','FileTypeName','SimpleTextFile')\n")
    s_file.write("scirun_set_module_state('WriteMatrix:1','Filename','n_comp.txt')\n")
    s_file.write("scirun_set_module_state('WriteMatrix:0','Filename','"+SR_output_file+"')\n")
    s_file.write("scirun_set_module_state('WriteMatrix:0','FileTypeName','SimpleTextFile')\n")
#    s_file.write("scirun_execute_all()\n")
    s_file.close()
    print( scirun_call+" -0 -S "+script_file_tmp)
    return script_file_tmp, SR_output_file, samples_file
    

def run_scirun(script_file_tmp,SR_output_file):
    output=os.system(scirun_call+" -0 -S "+script_file_tmp)
    return True
    
def load_SR_data(SR_output_file):
    pot_solution = np.loadtxt(SR_output_file)
    return pot_solution.T
    
def load_solutions(resultset, pacing, options):
    
    solutions = load_mat_data_from_dir(resultset, pacing, options)
    
    return solutions
    
def load_mat_data_from_dir(resultset, pacing, options, varname_def=[]):

    if pacing == "default":
        raise ValueError("no solutions for "+pacing+" pacing")
    files = os.listdir(results_dir)
    files.sort()
    
    print(len(files))
        
    ofiles = []
    solutions = []
    data=np.array([])
    
    fid = fid_map[options.pacing]
    plot_fids = False
    if options.make_FIDs or (not fid and not "LAT" in options.resultset):
#        raise ValueError(pacing+" not fully implemented for ECGs. Run --make_FIDs option")
        # only QRS implemented right now
        fid = {"sf": 1, "qon": 0, "qoff": - 1, "ton": 0, "toff": - 1}
        plot_fids = True
    
    for f in files:
        if pacing+"_stim" in f:
#            print(f)
            tmp = scipy.io.loadmat(os.path.join(results_dir,f))
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
            ofiles.append(f)
            
#            qon = int(fid["qon"]/fid["sf"])
#            toff = int(fid["toff"]/fid["sf"])
            
            if 'pseudoECG' == options.resultset:
                qon = int(fid["qon"]/fid["sf"])
                toff = int(fid["toff"]/fid["sf"])
                data = np.transpose(tmp[varname])[:, qon:toff]
                solutions.append(data.flatten()/(4*np.pi))
            elif 'NGpseudoECG' == options.resultset:
                qon = int(fid["qon"]/fid["sf"])
                toff = int(fid["toff"]/fid["sf"])
                data = tmp[varname][:, qon:toff]
                solutions.append(data.flatten()/(4*np.pi))
            elif 'ECGsimBSP' == options.resultset or 'ECG' == options.resultset:
#                data = tmp[varname]
                qon = int(fid["qon"]/fid["sf"])
                toff = int(fid["toff"]/fid["sf"])
                tmp[varname]
                data = tmp[varname][:, qon:toff]
                solutions.append(data.flatten())
            elif 'ECGsimLAT' == options.resultset or 'inriaTMP' == options.resultset:
                data = tmp[varname]
                solutions.append(data.flatten())
            elif 'inriaLAT' == options.resultset:
                data = tmp[varname]
                solutions.append(np.transpose(data).flatten())
            else:
                print("solution file could not be loaded ", os.path.join(results_dir,f))
                return False
                    
    if plot_fids:
        choose_FIDS(solutions, data.shape)
        plt.show()
        return False
        
    if not solutions:
        print("files matching "+pacing+" not found in "+results_dir)
        return False
    
                
#    if 'ECGsim' in results_dir or 'inriaTMP' == options.resultset:
#        sol_size = np.shape(tmp[varname])
#    else:
#        sol_size = np.shape(np.transpose(tmp[varname]))
            
#    print(len(solutions))
#    print(np.vstack(solutions).shape)
#    print(ofiles)
    
    return solutions, data.shape, ofiles
    
def compute_surface_variance(domain):
    eig_vect_files = ["eigenvectors0.eval", "eigenvectors1.eval", "eigenvectors2.eval", "eigenvectors3.eval", "eigenvectors4.eval"]
    print("computing surface variance")
    magns = []
    for k in range(len(domain[0])):
        evec = np.loadtxt(os.path.join(shape_data_dir, eig_vect_files[k]))*domain[0][k]
        
        mag = np.sqrt(np.sum(evec*evec, axis=1))
        
        print(mag.shape)
        magns.append(mag)
        
    magnitude = np.vstack(magns)
    print(magnitude.shape)
    
    var =np.sum(magnitude*magnitude,axis = 0)
    
    print(var.shape)
    file = os.path.join(output_dir, "surface_variance.mat")
    print(file)
    scipy.io.savemat(file, dict(var = var))
    
    return True
    
    
def correct_samples(samples,pce,filenames):
    print(samples.shape)
#    print(filenames)
    
    ind = []
    
    for f in filenames:
        ind.append(int(re.search(r'\d+', f).group()))
        
#    print(ind)
    print(len(ind))
    length_index = len(pce.index_set.get_indices())
    if len(ind)<=(length_index):
        print("number of actual samples less than number of indices.  a list of new parameters to run:")
        pce.generate_samples()
        new_samples = pce.samples
        print(type(new_samples))
        print(len(new_samples))
        r_ind = random.sample(range(len(new_samples)),length_index-len(ind)+10)
        print(new_samples[r_ind,:])
        print("Append to samples file? [y/N]")
        str_input = str(input())
        if str_input.lower()=='y' or str_input.lower()=='yes':
            print("saving to "+samples_file)
            scipy.io.savemat(samples_file, dict(samples=np.vstack([samples,new_samples[r_ind,:]])))
            print("exiting")
            return []
        else:
            print("exiting")
            return []

    
    return samples[ind,:]
    
def sample_workaround(samples, domain):
    print('executing workaround')
    samples_new = samples/np.dot(np.ones((np.shape(samples)[0],1)),domain[0,:].reshape((1,5)))
    return samples_new
    
#def cleanup_run(script_file_tmp,SR_output_file
#    #    os.remove(script_file_tmp)
#    #    os.remove(SR_output_file)
#    return True

def findMinDVDT(signals, samp, window=20, deg=3, ingnore_first = 100, fids = {}):
    [M,T] = signals.shape
    tau = np.zeros((M))
    dy = np.zeros(signals.shape)
    
    ig_factor = 0.5
    if fids:
        qoff = int(fids["qoff"]*fids["sf"][0,0]/1000 - ignore_first*ig_factor)
#        print(fids["qoff"], qoff)
    else:
        qoff = int(ignore_first*ig_factor)
    for k in range(M):
        cent = math.ceil(window/2)
        X = np.zeros((window,deg+1))
        L = np.arange(-cent,cent)
        for p in range(1,deg+2):
            X[:,p-1] = L**((deg+1)-p)
        E = np.dot(np.linalg.inv(np.dot(X.T,X)),X.T)
        temp = np.hstack((signals[k,:],signals[k,-1]*np.ones((cent))))
        a = sig.lfilter(E[deg-1,np.arange(window-1,-1,-1)],1,temp)
#        plt.plot(temp)
#        plt.plot(a)
#        plt.show()
        dy[k,:] = a[cent-1:-1]/samp
        tau[k] = (np.argmin(dy[k,ingnore_first:])+ingnore_first)/samp
    return tau, dy
    
def phaseAnalysis(pots, sample_rate = 1):
    #TODO
    #ref:
    #https://www.ahajournals.org/doi/full/10.1161/circep.110.853804


def find_LAT(pots, sample_rate=1):
    #512 nodes x 755
    lats = []
    dys = []
    for p in pots:
        p = np.array(p)
        signals = np.resize(p,(512,tp_sz[1]))
        [lat,dy]=findMinDVDT(signals, sample_rate, 80, 3, 50)
        lats.append(lat)
        dys.append(dy)
    return np.array(lats), np.array(dys)*sample_rate
    
def model_RMS(model_output, pot_size, axis = 0):
    RMS = np.zeros((len(model_output),pot_size[axis-1]))
    for k in range(len(model_output)):
        m = np.resize(model_output[k], pot_size)
        RMS[k,:] = np.sqrt(np.sum(m*m/pot_size[axis],axis=axis))
    return RMS
        
    
def plot_test(model_data, lat, dys, channel):
#    signals = np.resize(model_output[3],(512,tp_sz[1]))
#    [lat,dy]=findMinDVDT(signals, 2000, 10, 2)
#    plt.plot(np.arange(0,tp_sz[1]/2000,1/2000),signals[channel,:])
#    plt.plot([lat[channel],lat[channel]],[np.min(signals[channel,:]), np.max(signals[channel,:])])
#    plt.plot(np.arange(0,tp_sz[1]/2000,1/2000),dy[channel,:])
#    plt.show()
    return
    
def plot_ecg_uq(time, median, quantiles, ecgs, ecg_size, num_samples= 10, title = ""):

    if len(time)==0:
        time = list(range(ecg_size[1]))
    rand_ind = random.sample(range(len(ecgs)),num_samples)
#    print(rand_ind)
    signals = []
    for r in rand_ind:
        signals.append(np.resize(ecgs[r],ecg_size))
    
    lead_nums  = [21, 35, 43, 50, 64, 67]

#    band_mass = 1/(2*(Q+1))
#    for ind in range(Q):
#        alpha = (Q-ind) * 1/Q - (1/(2*Q))
#        if ind == 0:
#            plt.fill_between(time,quantiles_RMS[ind,:], quantiles_RMS[Q+ind, :],
#                             interpolate=True, facecolor='red', alpha=alpha,
#                             label='{0:1.2f} probability mass (each band)'.format(band_mass))
#        else:
#            plt.fill_between(time,quantiles_RMS[ind, :], quantiles_RMS[Q+ind, :], interpolate=True, facecolor='red', alpha=alpha)
            
#    print(quantiles.shape)
    Q=4
    band_mass = 1/(2*(Q+1))
    
    time = list(range(ecg_size[1]))
   
    fig, ax = plt.subplots(6,1)
    
    for l in range(len(lead_nums)):
        ax[l].plot(time, median[lead_nums[l],:], '-k',label='median')
        ax[l].set(ylabel='v'+str(l)+'(mV)')
        for ind in range(Q):
            alpha = (Q-ind) * 1/Q - (1/(2*Q))
            if ind == 0:
#                print(quantiles[ind,lead_nums[l],:]-quantiles[Q+ind,lead_nums[l],:])
                ax[l].fill_between(time,quantiles[ind,lead_nums[l],:],quantiles[Q+ind,lead_nums[l],:],interpolate=True, facecolor='red', alpha=alpha, label='{0:1.2f} probability mass (each band)'.format(band_mass))
            else:
                ax[l].fill_between(time, quantiles[ind,lead_nums[l],:],quantiles[Q+ind,lead_nums[l],:],interpolate=True, facecolor='red', alpha=alpha)
        for s in range(num_samples):
            tmp = signals[s]
            ax[l].plot(time,tmp[lead_nums[l],:], ':b', linewidth=0.5)
    ax[l].set(xlabel = "time (ms)")
    fig.suptitle(title)

#        ax[l].plot(np.arange(0,mat_size[1]/2,0.5),signals[k+l,:])
#            ax[].plot([lat[sample_num,k+l]*1000,lat[sample_num,k+l]*1000],[np.min(signals[k+l,:]), np.max(signals[k+l,:])])
#            ax[r,c].plot(np.arange(0,mat_size[1]/2,0.5),dy[k+l,:])
#            ax[r,c].set(ylabel="lead "+str(k+l))
#        fig.suptitle("sample "+str(sample_num)+", leads "+str(k)+" - "+str(np.min(np.min([k+num_plots, mat_size[0]]))))
##        plt.show()
    return
    
def plot_all_solutions(model_output, lat, dys, mat_size, channel_nums):
    model_output = model_output.reshape((len(model_output),mat_size[0],mat_size[1]))
    if isinstance(channel_nums, (list, tuple, np.ndarray)):
        n_row = 10
        n_col = 5
        num_plots = n_row*n_col
        fig_size = (18, 10)
        fg, ax = plt.subplots(n_row,n_col,figsize = fig_size)
        for l in range(num_plots):
            if l>=mat_size[0]:
                break
            c = int(np.floor(l/n_row))
            r = int(l-c*(n_row))
            plot_all_solutions_one(model_output[:,channel_nums[l],:], lat[:,channel_nums[l]], dys[:,channel_nums[l],:], channel_nums[l], ax[r,c])
    else:
        fg,ax = plt.subplots(1)
        plot_all_solutions_one(model_output[:,channel_nums,:], lat[:,channel_nums], dys[:,channel_nums,:], channel_nums, ax)
    fg.suptitle("all solutions")
    return
     
def plot_all_solutions_one(signals, lat, dys, channel, ax):
    sig_size = signals.shape
    ax.plot(np.arange(0,sig_size[1]/2,0.5),dys.T, 'r')
    ax.plot(np.arange(0,sig_size[1]/2,0.5),signals.T, 'b')
    ax.vlines(lat*1000,np.min(signals), np.max(signals))
    ax.set(ylabel="lead "+str(channel)+" (mV)")
    ax.set(xlabel="time (ms)")
    return
    
def choose_FIDS(model_output,model_size):

    RMS = model_RMS(model_output, model_size)
    
    sol_num = 1
    lead_num = 5
    
    fg, ax = plt.subplots(2,1)
    print(model_size)
    print(type(model_output))
    print(len(model_output))
    print(model_output[0].shape)
    print(RMS.shape)
    print(np.vstack(model_output).shape)
    ax[0].plot(np.arange(0,model_size[1]/2,0.5),RMS[sol_num,:], 'r')
    model_output = np.vstack(model_output).reshape((len(model_output),model_size[0],model_size[1]))
    ax[1].plot(np.arange(0,model_size[1]/2,0.5),model_output[sol_num,:,:].T)
        
    return
    
        
    
    
            
            
    




#plt.figure()
#
##plt.plot(model_output[:V, :].T, 'k', alpha=0.8, linewidth=0.2)
#plt.plot(model_output[0:383],median[lead_num,:], 'b', label='PCE median')
#
#band_mass = 1/(2*(Q+1))
#
#for ind in range(Q):
#    alpha = (Q-ind) * 1/Q - (1/(2*Q))
#    if ind == 0:
#        plt.fill_between(np.array(range(383)),quantiles[ind,lead_num, :], quantiles[Q+ind,lead_num, :],
#                         interpolate=True, facecolor='red', alpha=alpha,
#                         label='{0:1.2f} probability mass (each band)'.format(band_mass))
#    else:
#        plt.fill_between(np.array(range(383)),quantiles[ind,lead_num, :], quantiles[Q+ind,lead_num, :], interpolate=True, facecolor='red', alpha=alpha)
#        
#plt.show()


#[[-110.858, 15.0996, -10.661,    16.2467,    20.3813,    -0.0670236,    13.5023,   -10.3614,   31.0718],
#[56.5006,    27.2608,    -27.2569,    -10.138,    11.8904,    -7.73922,    -22.1777, -13.3565, -18.8841],
#[-92.5537, -61.5309, -11.5002,    -19.1763,    -8.74692,    -0.827493,    -14.9879,    7.691,    -6.87735],
#[76.4543, -55.1819, -35.4948,  7.01125,    -2.81745,    18.0575,    14.0494,    -1.68374,    -5.01133],
#[25.6006,    -8.0849,    -1.19431,    20.4601,    -15.6585,    -34.1141,    6.46472,    3.69296,    -4.5335],
#[-12.2401,    81.3567,    -21.6332,    -13.1139,    -16.7665,    9.01392,    6.22342,    10.9451,    -17.585],
#[6.13693,    -3.04349,    48.4903,    -18.6867,    -18.8935,    4.69597,    6.53618,    -15.3208,    9.98601],
#[42.3389,    -4.6396,    28.6189,    -19.4726,    30.8781,    -4.95629,    7.17589,    12.4997,    16.3259],
#[8.61998,    8.7637,    30.6312,    36.8695,    -0.266909,    15.9367,    -16.7863,    5.89374,    -4.49245]]

#plt.hist(epi_coeffs[:,0], bins = 3), plt.show()

def run_pipeline(resultset, pacing, sample_params, options):

    set_dirs(options.data_dir, resultset, pacing)
    
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
    options = parser.parse_args()

#    global resultset, pacing
        
    pacing = options.pacing
    resultset = options.resultset
    
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
