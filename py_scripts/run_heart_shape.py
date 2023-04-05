from itertools import chain, combinations
import scipy.io
import os
import argparse

import numpy as np
import math
import scipy.signal as sig

from UncertainSCI.distributions import BetaDistribution
from UncertainSCI.indexing import TotalDegreeSet
from UncertainSCI.pce import PolynomialChaosExpansion

import matplotlib.pyplot as plt


#place for UQ files
output_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data"
torso_pots_dir = "/Users/jess/CIBC/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/Interventions/BSPM_pf/"
torso_pots_fname = "6105d439_43.120avg-pf.mat"
#torso_pots_fname = "6105829a_01.120avg-pf.mat" # sinus
torso_pots_fname = "61058472_03.120avg-pf.mat" # LV paced
#torso_pots_fname = "6105b06d_10.120avg-pf.mat" # pace (phrenic)
#torso_pots_fname = "6105b134_16.120avg-pf.mat" # RV paced

torso_pots = os.path.join(torso_pots_dir,torso_pots_fname)

tmp_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/tmp/"
SR_output_file = os.path.join(tmp_dir, torso_pots_fname[:-4]+"_UQ_heart_SR_solutions.txt")
samples_file = os.path.join(tmp_dir, "UQ_heart_samples.mat")

# get matrix size
tp = scipy.io.loadmat(torso_pots)
tp_sz = tp['ts']['potvals'][0][0].shape
N = 512*tp_sz[1]

pot_size = (512, tp_sz[1])




domain = np.array([[-125, 125], [-85, 85], [-60, 60], [-40, 40]]).T
sample_params = { "dimension": 4, "alpha": 1, "beta": 1, "domain": domain}



def build_parser():
    parser = argparse.ArgumentParser()

    # This will be implemented as rollout broadens
    parser.add_argument('--compute_samples', type=bool,
                        dest='compute_samples', help='compute samples to run in model.  Otherwise load from file.',
                        metavar='', required=False, default = False)
                        
    parser.add_argument('--run_model', type=bool,
                        dest='run_model', help='run the model with the samples.  Otherwise load from file.',
                        metavar='', required=False, default = False)
                        
    parser.add_argument('--run_pce', type=bool,
                        dest='run_pce', help='run the build PCE from samples and model solutions.  Otherwise load from file.',
                        metavar='', required=False, default = False)
    parser.add_argument('--run_model_values', type=bool,
                        dest='run_model_values', help='run the secondary model values for the model, LATs and RMS curves.  Otherwise load from file.',
                        metavar='', required=False, default = False)
                        
    parser.add_argument('--plot_solutions', type=bool,
                        dest='plot_solutions', help='plot the signals, dv/dt, and activation times for solutions.',
                        metavar='', required=False, default = False)
                        
    parser.add_argument('--plot_channels', type=int,
                        dest='plot_channels', help='plot all channels a given solution.  The signals, dv/dt, and activation times are plotted.',
                        metavar='', required=False, default = -1)
    return parser
    
    
def make_samples(Filename, pce):
    pce.generate_samples()
    scipy.io.savemat(Filename, dict(samples=pce.samples))
    return True
    
def set_distribution(sample_params):
    dimension = sample_params['dimension']
    alpha = sample_params['alpha']
    beta = sample_params['beta']
    domain = sample_params['domain']
    dist = BetaDistribution(alpha=alpha, beta=beta, dim=dimension, domain=domain)
    # # Expressivity setup
    order = 5
    indices = TotalDegreeSet(dim=dimension, order=order)
    pce = PolynomialChaosExpansion(indices, dist)
    return pce
   
    
def make_runscript():
        # convert samples into transformation matrices, save to disk, load in SCIRun, run model, save solutions to disk, load back to script.
#    scirun_call = "/Users/jess/software/SCIRun/bin_clean/SCIRun/SCIRun_test"
    scirun_call = "/Applications/SCIRun.app/Contents/MacOS/SCIRun"
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
    
#def cleanup_run(script_file_tmp,SR_output_file
#    #    os.remove(script_file_tmp)
#    #    os.remove(SR_output_file)
#    return True

def findMinDVDT(signals, samp, window=20, deg=3, ingnore_first = 100 ):
    [M,T] = signals.shape
    tau = np.zeros((M))
    dy = np.zeros(signals.shape)
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
    
def model_RMS(model_output, pot_size):
    RMS = np.zeros((len(model_output),pot_size[1]))
    for k in range(len(model_output)):
        m = np.resize(model_output[k], pot_size)
        RMS[k,:] = np.sqrt(np.sum(m*m,axis=0))
    return RMS
        
    
def plot_test(model_data, lat, dys, channel):
#    signals = np.resize(model_output[3],(512,tp_sz[1]))
#    [lat,dy]=findMinDVDT(signals, 2000, 10, 2)
#    plt.plot(np.arange(0,tp_sz[1]/2000,1/2000),signals[channel,:])
#    plt.plot([lat[channel],lat[channel]],[np.min(signals[channel,:]), np.max(signals[channel,:])])
#    plt.plot(np.arange(0,tp_sz[1]/2000,1/2000),dy[channel,:])
#    plt.show()
    return
    
def plot_all_channels(model_output, lat, dys, mat_size, sample_num):
    signals = np.resize(model_output[sample_num],(mat_size[0],mat_size[1]))
    dy = dys[sample_num,:,:]
    n_row = 10
    n_col = 5
    num_plots = n_row*n_col
    fig_size = (18, 10)
    for k in range(0,mat_size[0],num_plots):
        fig, ax = plt.subplots(n_row,n_col,figsize = fig_size)
        for l in range(num_plots):
            if k+l>=mat_size[0]:
                break
            c = int(np.floor(l/n_row))
            r = int(l-c*(n_row))
            ax[r,c].plot(np.arange(0,mat_size[1]/2,0.5),signals[k+l,:])
            ax[r,c].plot([lat[sample_num,k+l]*1000,lat[sample_num,k+l]*1000],[np.min(signals[k+l,:]), np.max(signals[k+l,:])])
            ax[r,c].plot(np.arange(0,mat_size[1]/2,0.5),dy[k+l,:])
            ax[r,c].set(ylabel="lead "+str(k+l))
        fig.suptitle("sample "+str(sample_num)+", leads "+str(k)+" - "+str(np.min(np.min([k+num_plots, mat_size[0]]))))
#        plt.show()
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
#    fg.suptitle("all solutions")
    return
     
def plot_all_solutions_one(signals, lat, dys, channel, ax):
    sig_size = signals.shape
    ax.plot(np.arange(0,sig_size[1]/2,0.5),dys.T, 'r')
    ax.plot(np.arange(0,sig_size[1]/2,0.5),signals.T, 'b')
    ax.vlines(lat*1000,np.min(signals), np.max(signals))
    ax.set(ylabel="lead "+str(channel))
    ax.set(xlabel="time (ms)")
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


def main():

    parser = build_parser()
    options = parser.parse_args()
    
    pce = set_distribution(sample_params)
    
    if options.compute_samples:
        make_samples(samples_file, pce)
        print("Samples compute.  Make sure to run the model")
        return
        
        
        
        
    if options.run_model:
        prep_run=make_runscript()
        input("press enter when SCIRun net is finished.")
        model_output = load_SR_data(SR_output_file)
        scipy.io.savemat(SR_output_file[:-4]+".mat",{"model_solutions" : model_output, "mat_size": pot_size})
        options.run_pce = True
        options.run_model_values = True
    else:
        tmp = scipy.io.loadmat(SR_output_file[:-4]+".mat")
        model_output  = tmp["model_solutions"]
        del tmp
        
#    model_output = load_SR_data(SR_output_file)
#    scipy.io.savemat(SR_output_file[:-4]+".mat",{"model_solutions" : model_output, "mat_size": pot_size})

    
    print("SCIRun solutions loaded")
    
    tmp_samp = scipy.io.loadmat(samples_file)
    #    print(tmp_samp )
    samples = tmp_samp['samples']
    del tmp_samp
    print("loaded samples")

    

    LAT_file = os.path.join(tmp_dir, torso_pots_fname[:-4]+"_UQ_heart_LAT_solutions.mat")

    if options.run_model_values:
        LATs, dys = find_LAT(model_output,2000)
        scipy.io.savemat(LAT_file, dict(LATs = LATs, dys = dys, sample_rate = 2000))
        options.run_pce = True
        print("LATs computed")
        
    else:
        tmp = scipy.io.loadmat(LAT_file)
        LATs = tmp["LATs"]
        dys = tmp["dys"]
        del tmp
        print("LATs loaded")
    



    UQ_file = os.path.join(output_dir, torso_pots_fname[:-4]+"_heart_shape_UQ_values.mat")
    LAT_UQ_file = os.path.join(output_dir, torso_pots_fname[:-4]+"_heart_shape_LAT_UQ_values.mat")
    RMS_UQ_file = os.path.join(output_dir, torso_pots_fname[:-4]+"_heart_shape_RMS_UQ_values.mat")
    
    Q = 4  # Number of quantile bands to plot

    
    
    if options.run_pce:
        RMS = model_RMS(model_output, pot_size)
        pce.build(model_output=np.hstack((model_output,RMS,LATs)), samples = samples)
        print("PCE model built")


        mean = pce.mean()[:N]
        stdev = pce.stdev()[:N]
        mean_RMS = pce.mean()[N:N+pot_size[1]]
        stdev_RMS = pce.stdev()[N:N+pot_size[1]]
        mean_LAT = pce.mean()[N+pot_size[1]:]
        stdev_LAT = pce.stdev()[N+pot_size[1]:]



    # Power set of [0, 1, ..., dimension-1]
        variable_interactions = list(chain.from_iterable(combinations(range(sample_params['dimension']), r) for r in range(1, sample_params['dimension']+1)))

    # "Total sensitivity" is a non-partitive relative sensitivity measure per parameter.
        total_sensitivity = pce.total_sensitivity()[:,:N]
        total_sensitivity_RMS = pce.total_sensitivity()[:,N:N+pot_size[1]]
        total_sensitivity_LAT = pce.total_sensitivity()[:,N+pot_size[1]:]

        # "Global sensitivity" is a partitive relative sensitivity measure per set of parameters.
        global_sensitivity = pce.global_sensitivity(variable_interactions)[:,:N]
        global_sensitivity_RMS = pce.global_sensitivity(variable_interactions)[:,N:N+pot_size[1]]
        global_sensitivity_LAT = pce.global_sensitivity(variable_interactions)[:,N+pot_size[1]:]

        dq = 0.5/(Q+1)
        q_lower = np.arange(dq, 0.5-1e-7, dq)[::-1]
        q_upper = np.arange(0.5 + dq, 1.0-1e-7, dq)
        quantile_levels = np.append(np.concatenate((q_lower, q_upper)), 0.5)

        quantiles = pce.quantile(quantile_levels, M=int(2e3))
        quantiles_LAT = quantiles[:,N+pot_size[1]:]
        quantiles_RMS = quantiles[:,N:N+pot_size[1]]
        quantiles = quantiles[:,:N]

        median = pce.quantile(0.5, M=int(1e3))[0, :]
        median_LAT = median[N+pot_size[1]:]
        median_RMS = median[N:N+pot_size[1]]
        median = median[:N]

        print("metrics computed")
        
        print(LATs.shape)
        print(np.max(LATs,axis=0)-np.min(LATs,axis=0))

        scipy.io.savemat(UQ_file, dict(mean = mean.T, stdev = stdev.T, tot_sensitivity = total_sensitivity.T, glob_sensitivity = global_sensitivity.T, quantiles = quantiles.T, median = median.T, pot_size = np.array(pot_size)))

        scipy.io.savemat(LAT_UQ_file, dict(mean = mean_LAT.T, stdev = stdev_LAT.T, tot_sensitivity = total_sensitivity_LAT.T, glob_sensitivity = global_sensitivity_LAT.T, quantiles = quantiles_LAT.T, median = median_LAT.T, LAT_min = np.min(LATs,axis=0), LAT_max = np.max(LATs,axis=0)))

        scipy.io.savemat(RMS_UQ_file, dict(mean = mean_RMS.T, stdev = stdev_RMS.T, tot_sensitivity = total_sensitivity_RMS.T, glob_sensitivity = global_sensitivity_RMS.T, quantiles = quantiles_RMS.T, median = median_RMS.T))
        
    else:
        print("loading pce stats from disk")
        tmp = scipy.io.loadmat(UQ_file)
        mean = tmp["mean"].T
        stdev = tmp["stdev"].T
        median =tmp["median"].T
        total_sensitivity =tmp["tot_sensitivity"].T
        glob_sensitivity = tmp["glob_sensitivity"].T
        quantiles = tmp["quantiles"].T
        
        tmp = scipy.io.loadmat(LAT_UQ_file)
        mean_LAT = tmp["mean"].T
        stdev_LAT = tmp["stdev"].T
        median_LAT =tmp["median"].T
        total_sensitivity_LAT =tmp["tot_sensitivity"].T
        glob_sensitivity_LAT = tmp["glob_sensitivity"].T
        quantiles_LAT = tmp["quantiles"].T
        
        
        tmp = scipy.io.loadmat(RMS_UQ_file)
        mean_RMS = tmp["mean"].T
        stdev_RMS = tmp["stdev"].T
        median_RMS =tmp["median"].T
        total_sensitivity_RMS =tmp["tot_sensitivity"].T
        glob_sensitivity_RMS = tmp["glob_sensitivity"].T
        quantiles_RMS = tmp["quantiles"].T
        print("loaded")
        


    dt = 0.5

    lead_num = 404

    median = median.reshape(pot_size)
    quantiles = quantiles.reshape((quantiles.shape[0],pot_size[0],pot_size[1]))

    time = np.arange(0,pot_size[1]*dt,dt)
    RMS_med = np.sqrt(np.sum(median*median,axis=0))
    RMS_quant = np.sqrt(np.sum(quantiles*quantiles,axis=1))

    if options.plot_channels>=0:
        plot_all_channels(model_output, LATs, dys, pot_size, options.plot_channels)
     
    if options.plot_solutions:
        channels = (stdev_LAT>3*np.mean(stdev_LAT)).nonzero()[0]
        plot_all_solutions(model_output, LATs, dys*50, pot_size, channels)
        


    plt.figure()
    plt.plot(time, median_RMS, 'b', label='RMS median')

    band_mass = 1/(2*(Q+1))
    for ind in range(Q):
        alpha = (Q-ind) * 1/Q - (1/(2*Q))
        if ind == 0:
            plt.fill_between(time,quantiles_RMS[ind,:], quantiles_RMS[Q+ind, :],
                             interpolate=True, facecolor='red', alpha=alpha,
                             label='{0:1.2f} probability mass (each band)'.format(band_mass))
        else:
            plt.fill_between(time,quantiles_RMS[ind, :], quantiles_RMS[Q+ind, :], interpolate=True, facecolor='red', alpha=alpha)


    plt.show()
    
    return
    

if __name__ == "__main__":
    main()
