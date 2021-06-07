from itertools import chain, combinations
import scipy.io
import os

import numpy as np
import math
import scipy.signal as sig



import matplotlib.pyplot as plt


#place for UQ files
output_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data"
torso_pots_dir = "/Users/jess/CIBC/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/Interventions/BSPM_pf/"
#torso_pots_fname = "6105d439_43.120avg-pf.mat"
# sinus, paced, phrenic paced, rv paced
torso_pots_fname = ["6105829a_01.120avg-pf.mat",  "61058472_03.120avg-pf.mat", "6105b06d_10.120avg-pf.mat", "6105b134_16.120avg-pf.mat" ]
labels = ["Sinus", "LV Paced", "Apically Paced", "RV Paced"]

#torso_pots = os.path.join(torso_pots_dir,torso_pots_fname)
tmp_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/tmp/"


#UQ_file = os.path.join(output_dir, torso_pots_fname[:-4]+"_heart_shape_UQ_values.mat")
#LAT_UQ_file = os.path.join(output_dir, torso_pots_fname[:-4]+"_heart_shape_LAT_UQ_values.mat")
#RMS_UQ_file = os.path.join(output_dir, torso_pots_fname[:-4]+"_heart_shape_RMS_UQ_values.mat")


# get matrix size
tp = scipy.io.loadmat(os.path.join(torso_pots_dir,torso_pots_fname[0]))
tp_sz = tp['ts']['potvals'][0][0].shape
            
    
pot_size = (512, tp_sz[1])
    
        
# run model to get solutions

RMS_UQ_file = os.path.join(output_dir, torso_pots_fname[0][:-4]+"_heart_shape_RMS_UQ_values.mat")
tmp = scipy.io.loadmat(RMS_UQ_file)
median_RMS = tmp["median"][None]
quantiles_RMS = tmp["quantiles"].T

Q = int((len(quantiles_RMS)-1)/2)
band_mass = 1/(2*(Q+1))

dt = 0.5







fig, ax = plt.subplots(1,len(torso_pots_fname))

fig2, ax2 = plt.subplots(1,len(torso_pots_fname))





for k in range(len(torso_pots_fname)):
    UQ_file = os.path.join(output_dir, torso_pots_fname[k][:-4]+"_heart_shape_UQ_values.mat")
    tmp = scipy.io.loadmat(RMS_UQ_file)
    global_sensitivity = tmp["glob_sensitivity"]
    print(np.mean(global_sensitivity,axis=0))
    RMS_UQ_file = os.path.join(output_dir, torso_pots_fname[k][:-4]+"_heart_shape_RMS_UQ_values.mat")
    tmp = scipy.io.loadmat(RMS_UQ_file)
    median_RMS = tmp["median"].reshape(len(tmp["median"].T))
    quantiles_RMS = tmp["quantiles"].T
    
    time = np.arange(0,len(median_RMS)*dt,dt)
    
    ax[k].plot(time, median_RMS/1000, 'b', label='RMS median')
    

    for ind in range(Q):
        alpha = (Q-ind) * 1/Q - (1/(2*Q))
        if ind == 0:
            ax[k].fill_between(time,quantiles_RMS[ind,:]/1000, quantiles_RMS[Q+ind, :]/1000,
                             interpolate=True, facecolor='red', alpha=alpha,
                             label='{0:1.2f} probability mass (each band)'.format(band_mass))
        else:
            ax[k].fill_between(time,quantiles_RMS[ind, :]/1000, quantiles_RMS[Q+ind, :]/1000, interpolate=True, facecolor='red', alpha=alpha)
            
    ax[k].set(ylabel = labels[k]+" RMS (mV)")
    ax[k].set(xlabel = "time (ms)")
#ax[-1].set(xlabel = "time (ms)")

    global_sensitivity_RMS = tmp["glob_sensitivity"]
    
#    print(np.median(global_sensitivity_RMS,axis=0))
    
    ax2[k].plot(time, global_sensitivity_RMS)
    
    ax2[k].set(ylabel = labels[k]+" RMS (mV)")
    ax2[k].set(xlabel = "time (ms)")
ax2[-1].legend(list(range(len(ax2[1].lines))), bbox_to_anchor= (1.1,1), ncol = 2)

#bbox_to_anchor= (1.5,-0.2),

plt.show()

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
