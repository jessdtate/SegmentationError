from itertools import chain, combinations
import scipy.io
import os

import numpy as np

from UncertainSCI.distributions import BetaDistribution
from UncertainSCI.indexing import TotalDegreeSet
from UncertainSCI.pce import PolynomialChaosExpansion

#place for UQ files
output_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data"

# Number of parameters
dimension = 4

# Specifies 1D distribution on [0,1] (alpha=beta=1 ---> uniform)
alpha = 1
beta = 1.
# domain is the range of the hypercube
domain = np.array([[-125, 125], [-85, 85], [-60, 60], [-40, 40]]).T
dist = BetaDistribution(alpha=alpha, beta=beta, dim=dimension, domain=domain)

# # Expressivity setup
order = 3
indices = TotalDegreeSet(dim=dimension, order=order)
pce = PolynomialChaosExpansion(indices, dist)
pce.generate_samples()

#samples = pce.samples
#print(samples)


# solution points (512 nodes x 755 time points)
N = 386560

    
    
def make_runscript(samples,N):
        # convert samples into transformation matrices, save to disk, load in SCIRun, run model, save solutions to disk, load back to script.
#    scirun_call = "/Users/jess/software/SCIRun/bin_clean/SCIRun/SCIRun_test"
    scirun_call = "/Applications/SCIRun.app/Contents/MacOS/SCIRun"
    network_file = "/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/run_UQ_model_all.srn5"
    tmp_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/tmp/"
    # experiment files.  It is here for now to make some things easier
    torso_pots = "/Users/jess/CIBC/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/Interventions/BSPM_pf/6105d439_43.120avg-pf.mat"
    trans_file = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/staple_submission/geom/trans_to_seg.mat"
    man_trans_file = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/epi_centered/9samples/manual_trans.mat"
    heart_geom = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/epi_centered/9samples/epi_9_0-10.mat"
    torso_geom = "/Users/jess/CIBC/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/Meshes/Dalhousie_torso_morphed.mat"
    torso_elec_geom = "/Users/jess/CIBC/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/Meshes/Dalhousie_torso_subset.mat"
    script_file_tmp = os.path.join(tmp_dir, "UQ_heart_tmp.py")
    SR_output_file = os.path.join(tmp_dir, "UQ_heart_SR_solutions.txt")
    samples_file = os.path.join(tmp_dir, "UQ_heart_samples.mat")
    scipy.io.savemat(samples_file, dict(samples=samples))
    s_file=open(script_file_tmp,'w+')
    s_file.write("scirun_load_network('"+network_file+"')\n")
    s_file.write("scirun_set_module_state('ReadField:0','Filename','"+heart_geom+"')\n")
    s_file.write("scirun_set_module_state('ReadField:1','Filename','"+torso_geom+"')\n")
    s_file.write("scirun_set_module_state('ReadField:2','Filename','"+torso_elec_geom+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:0','Filename','"+man_trans_file+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:6','Filename','"+trans_file+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:7','Filename','"+torso_pots+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:8','Filename','"+samples_file+"')\n")
    s_file.write("scirun_set_module_state('WriteMatrix:1','FileTypeName','SimpleTextFile')\n")
    s_file.write("scirun_set_module_state('WriteMatrix:1','Filename','n_comp.txt')\n")
    s_file.write("scirun_set_module_state('WriteMatrix:0','Filename','"+SR_output_file+"')\n")
    s_file.write("scirun_set_module_state('WriteMatrix:0','FileTypeName','SimpleTextFile')\n")
#    s_file.write("scirun_execute_all()\n")
    s_file.close()
    print( scirun_call+" -0 -S "+script_file_tmp)
    return script_file_tmp, SR_output_file
    

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
    
    
    
    
    
    

    
# run model to get solutions

prep_run = make_runscript(pce.samples,N)

input("press enter when SCIRun net is finished.")

model_output = load_SR_data(prep_run[1])

pce.build(model_output=model_output)

mean = pce.mean()
stdev = pce.stdev()

# Power set of [0, 1, ..., dimension-1]
variable_interactions = list(chain.from_iterable(combinations(range(dimension), r) for r in range(1, dimension+1)))

# "Total sensitivity" is a non-partitive relative sensitivity measure per parameter.
total_sensitivity = pce.total_sensitivity()

# "Global sensitivity" is a partitive relative sensitivity measure per set of parameters.
global_sensitivity = pce.global_sensitivity(variable_interactions)

Q = 3  # Number of quantile bands to plot
dq = 0.5/(Q+1)
q_lower = np.arange(dq, 0.5-1e-7, dq)[::-1]
q_upper = np.arange(0.5 + dq, 1.0-1e-7, dq)
quantile_levels = np.append(np.concatenate((q_lower, q_upper)), 0.5)

quantiles = pce.quantile(quantile_levels, M=int(2e3))

median = pce.quantile(0.5, M=int(1e3))[0, :]


UQ_file = os.path.join(output_dir, "hear_shape_UQ_values.mat")

scipy.io.savemat(UQ_file, dict(mean = mean.T, stdev = stdev.T, tot_sensitivity = total_sensitivity.T, glob_sensitivity = global_sensitivity.T, quantiles = quantiles.T, median = median.T))
    

lead_num = 404

median = median.reshape((512,755))
quantiles = quantiles.reshape((quantiles.shape[0],512,755))

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
