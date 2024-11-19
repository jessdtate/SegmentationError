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


from ECGtools import *
from EPUQtools import *

#place for UQ files
exec(open("initValues.py").read())
    
    

def build_parser():
  parser = argparse.ArgumentParser()

  # This will be implemented as rollout broadens
  parser.add_argument('--data_dir', required=False,
                      help='directory where the shape data is located, relative to the shape_data dir',
                      default = default_datadir)
  parser.add_argument('--pacing', required=False,
                      help='pacing profile to run.  {[sinus], apex, RV, LV, septal, RVV, LVV, Rvot, Lvot, RVB, all, old, or new}',
                      default = "sinus")
  parser.add_argument('--distrib', required=False,
                      help='assumed parameter distributions.  {[uniform], gaussian}',
                      default = "uniform")
  parser.add_argument('--compute_samples',
                      help='compute samples to run in model.  Otherwise load from file.',
                      required=False, action = "store_true")
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
  parser.add_argument('--save',
                      help='do not save plots and tables to disk',
                      required=False, action = "store_false")
  return parser

def set_dirs(options):

  data_dir = options.data_dir
  resultset = options.resultset
  pacing = options.pacing

  if pacing == "all" or pacing == "old":
    pacing = "sinus"
  elif pacing == "new":
    pacing = "RVV"
    
  if resultset == "all":
    resultset = "pseudoECG"
    
  output_dir_root=""
  if options.pipeline == "Forward":
    output_dir_root=Forward_dir_root
  elif options.pipeline == "Inverse":
    output_dir_root=Inverse_dir_root
  else:
    raise ValueError("unrecognized pipeline  "+ options.pipeline)
    

  output_dir = os.path.join(output_dir_root,data_dir)
  solutions_dir = os.path.join(output_dir,"solutions")

  samples_file = os.path.join(output_dir, "UQ_heart_forward_samples.mat")
    
    
  shape_data_dir = os.path.join(shapedata_dir_root,data_dir)
  shape_data_outdir = os.path.join(output_dir,"shape_models/")

  pacing_map = {"sinus" : "6105829a_01.120avg-pf.mat",
                "LV" : "61058472_03.120avg-pf.mat",
                "apex" : "6105b06d_10.120avg-pf.mat",
                "RV": "6105b134_16.120avg-pf.mat",
                "septal": "6105d439_43.120avg-pf.mat",
                "RVV" : "", "LVV" : "", "Rvot" : "", "Lvot" : "", "RVB" : "" }

  resultdir_map = {"pseudoECG": "Solutions_inria_iso/inria_iso_PseudoECG",
                  "NGpseudoECG": "Solutions_inria_iso/inria_iso_NG_PseudoECG",
                  "inriaLAT" : "Solutions_inria_iso/inria_iso_hires_LAT_sampled" ,
                  "ECGsimLAT" : "Solutions_ECGsim/activation_times_sampled" ,
                  "ECGsimBSP" : "Solutions_ECGsim/BSP_signals",
                  "ECG" : "Solutions_ECGsim/triECG" ,
                  "inriaTMP" : "Solutions_inria_iso/inria_iso_hires_sampled",
                  "inriaReEntryECG" : "Solutions_inria_iso/inria_iso_reentry_PseudoECG",
                  "inriaReEntryTMP" : "Solutions_inria_iso/inria_iso_reentry_hires_sampled",
                  "inriaReEntryLAT" : "Solutions_inria_iso/inria_iso_reentry_hires_LAT_sampled",
                  "inriaReEntryDF" : "Solutions_inria_iso/inria_iso_reentry_hires_DF_sampled",
                  "inriaReEntryPM" : "Solutions_inria_iso/inria_iso_reentry_hires_PhaseMaps_sampled",
                  "inriaReEntryHT" : "Solutions_inria_iso/inria_iso_reentry_hires_HTi_sampled"
  }
  results_dir = os.path.join(solutions_dir,resultdir_map[resultset])

  #    torso_pots_fname = "6105d439_43.120avg-pf.mat"
  #torso_pots_fname = "6105829a_01.120avg-pf.mat" # sinus
  #    torso_pots_fname = "61058472_03.120avg-pf.mat" # paced
  #torso_pots_fname = "6105b06d_10.120avg-pf.mat" # pace (phrenic)
  #    torso_pots_fname = "6105b134_16.120avg-pf.mat" # RV paced

  torso_pots_fname = pacing_map[pacing]
  torso_pots = os.path.join(torso_pots_dir,torso_pots_fname)
  
  image_dir = os.path.join(output_dir, "images")
  tmp_dir = os.path.join(output_dir,"tmp")
  
  SR_output_file = ""
  if options.pipeline == "Forward":
    SR_output_file = os.path.join(tmp_dir, torso_pots_fname[:-4]+"_UQ_heart_forward_SR_solutions.txt")
  elif options.pipeline == "Inverse":
    SR_output_file = os.path.join(tmp_dir, torso_pots_fname[:-4]+"_UQ_heart_inverse_SR_solutions.txt")
  else:
    raise ValueError("unrecognized pipeline  "+ options.pipeline)
  
  
  pipePaths = { "shape_data_dir" : shape_data_dir,
                "shape_data_outdir" : shape_data_outdir,
                "output_dir" : output_dir,
                "tmp_dir" : tmp_dir,
                "SR_output_file" : SR_output_file,
                "samples_file" : samples_file,
                "torso_pots_fname" : torso_pots_fname,
                "torso_pots": torso_pots,
                "solutions_dir" : solutions_dir,
                "results_dir" : results_dir,
                "image_dir" : image_dir
                }

  if not torso_pots_fname:
    print(pacing+" not mapped to recorded ECGs")
    return pipePaths
    

  

  if not os.path.exists(shape_data_outdir):
    os.makedirs(shape_data_outdir)



  if not os.path.exists(output_dir):
    os.makedirs(output_dir)

  if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)
    
  if not os.path.exists(image_dir):
    os.makedirs(image_dir)
    
  if not os.path.exists(output_dir+"UQ_data"):
    os.makedirs(output_dir+"UQ_data")
    
    
  


  # get matrix size
  tp = scipy.io.loadmat(torso_pots)
  tp_sz = tp['ts']['potvals'][0][0].shape
  #    N = 1024*tp_sz[1]
  #
  #    pot_size = (1024, tp_sz[1])

  return pipePaths
  
  

def get_params_from_scores(params, pipePaths):

  dim = params["dimension"]

  csv_files = glob.glob(pipePaths["shape_data_dir"]+"/*.csv" )

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
  #    print(type(data))
  #    print(type(data[0]))

  scores = np.array([ [float(dcr) for dcr in dr[1:] ] for dr in data[1:]])


  #    print(scores)
  std_scores = np.std(scores[:,0:dim],axis = 0)
  print("variance")
  print(std_scores*std_scores)

  mean_scores = np.mean(scores[:,0:dim],axis = 0)
  print("mean")
  print(mean_scores)
  cov_matrix = np.cov(scores[:,0:dim].T)
  print("covariance matrix")
  print(cov_matrix.shape)
  print(cov_matrix)

  range = np.vstack((np.max( scores[:,0:dim], axis = 0), np.min( scores[:,0:dim], axis = 0)))

  print(range)
    
  params["range"] = range
  params["mean"] =  mean_scores
  params["std" ] = std_scores
  params["covariance"] = cov_matrix

  return params
  
def make_samples(Filename, pce):
    # TODO needs to be fixed with different distribs
#    pce.generate_samples(oversampling = 20)
#    scipy.io.savemat(Filename, dict(samples=pce.samples))
  return False
  
  
def set_distribution(params,  **kwargs):

  switcher = {
      "gaussian": set_NormalDistribution,
      "uniform" : set_UniformDistribution
  }
  distrib = params["distribution"]
  print(distrib)

  func = switcher.get(distrib, lambda: "Invalid method")
  dist = func(params,  **kwargs)
  
  print("distribution made")
  
  order = params["order"]
  dim = params["dimension"]
  indices = TotalDegreeSet(dim=dim, order=order)
  
  print("number of indices = ",indices.get_indices().shape)
  pce = PolynomialChaosExpansion(indices, dist)
  
  print("pce initialized")
  return pce
  
def set_NormalDistribution(params,  **kwargs):
  print("making normal distribution")
  
  dimension = params["dimension"]
  cov = params["covariance"]
  mean = params["mean"]
  
  return NormalDistribution(cov=cov, mean=mean, dim=dimension)

def set_UniformDistribution(params,  **kwargs):
  print("making uniform distribution")
  dimension = params["dimension"]
  domain = params["domain"]

  dist = BetaDistribution(alpha=1, beta=1, dim=dimension, domain=domain)

  return dist
  
  
def make_meshes(indicator_files):
    # better done on the server
    # replace with script
  B = 0.01
  F= 2
  ofiles = []
  for files in indicator_files:
    tmp_file = os.path.split(files[0])
    ofile_rt = tmp_file[1][:-4] + "_MD_F" + str(F)
    call_text = cleaver_call + "-I -i " + files+" -j -B "+ str(B) + " -f matlab -o " + tmp_file[0] + " -n " + ofile_rt + " -v -F " + str(F)
    
    os.system(call_text)
    
    ofiles.append(os.path.join(tmp_file[0],ofile_rt+".mat"))
    
  return ofiles
        
    
def make_shape_points(samples):
    
  mean_shape_file = "mode_0-9.pts"
  eig_vect_files = ["eigenvectors0.eval", "eigenvectors1.eval", "eigenvectors2.eval", "eigenvectors3.eval", "eigenvectors4.eval"]

  eigvect = []
  for ef in eig_vect_files:
    eigvect.append(np.loadtxt(shape_data_dir+ef))

  mean_points = np.loadtxt(shape_data_dir+mean_shape_file)
  filenames=[]
  for l in range(samples.shape[0]):
  #        print(samples[l])
    displacement = np.zeros(mean_points.shape)
    for k in range(samples.shape[1]):
      displacement+=samples[l,k]*eigvect[k]
        
    shape_points = mean_points+displacement
    
    filename = shape_data_outdir+"model_params_"+f"{l:03}"+".pts"
    np.savetxt(filename,shape_points)
    filenames.append(filename)
    
  return filenames
            
def make_indicator_functions(shape_points_files):
    # run seg3D to make indicator functions
  seg3d_script = "/Users/jess/CIBC/FP/segmentation_error/seg_error/bash_scripts/make_seg3d_ind_function_MD.sh"
  
  outfiles = []
  for f in shape_points_files:
      os.system(seg3d_script+" "+f)
      outfiles.append([f[:-4]+"_surface_vent_vol.nrrd",f[:-4]+"_surface_lv_vol.nrrd",f[:-4]+"_surface_rv_vol.nrrd",f[:-4]+"_surface_air_vol.nrrd"])
  
  return outfiles
  
  
def make_runscript_forward():
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
  
  
def load_solutions(options, **kwargs):
  defaultKwargs = { "varname_def" : [],
  }
  kwargs = { **defaultKwargs, **kwargs }
    
  solutions = load_mat_data_from_dir(options, **kwargs)
  
  return solutions
    
def load_mat_data_from_dir(options, **kwargs):

  print(kwargs["varname_def"])
  print(kwargs.keys())
  
  results_dir = kwargs["results_dir"]

  if options.pacing == "default":
    raise ValueError("no solutions for "+options.pacing+" pacing")
  files = os.listdir(results_dir)
  files.sort()
  
  

  print(options.resultset)
  print(options.pacing+"_stim")
  print(len(files))
    
  ofiles = []
  solutions = []
  data=np.array([])

  fid = fid_map[options.pacing]
  plot_fids = False
  if options.make_FIDs or (not fid and not ("LAT" in options.resultset or "ReEntry" in options.resultset)):
  #        raise ValueError(pacing+" not fully implemented for ECGs. Run --make_FIDs option")
    # only QRS implemented right now
    fid = {"sf": 1, "qon": 0, "qoff": - 1, "ton": 0, "toff": - 1}
    plot_fids = True

  for f in files:
  #        print(f)
    if options.pacing+"_stim" in f:
  #            print(f)
      temp_mat = loadMatlabMatrix(os.path.join(results_dir,f), kwargs["varname_def"])

      ofiles.append(f)

      #            qon = int(fid["qon"]/fid["sf"])
      #            toff = int(fid["toff"]/fid["sf"])

      #            print(temp_mat.shape)

      if 'pseudoECG' == options.resultset:
        qon = int(fid["qon"]/fid["sf"])
        toff = int(fid["toff"]/fid["sf"])
        data = np.transpose(temp_mat)[:, qon:toff]
        solutions.append(data.flatten()/(4*np.pi))
      elif 'NGpseudoECG' == options.resultset:
        qon = int(fid["qon"]/fid["sf"])
        toff = int(fid["toff"]/fid["sf"])
        data = temp_mat[:, qon:toff]
      #                solutions.append(data.flatten()/(4*np.pi))
        solutions.append(data.flatten())
      elif 'ECGsimBSP' == options.resultset or 'ECG' == options.resultset:
      #                data = temp_mat
        qon = int(fid["qon"]/fid["sf"])
        toff = int(fid["toff"]/fid["sf"])
        temp_mat
        data = temp_mat[:, qon:toff]
        solutions.append(data.flatten())
      elif 'ECGsimLAT' == options.resultset or 'inriaTMP' == options.resultset:
        data = temp_mat
        solutions.append(data.flatten())
      elif 'inriaLAT' == options.resultset:
        data = temp_mat
        solutions.append(np.transpose(data).flatten())
      elif "inriaReEntryECG" == options.resultset:
        qon = 0
        toff = -1
        data = np.transpose(temp_mat)[:, qon:toff]
        solutions.append(data.flatten()/(4*np.pi))
      elif "inriaReEntryPM" == options.resultset or "inriaReEntryHT" == options.resultset or "inriaReEntryTMP" == options.resultset:
        qon = 0
        toff = -1
        data = temp_mat[:, qon:toff]
        solutions.append(data.flatten())
      elif "inriaReEntryLAT" == options.resultset or "inriaReEntryDF" == options.resultset:
        data = temp_mat
#        solutions.append(np.transpose(data).flatten())
        solutions.append(data.flatten())
      else:
        print("solution file could not be loaded ", os.path.join(results_dir,f))
        return False
                
  if plot_fids:
    choose_FIDS(solutions, data.shape)
    plt.show()
    return False
    
  if not solutions:
    print("files matching "+options.pacing+" not found in "+results_dir)
    return False
    
                
#    if 'ECGsim' in results_dir or 'inriaTMP' == options.resultset:
#        sol_size = np.shape(temp_mat)
#    else:
#        sol_size = np.shape(np.transpose(temp_mat))
            
#    print(len(solutions))
#    print(np.vstack(solutions).shape)
#    print(ofiles)
    
  return solutions, data.shape, ofiles
  
def compute_surface_variance(std_devs, pipePaths):
  eig_vect_files = ["eigenvectors0.eval", "eigenvectors1.eval", "eigenvectors2.eval", "eigenvectors3.eval", "eigenvectors4.eval"]
  print("computing surface variance")
  magns = []

  for k in range(len(std_devs)):
    evec = np.loadtxt(os.path.join(pipePaths["shape_data_dir"], eig_vect_files[k]))*std_devs[k]
    
    mag = np.sqrt(np.sum(evec*evec, axis=1))
    
    print(mag.shape)
    magns.append(mag)
    
  magnitude = np.vstack(magns)
  print(magnitude.shape)

  var =np.sum(magnitude*magnitude,axis = 0)

  print(var.shape)
  file = os.path.join(pipePaths["output_dir"], "surface_variance.mat")
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

def find_LATs(model_output, model_size, **kwargs):
  defaultKwargs = { "sample_rate" : 1 }
  kwargs = { **defaultKwargs, **kwargs }

  lats = []
  dys = []
  for pots in model_output:
    pots = np.array(pots)
    signals = np.resize(pots,model_size)
    [lat,dy]=findLAT(signals, **kwargs )
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
  
def plot_ecg_uq(time, median, quantiles, ecgs, ecg_size, filename="", num_samples= 10, title = ""):

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

  fig, ax = plt.subplots(6,1, sharex = True, sharey = True,
                        gridspec_kw={"wspace": 0.08, "hspace": 0.1,
                        "left" : 0.2, "bottom" : 0.1, "right" : 0.95, "top" : 0.95} )
  fig.set_figwidth(4)
  fig.set_figheight(6)

  for l in range(len(lead_nums)):
    ax[l].plot(time, median[lead_nums[l],:], '-k',label='median')
    ax[l].set(ylabel='v'+str(l)+'(mV)')
    ax[l].tick_params(direction='in')
    ax[l].spines[['right', 'top']].set_visible(False)
  #        if l>0:
  #            ax[l].sharex(ax[0])
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

  if len(filename)>0:
    plt.savefig(filename)
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
  
def reentry_processing(*args, **kwargs):
  # arg order
  args_str = ["model_output", "model_size", "sample_params", "options", "output_dir", "samples", "ofiles"]
  
  if not len(args) == len(args_str):
    raise ValueError("need "+len(args_str)+" args: "+ args_str)

  switcher = {"inriaReEntryECG" : reentryECG,
              "inriaReEntryTMP" : reentryTMP,
              "inriaReEntryLAT" : reentryLAT,
              "inriaReEntryDF"  : reentryDF,
              "inriaReEntryPM"  : reentryPM,
              "inriaReEntryHT"  : reentryHT
  }
  
  
  resultset = args[3].resultset
  
  print(kwargs["output_dir"])
  UQ_file = make_UQ_filenames(args[4], args[2], args[3], tag = "UQ_values")
  print(UQ_file)
  RMS_UQ_file = make_UQ_filenames(args[4], args[2], args[3], tag = "RMS_UQ_values")
  
  func = switcher.get(resultset, lambda: "Invalid method")
  some_outputs = func(*args,  **kwargs)
  
  return
  
def reentryECG(*args, **kwargs):
  print("TODO")
  return
  
def reentryTMP(*args, **kwargs):
  print("TODO")
  return
  
def reentryLAT(*args, **kwargs):

  pce = set_distribution(args[2])
  
  if not check_Model_outputs(args[5], args[0], pce, args[6], args[3]):
    return
  
  mod_out = []
  print(args[1])
  print(len(args[0]))
  
  for mo in args[0]:
    n_pks = int(len(mo)/args[1][0])
    print(len(mo))
    print(n_pks)
    mod_out.append(cycle_lengths(mo.reshape(1024, n_pks)))
    
  mod_size = (args[1][0], 1)
    
    
    
  UQ_file = make_UQ_filenames(args[4], args[2], args[3], tag = "LAT_CL_UQ_values")
  

  Q = 4  # Number of quantile bands to plot

  # only two dimensions currently
  N = np.prod(mod_size)
  print(N)
  print(mod_size)
    
    
  if not os.path.exists(UQ_file):
    if args[3].hardchecks and not args[3].run_pce:
      print(" cannot run pce implicitly with hardchecks on, or in batch mode.  Skipping "+args[3].pacing)
      return
    args[3].run_pce = True
    
    
  if args[3].run_pce:
   
    pce.set_samples(args[5])
    print(np.vstack(mod_out).shape)
    print(len(mod_out))
    run_pce_data(pce, np.vstack(mod_out), mod_size, args[2], args[4], UQ_file, Q)
    print("PCE model built")
  
  cycles_UQ = load_pce_data(UQ_file)
  print(" pce stats loaded from disk")

  summary_stats = make_summary_stats(cycles_UQ, len(mod_out))
  print("cycle stats")
  #    print(summary_stats)
  print(tabulate(summary_stats, headers = "firstrow"))

  if args[3].save:
    save_stats(UQ_file, summary_stats)
  return
  
def reentryDF(*args, **kwargs):
  print("TODO")
  return
  
def reentryPM(*args, **kwargs):
  
  sol_input = load_solutions(args[3], varname_def = "cycles")
  
  print(len(sol_input[0]))
  
#  for sol in sol_input[0]:
#    print(len(sol))
    
  if sol_input:
    model_output = sol_input[0]
    model_size = sol_input[1]
    ofiles = sol_input[2]
  else:
    print("cannot load "+pacing+" data.  skipping ")
    return

  pce = set_distribution(args[2])
  
  if not check_Model_outputs(args[5], model_output, pce, ofiles, args[3]):
    return
  
  mod_out = []
  for mo in model_output:
    n_pks = int(len(mo)/model_size[0])
    mod_out.append(cycle_lengths(mo.reshape(1024, n_pks)))
    
  mod_size = (model_size[0], 1)
    
    
    
  UQ_file = make_UQ_filenames(args[4], args[2], args[3], tag = "cycle_UQ_values")
  

  Q = 4  # Number of quantile bands to plot

  # only two dimensions currently
  N = np.prod(mod_size)
  print(N)
  print(mod_size)
    
    
  if not os.path.exists(UQ_file):
    if args[3].hardchecks and not args[3].run_pce:
      print(" cannot run pce implicitly with hardchecks on, or in batch mode.  Skipping "+args[3].pacing)
      return
    args[3].run_pce = True
    
    
  if args[3].run_pce:
   
    pce.set_samples(args[5])
    print(np.vstack(mod_out).shape)
    print(len(mod_out))
    run_pce_data(pce, np.vstack(mod_out), mod_size, args[2], args[4], UQ_file, Q)
    print("PCE model built")
  
  cycles_UQ = load_pce_data(UQ_file)
  print(" pce stats loaded from disk")

  summary_stats = make_summary_stats(cycles_UQ, len(mod_out))
  print("cycle stats")
  #    print(summary_stats)
  print(tabulate(summary_stats, headers = "firstrow"))

  if args[3].save:
    save_stats(UQ_file, summary_stats)
  
  return
  
def reentryHT(*args, **kwargs):
  print("TODO")
  return

    
    
def run_all_results(result_lists, sample_params, options, **kwargs):

  pacing = options.pacing

  for r in result_lists:
    print("running :", r)
    print("===============")
    options.resultset = r

    if pacing == "all" or pacing == "new" or pacing == "old":
      options.compute_samples=False
      options.run_model=False
      options.make_FIDs=False
      options.hardchecks=True
      
      for p in pacing_lists[pacing]:
        options.pacing = p
        run_pipeline(sample_params, options, **kwargs)
        
    else:
      run_pipeline(sample_params, options, **kwargs)


  return
    

def run_pipeline(sample_params, options, **kwargs):

  pacing = options.pacing
  resultset = options.resultset
  distrib = options.distrib

  print("running :", pacing)
  print("===============")

  pipePaths = set_dirs(options)
  
  print(type(pipePaths))
  
  samples_file = pipePaths["samples_file"]

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
            
  #        make_samples(samples_file, pce)
  #        print("Samples compute.  Make sure to run the model")
    print("cannot compute and save samples at this time")
    return False
  #        return
    samples = pce.samples
    options.run_model = True
  else:
    tmp = scipy.io.loadmat(samples_file)
    samples = tmp["samples"]
    
  if options.run_model:
    # Likely out of date.  Little need to run this once done
    files = make_shape_points(samples)
    ind_func_files = make_indicator_functions(files)
    print("make cleaver meshes")
    
    return
    
    
  if "pce" in locals():
    raise ValueError("PCE not cleared")
    
  pce = set_distribution(sample_params)

  fid = fid_map[options.pacing]

  if not fid and options.hardchecks:
    raise ValueError(" FID missing and cannot be run with hardchecks on or in batch mode ")

  if not fid and not ("LAT" in options.resultset or "ReEntry" in options.resultset):
    options.make_FIDs=True


  sol_input = load_solutions(options, **kwargs)

  if "model_output" in locals():
    raise ValueError("model_output not cleared")

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
  
  if not check_Model_outputs(samples, model_output, pce, ofiles, options):
    return
  
  #    elif  (not np.shape(samples)[0]==len(model_output)) and options.hardchecks:
  #        raise ValueError("sample and model mismatch that cannot be resolved in batch mode or with hardchecks")

  #    samples = sample_workaround(samples, domain)

  output_dir = kwargs["output_dir"]
  image_dir = kwargs["image_dir"]
  
  
  
  print(output_dir)
  UQ_file = make_UQ_filenames(output_dir, sample_params, options, tag = "UQ_values")
  print(UQ_file)
  RMS_UQ_file = make_UQ_filenames(output_dir, sample_params, options, tag = "RMS_UQ_values")

  Q = 4  # Number of quantile bands to plot

  # only two dimensions currently
  N = np.prod(model_size)
  print(N)
  print(model_size)
  do_RMS=not np.any(N==model_size)

  if not os.path.exists(UQ_file):
    if options.hardchecks and not options.run_pce:
      print(" cannot run pce implicitly with hardchecks on, or in batch mode.  Skipping "+pacing)
      return
    options.run_pce = True
    
  if not "ReEntryLAT" in resultset:
    
    if options.run_pce:
     
      pce.set_samples(samples)
      
      print(np.vstack(model_output).shape)
      print(len(model_output))
      
      run_pce_data(pce, np.vstack(model_output), model_size, sample_params, output_dir, UQ_file, Q)
      print("PCE model built")
      
    if do_RMS:
      print("running RMS")
      RMS = model_RMS(model_output, model_size)
      
      
      print(RMS.shape)
    #        print(RMS)
    
      if not os.path.exists(RMS_UQ_file):
        if options.hardchecks and not options.run_pce:
          print(" cannot run pce implicitly with hardchecks on, or in batch mode.  Skipping "+pacing)
          return
        options.run_pce = True
          
      if options.run_pce:
        print("running RMS PCE")
        if "pce_rms" in locals():
          raise ValueError("PCE_rms not cleared")
        pce_rms = set_distribution(sample_params)
        pce_rms.set_samples(samples)
        run_pce_data(pce_rms, RMS, (1,model_size[0]),sample_params, output_dir, RMS_UQ_file, Q)
      
    pot_UQ = load_pce_data(UQ_file)

    print(" pce stats loaded from disk")

    summary_stats = make_summary_stats(pot_UQ, len(model_output))
    print("pot stats")
    #    print(summary_stats)
    print(tabulate(summary_stats, headers = "firstrow"))

    if options.save:
      save_stats(UQ_file, summary_stats)
      
      
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
      ECG_plot_file = os.path.join(image_dir,  resultset+"_"+pacing+"_"+distrib+"_UQ_ECG_plot.png")
      plot_ecg_uq(time, pot_UQ["median"], quantiles, model_output, model_size, ECG_plot_file, 0,  pacing)

    #    if options.plot_channels>=0:
    #        plot_all_channels(model_output, LATs, dys, pot_size, options.plot_channels)
       
    #    if options.plot_solutions:
    #        channels = (stdev_LAT>3*np.mean(stdev_LAT)).nonzero()[0]
    #        plot_all_solutions(model_output, LATs, dys*50, pot_size, channels)
          
      print("rms_UQ")
      print(RMS.shape)
      print(rms_UQ["median"].shape)
      RMS_plot_file = os.path.join(image_dir,  resultset+"_"+pacing+"_"+distrib+"_UQ_RMS_plot.png")
      plot_quantile_bands(time, rms_UQ["median"], rms_UQ["quantiles"], Q, RMS_plot_file,  title =  pacing, ylabel = "RMS (mV)", xlabel = "time (ms)")
      
  if "ReEntry" in resultset:
    reentry_processing(model_output, model_size, sample_params, options, output_dir, samples, ofiles, **kwargs)
    
  return
  


################################
###############################

    

def main():

    parser = build_parser()
    options = parser.parse_args()
    
    
    
    return
    

if __name__ == "__main__":
    main()
