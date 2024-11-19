import numpy as np
import socket
import sys

hostname = socket.gethostname()
print(hostname)


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

resultset_lists = { "all" : [ "pseudoECG", "NGpseudoECG", "inriaLAT", "ECGsimLAT", "ECGsimBSP", "ECG"],
                    "reentry" : ["inriaReEntryECG", "inriaReEntryDF", "inriaReEntryPM"]
}

#domain = np.array([[-165, 165], [-112, 112], [-85, 85], [-55, 55], [-45,45]]).T
# 1.5 sigma
default_domain = np.array([[-101, 101], [-55, 55], [-40, 40], [-35, 35], [-25,25]]).T
default_cov = np.array( [[ 1.0e+04, -3.5e-05, -5.9e-04, -6.1e-04, 4.3e-04],
                         [-3.5e-05,  3.0e+03, -1.0e-04, -1.6e-03, 8.6e-04],
                         [-5.9e-04, -1.0e-04,  1.7e+03,  2.1e-04, -3.4e-04],
                         [-6.1e-04, -1.6e-03,  2.1e-04,  1.2e+03, -2.1e-05],
                         [ 4.3e-04,  8.6e-04, -3.4e-04, -2.1e-05, 6.2e+02]])
default_mean = np.array( [0,0,0,0,0])
                         
# parameter set for distribution options
default_sample_params = { "order" : 5, "dimension": 5, "alpha": 1, "beta": 1}

default_datadir = "vent_MD_cheat/11samples/"




if "metis" in hostname or "local" in hostname:
#if True:
  sys.path.append("/Users/jess/CIBC/FP/UQ/py_scripts")
  Forward_dir_root = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/"
  Inverse_dir_root = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Inverse/"
  shapedata_dir_root = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/"


  torso_pots_dir = "/Users/jess/CIBC/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/Interventions/BSPM_pf/"


  scirun_call = "/Users/jess/software/SCIRun_testing2/bin_515/SCIRun/SCIRun_test"
  #scirun_call = "/Applications/SCIRun.app/Contents/MacOS/SCIRun"
  seg3D_call = "/Applications/Seg3D2.app/Contents/MacOS/Seg3D2"
  cleaver_call = "/Users/jess/software/cleaver2/build/bin/cleaver-cli"
    
elif "cibc" in hostname or "shell" in hostname:

  sys.path.append("/home/sci/jess/UQ/py_scripts")
  Forward_dir_root = "/home/sci/segmentation_error/Dalhousie_seg/UQ_data/Forward/"
  Inverse_dir_root = "/home/sci/segmentation_error/Dalhousie_seg/UQ_data/Inverse/"
  shapedata_dir_root = "/home/sci/segmentation_error/Dalhousie_seg/shape_data/"


  torso_pots_dir = "/usr/sci/projects/Edgar/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/Interventions/BSPM_pf/"
    
    
    
    
