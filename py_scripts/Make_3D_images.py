import scipy.io
import os
import argparse
import socket

import numpy as np

import sys


fid_map = {
            "apex" : {"sf": 0.5, "qon" : 0, "qoff" : 54, "stoff" : 110, "toff": 220},
            "LV" : {"sf": 0.5, "qon" : 0, "qoff" : 65, "stoff" : 114, "toff": 230},
            "RV" : {"sf": 0.5, "qon" : 0, "qoff" : 66, "stoff" : 116, "toff": 230},
            "septal" : {"sf": 0.5, "qon" : 0, "qoff" : 50, "stoff" : 108, "toff": 220},
            "RVV" : {"sf": 0.5, "qon" : 0, "qoff" : 60, "stoff" : 110, "toff": 225},
            "LVV" : {"sf": 0.5, "qon" : 0, "qoff" : 60, "stoff" : 115, "toff": 220},
            "Rvot" : {"sf": 0.5, "qon" : 0, "qoff" : 60, "stoff" : 105, "toff": 220},
            "Lvot" : {"sf": 0.5, "qon" : 0, "qoff" : 60, "stoff" : 105, "toff": 220},
            "RVB" : {"sf": 0.5, "qon" : 2, "qoff" : 65, "stoff" : 105, "toff": 225},
            "sinus" : {"sf": 0.5, "qon" : 0, "qoff" : 37, "stoff" : 100, "toff": 200}
}

pacing_lists = {"all" : [ k for  k in fid_map.keys()],
                "old" : ["sinus", "apex", "LV", "RV", "septal"],
                "new" : ["RVV", "LVV", "Rvot", "Lvot", "RVB"]
}

#pacing_lists = {"all" : ["RV"]
#}

resultset_lists = { "all" : [ "NGpseudoECG", "inriaLAT", "ECGsimLAT", "ECGsimBSP", "ECG"]}


#distrib_lists = {"all" : [ "uniform", "gaussian"]}
distrib_lists = {"all" : [ "uniform"]}

environ_vars = { "data_dir" : "DATADIR", "resultset_h" : "RESULTSET_HEART", "resultset_t" : "RESULTSET_TORSO", "pacing": "PACING", "distrib" : "DISTRIBUTION" }


scirun_call = "/Users/jess/software/SCIRun_testing2/bin_515/SCIRun/SCIRun_test"
scirun_net = "/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/view_UQ_data.srn5"
  
resultset_pairs = [  ("ECGsimLAT", "ECGsimBSP"), ("inriaLAT", "NGpseudoECG") ]
#resultset_pairs = [   ("inriaLAT", "NGpseudoECG") ]


data_dirs = ["/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/"]


for dd in data_dirs:
  os.environ["DATADIR"] = dd
  
  for dist in distrib_lists["all"]:
    os.environ["DISTRIBUTION"] = dist
    
    for rp in resultset_pairs:
      os.environ["RESULTSET_HEART"] = rp[0]
      os.environ["RESULTSET_TORSO"] = rp[1]
      
      for p in pacing_lists["all"]:
        if p == "sinus" and rp[1] == "NGpseudoECG":
          continue
          print("swapping NPpseudoECG for pseudoECG for sinus pacing")
          os.environ["RESULTSET_TORSO"] = "pseudoECG"
          
        os.environ["PACING"] = p

        call = scirun_call+" -0 -E "+scirun_net
        print(call)
        os.system(call)

