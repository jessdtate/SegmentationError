#!/Users/jess/anaconda3/bin/python

import os
import threading
import numpy as np
import scipy.io
import sys, getopt
import glob


def main(argv):
    in_set = False
    out_set = False
    inputdir = ''
    outputdir = ''
    
    opts, args = getopt.getopt(argv,"hi:o:",["idir=","odir="])
    
    print("opts:")
    print(opts)
    print("args:")
    print(args)
    
#        print("correct_image.py -i <inputdir> -o <outputdir>")
#        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("check_samples.py -i <inputdir> -o <outputdir>")
            sys.exit()
        elif opt in ("-i", "--idir"):
            inputdir = arg
            in_set = True
        elif opt in ("-o", "--odir"):
            outputdir = arg
            out_set = True

    if not in_set:
        print("input directory required")
        sys.exit(2)
    if not out_set:
        outputdir = os.path.join(inputdir,'../')
        
    mean_filename = os.path.join(inputdir, '../../../../../shape_data/vent_MD_cheat/11samples/', 'mode_0-9.pts')
    
    output_filename = "UQ_samples_reverse.mat"
    
    num_axes = 5
    evals=[]
    for k in range(5):
        eval_name = os.path.join(inputdir, '../../../../../shape_data/vent_MD_cheat/11samples/', 'eigenvectors'+str(k)+'.eval')
        evals.append(np.loadtxt(eval_name))
        
    
    
    files = sorted(glob.glob(os.path.join(inputdir, "model_params*.pts")))
    
    if len(files) == 0:
        print(" input directory does not have any valid files ('model_params*.pts')")
        
    if not os.path.exists(outputdir):
        print("making dir: " + outputdir)
        os.makedirs(outputdir)
        
    params =[]
    mean_pts = np.loadtxt(mean_filename)
    
    for f in files:
        
        if  "surface" in f:
            continue
            
        print(f)
        param_pts = np.loadtxt(f)
        
        params.append(getParams(mean_pts, evals, param_pts))
        
    scipy.io.savemat(os.path.join(outputdir,output_filename),{'samples' : params})
        
    
def getParams(mean_pts, evals, param_pts):
    threshold = 1e-6
    diff = param_pts-mean_pts
    params=[]

    evst = []
    d = diff.flatten('F')
    for ev in evals:
        evst.append(np.hstack([ev[:,0],ev[:,1],ev[:,2]]))
    A = np.vstack(evst)
    ll = np.linalg.lstsq(A.T, d)
    if ll[1]>threshold:
        print('Warning: residual of :',11[1])
    
    return ll[0].tolist()

if __name__ == "__main__":
   main(sys.argv[1:])


    

