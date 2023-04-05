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
            print("correct_image.py -i <inputdir> -o <outputdir>")
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
        outputdir = inputdir
    
    files = sorted(glob.glob(os.path.join(inputdir, "*.mat")))
    
    if len(files) == 0:
        print(" input directory does not have any valid files (.mat)")
        
    if not os.path.exists(outputdir):
        print("making dir: " + outputdir)
        os.makedirs(outputdir)
    
    for f in files:
        tmp = scipy.io.loadmat(f)
        out_filename_front = os.path.join(os.path.dirname(f), os.path.basename(f) + "_front.png")
        out_filename_back = os.path.join(os.path.dirname(f), os.path.basename(f) + "_back.png")
        
        print(tmp.keys())
        print(f)
        print(out_filename_front)
        print(out_filename_back)
        

if __name__ == "__main__":
   main(sys.argv[1:])


    

