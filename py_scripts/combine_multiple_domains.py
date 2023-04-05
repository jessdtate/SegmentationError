import os
import platform
import argparse
import subprocess
import sys
import glob
from shutil import copyfile


shapeworks_root = "/Users/jess/software/ShapeWorks/"
sw_pypath = shapeworks_root+"Examples/Python"

sys.path.append(sw_pypath)


import setupenv

from CommonUtils import *


shapeworks_bin_path = shapeworks_root+"build/bin/"

setupenv.setup_shapeworks_env(shapeworks_bin_dir=shapeworks_bin_path, verbose=False)

from vent_multidomain import concatenate_particle_files

datasetName = "vent_multidomain"

d_dataset = ["Test_epi_centered","Test_lv_centered","Test_rv_centered"]
d_tag = ["epi","LV_endo","RV_endo"]
d_size = [512,256,256]






parentDir = sw_pypath+"/Test_vent_multidomain_cheat/"



def copy_data():

    prep_outdir = parentDir+"PrepOutput/distance_transforms"
    
    if not os.path.exists(prep_outdir):
        os.makedirs(prep_outdir)
        
    point_outdir = parentDir+"PointFiles/256_256_512/"
        
    if not os.path.exists(point_outdir):
        os.makedirs(point_outdir)
    
    for k in range(3):
    
        outdir = parentDir+"PrepOutput/distance_transforms"
        
        print(k)
        files_d = glob.glob(sw_pypath+ "/"+ d_dataset[k] + "/PrepOutput/distance_transforms/*" + d_tag[k]+"*nrrd")
        print(files_d)

        
        for f in files_d:
            f_new = os.path.basename(f).replace("_"+d_tag[k],"_d"+str(k+1))
            print(f_new)
            copyfile(f, os.path.join(prep_outdir,f_new))
        
        files_d = glob.glob(sw_pypath+ "/"+ d_dataset[k] + "/PointFiles/" + str(d_size[k]) + "/*" + d_tag[k]+"*.particles")
        print(files_d)
        
        for f in files_d:
            f_new = os.path.basename(f).replace("_"+d_tag[k],"_d"+str(k+1))
            print(f_new)
            copyfile(f, os.path.join(point_outdir,f_new))
            
    return True


copy_data()

prepDir = parentDir + "PrepOutput/"
dtFiles = sorted(glob.glob( prepDir+"distance_transforms/*.nrrd" ))

pointDir = parentDir + "PointFiles/"

particleFolder = pointDir + "256_256_512/"

localPointFiles = sorted(glob.glob( particleFolder + "*local.particles" ))
worldPointFiles = sorted(glob.glob( particleFolder + "*world.particles" ))

file_type = "world"
concatenate_particle_files(file_type, 3, particleFolder, prepDir + "distance_transforms", 1, pointDir)

command = "ShapeWorksStudio "+ pointDir+ "multiple_domain_"+file_type+".xml"
os.system(command)
