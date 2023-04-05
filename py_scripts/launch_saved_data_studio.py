import os
import platform
import argparse
import subprocess
import sys
import glob
from shutil import copyfile

# default values
data_dir = "vent_centered_sub"
numpoints = 1024


if len(sys.argv) == 1:
    print("this script normally expects the in put directory and number of points (optional)")
    print("using default path: " + data_dir)
    print("using defaul number of points: " + str(numpoints))
elif len(sys.argv) == 2:
    data_dir = sys.argv[1]
    print("using defaul number of points: " + str(numpoints))
elif len(sys.argv) == 3:
    data_dir = sys.argv[1]
    numpoints = sys.argv[2]

shapeworks_root = "/Users/jess/software/ShapeWorks/"
sw_pypath = os.path.join(shapeworks_root,"Examples/Python/")

sys.path.append(sw_pypath)
print(sw_pypath)

import setupenv

shapeworks_bin_path = os.path.join(shapeworks_root,"build/bin/")

setupenv.setup_shapeworks_env(shapeworks_bin_dir=shapeworks_bin_path, verbose=False)

from CommonUtils import *

from GroomUtils import *
from OptimizeUtils import *
from AnalyzeUtils import *



datasetName = os.path.join(sw_pypath,data_dir)
parentDir = os.path.join(sw_pypath, "Test_" + data_dir + "_mesh")
print(sw_pypath)
print(parentDir)

meshFiles = sorted(glob.glob(os.path.join(datasetName, "meshes/*.vtk")))


prep_outdir = os.path.join(parentDir,"PrepOutput/distance_transforms")

particleFolder = os.path.join(parentDir,"PointFiles/"+str(numpoints))

print(particleFolder)

localPointFiles = sorted(glob.glob( os.path.join(particleFolder, "*local.particles" )))
worldPointFiles = sorted(glob.glob( os.path.join(particleFolder, "*world.particles" )))

print(meshFiles)
print(localPointFiles)
print(worldPointFiles)

launchShapeWorksStudio(os.path.join(parentDir,"PointFiles/"), meshFiles, localPointFiles, worldPointFiles)

