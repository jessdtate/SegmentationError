import threading
import os
import fnmatch
import sys

parent_dir = "/Users/jess/software/ShapeWorks/Examples/Python/Test_epi_centered/"
network = "/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/view_sw_corr_points.srn5"

# make sure to get the ISO.nrrd files
iso_dir=os.path.join(parent_dir,"PrepOutput/groom_and_meshes/")
# need to convert .particles to .txt. bash script available
point_dir=os.path.join(parent_dir,"PointFiles/512")

iso_files = fnmatch.filter(os.listdir(iso_dir), '*.ISO.nrrd')
point_files = fnmatch.filter(os.listdir(point_dir), '*.particles.txt')

if not len(iso_files)==len(point_files):
    print(iso_dir)
    print(iso_files)
    print(point_dir)
    print(point_files)
    raise ValueError("something wrong with files in directories, they don't match")

print('list of files have same lengths')

#scirun_load_network(network)
#print("network loaded")

for point_file in point_files:
    print(point_file)
    p_pnt = point_file.find('.')
    
    iso_file =fnmatch.filter(os.listdir(iso_dir), point_file[0:p_pnt]+'*.ISO.nrrd')
    print(iso_file)
    
    scirun_set_module_state('ReadField:1','Filename',os.path.join(point_dir, point_file))
    scirun_set_module_state('ReadField:0','Filename',os.path.join(iso_dir, iso_file[0]))
    scirun_execute_all()
