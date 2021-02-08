import os

# exec(open("/Users/jess/CIBC/FP/segmentation_error/seg_error/py_scripts/implicit_model.py").read())

sys.path.append("/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/")
#path to the volumes
vol_path = "/Users/jess/software/ShapeWorks/Examples/Python/Test_epi_centered/PrepOutput/groom_and_meshes"

vol_fname = os.path.join(vol_path, "BO_epi_centered.isores.center.pad.cropped.ISO.nrrd")

point_path = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/epi_centered/9samples/"

outputpath = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/epi_centered/9samples/SR_surfaces/"


name = "epi_9_0-10.pts"
point_fname = os.path.join(point_path, name)

# load script in python
p_name = name.find('.')

out_fname =  name[0:p_name]+"_surface"

exec(open("/Users/jess/CIBC/FP/segmentation_error/seg_error/py_scripts/implicit_function_From_points.py").read())

implicit_model(vol_fname, point_fname, outputpath, out_fname)






