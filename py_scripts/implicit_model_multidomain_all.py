import os
import fnmatch

# exec(open("/Users/jess/CIBC/FP/segmentation_error/seg_error/py_scripts/implicit_model_multidomain_all.py").read())

sys.path.append("/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/")
#path to the volumes
#vol_path = "/Users/jess/software/ShapeWorks/Examples/Python/Test_epi_centered/PrepOutput/groom_and_meshes"
#vol_fname = os.path.join(vol_path, "BO_epi_centered.isores.center.pad.cropped.ISO.nrrd")
#point_path = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/epi_centered/9samples/"
#outputpath = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/epi_centered/9samples/SR_surfaces/"


#vol_path = "/Users/jess/software/ShapeWorks/Examples/Python/Test_lv_centered/PrepOutput/groom_and_meshes"
#vol_fname = os.path.join(vol_path, "BO_lv_endo.isores.center.pad.cropped.ISO.nrrd")
#point_path = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/lv_centered/9samples/"
#outputpath = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/lv_centered/9samples/SR_surfaces/"

vol_path = "/Users/jess/software/ShapeWorks/Examples/Python/Test_vent_multidomain/PrepOutput/groom_and_meshes"
vol_fname = os.path.join(vol_path, "BO_d3.isores.center.pad.cropped.ISO.nrrd")
point_path = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/vent_MD_cheat/9samples/shape_modle"
outputpath = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/vent_MD_cheat/9samples/shape_model/"


exec(open("/Users/jess/CIBC/FP/segmentation_error/seg_error/py_scripts/implicit_function_From_points_multidomain.py").read())

files = fnmatch.filter(os.listdir(point_path), 'model_param*.pts')

for name in files:
    
    point_fname = os.path.join(point_path, name)

    # load script in python
    p_name = name.find('.')

    out_fname =  name[0:p_name]+"_surface"

    pnt = load_points(point_fname)
    
    print(name)

    implicit_model(vol_fname, pnt, outputpath, out_fname)






