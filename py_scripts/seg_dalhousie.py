sys.path.append("/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/")
#path to the scan
filename  = "/Users/jess/FP/segmentation_error/Dalhousie_seg/finished_segs/vent_staple.nrrd"
# output directory
outputpath = "/Users/jess/FP/segmentation_error/Dalhousie_seg/finished_segs/"
# load script in python

exec(open("/Users/jess/FP/segmentation_error/seg_error/py_scripts/upsample_segs.py").read())

upsample_heart(filename, outputpath, "vent_staple_hires")

filename  = "/Users/jess/FP/segmentation_error/Dalhousie_seg/finished_segs/torso_staple.nrrd"

torso_points(filename,outputpath,"torso_points")




