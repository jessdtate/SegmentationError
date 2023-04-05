import scipy.io
import sys
import os
import fnmatch


path = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/vent_MD_cheat/9samples/Models_srf_reduction/"


files = fnmatch.filter(os.listdir(path), 'model_params_*_reduTriMesh.mat')

for name in files:

#name = "model_params_010_reduTriMesh.mat"
    
    fname = os.path.join(path, name)
    outname = name[:-4]+'_fixed.mat'

    print(name)
    print(outname)

    tmp = scipy.io.loadmat(fname)

    print(tmp)



    field={'node' : tmp['reduNode'], 'face' : tmp['reduFace']}

    scipy.io.savemat(os.path.join(path, outname),{'field': field})

