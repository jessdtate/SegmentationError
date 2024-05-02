import scipy.io
import os
import fnmatch
import numpy as np
from skimage.io import imsave
import glob

#image_dir = "/Users/jess/software/ShapeWorks/Examples/Python/Test_vent_multidomain/image"
#image_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/epi_centered/9samples/SR_images/front/"
#image_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/vent_centered_sub/SR_images/back/"
#image_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/vent_MD/11samples/SR_images/front/"
#image_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/vent_MD_cheat/9samples/SR_images/back"

image_dir = "/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/images"

#image_files = sorted(glob.glob(os.path.join(image_dir, "*.txt")))
image_files = sorted(glob.glob(os.path.join(image_dir, "*.mat")))

for im_f in image_files:
    print(im_f)
#    im_flat = []
#    im_flat = np.loadtxt(im_f)
    tmp = scipy.io.loadmat(im_f)
    
    vars=[k for k in tmp.keys() if not k[0]=='_']
    if not vars:
        print("No data found in file",os.path.join(results_dir,f))
    elif len(vars)>1:
        print("multiple variables found in file, using: ",vars[0])
    varname=vars[0]
    
    im_flat = tmp[varname]
    
    o_shape = im_flat.shape
    print(o_shape)
#    n_shape = (int(o_shape[0]),int(o_shape[1]/3),3)
    n_shape = (int(o_shape[0]/3),int(o_shape[1]),3)
    print(n_shape)
    im_len = n_shape[0]*n_shape[1]
    
    
    
#    im = np.reshape( im_flat, n_shape)
#    im = np.stack([im_flat[:,0:n_shape[1]], im_flat[:,n_shape[1]:2*n_shape[1]], im_flat[:,2*n_shape[1]:]],2)
    im = np.stack([im_flat[0:n_shape[0],:], im_flat[n_shape[0]:2*n_shape[0],:], im_flat[2*n_shape[0]:,:]],2)
        
    print(im.shape)
    
    imsave(os.path.join(image_dir,os.path.basename(im_f)[:-4]+".png"),im)


