import scipy.io
import os
import fnmatch
import numpy as np
from skimage.io import imsave

image_dir = "/Users/jess/software/ShapeWorks/Examples/Python/Test_epi_centered/image"

image_files = fnmatch.filter(os.listdir(image_dir), '*.txt')

for im_f in image_files:
    print(im_f)
    
    im_flat = np.loadtxt(os.path.join(image_dir,im_f))
    
    o_shape = im_flat.shape
    n_shape = (int(o_shape[0]/3),o_shape[1],3)
    
    im_len = n_shape[0]*n_shape[1]
    
    
    
#    im = np.reshape( im_flat, n_shape)
    im = np.stack([im_flat[0:n_shape[0],:], im_flat[n_shape[0]:2*n_shape[0],:], im_flat[2*n_shape[0]:,:]],2)
    
    print(im.shape)
    
    imsave(os.path.join(image_dir,im_f+".png"),im)


