import os
import threading


def upsample_heart(Filename, outputfilename,name):

  layers = importlayer(filename=Filename, importer='[Teem Importer]', mode='label_mask')
  wait_for_layer(layers[0])
  
  dims = get(stateid = 'group_0::dimensions')

  r_im1=resample(layerids=[layers[0]],x=dims[0],y=dims[1],z=dims[2]*4,crop='false',kernel='box',replace='false')
  wait_for_layer(r_im1[0])
  #de_im = erodefilter(layerid = r_im1[0], radius = 1, replace = True)
  #wait_for_layer(de_im)
  de_im = dilatefilter(layerid = r_im1[0], radius = 4, replace = True)
  wait_for_layer(de_im)
  de_im = erodefilter(layerid = de_im, radius = 5, replace = True)
  wait_for_layer(de_im)
  de_im = dilatefilter(layerid = de_im, radius = 1, replace = True)
  wait_for_layer(de_im)

#r_im2=resample(layerids=[de_im],x=dims[0],y=dims[1],z=dims[2]*4,crop='false',kernel='box',replace='false')

#wait_for_layer(r_im2[0])
  
  #de2_im = erodefilter(layerid = r_im2[0], radius = 2, replace = True)
  #wait_for_layer(de2_im)
  #de2_im = dilatefilter(layerid = r_im2[0], radius = 2, replace = True)
  #wait_for_layer(de2_im)
  #de2_im = erodefilter(layerid = de2_im, radius = 4, replace = True)
  #wait_for_layer(de2_im)
  #de2_im = dilatefilter(layerid = de2_im, radius = 2, replace = True)
  #wait_for_layer(de2_im)

  c_im = crop(layerids=[de_im],origin='[-64.8343,-299.594,-314]',size='[189.818,189.947,159]',replace='false')
  wait_for_layer(c_im[0])
  set(stateid=c_im[0]+'::name',value=name)
  result=exportsegmentation(layers=c_im[0],file_path=outputfilename ,exporter='[NRRD Exporter]',extension='.nrrd')

def torso_points(Filename, outputfilename,name):

  layers = importlayer(filename=Filename, importer='[Teem Importer]', mode='label_mask')
  wait_for_layer(layers[0])

  computeisosurface(layerid=layers[0],quality_factor=0.25)
  set(stateid=layers[0]+'::name',value=name)
  exportisosurface(layer=layers[0],file_path=outputfilename+name+'.fac',binary='false')


def wait_for_layer(layer,MAX_ITERATIONS = 1000):
  
  #checks to make sure that the layer is available before trying operations on it.
  TIMEOUT=1
  c = threading.Condition()
  c.acquire()
  
  layerStatus = get(stateid=layer+"::data")
  counter = 1
  success = 1
  with c:
    while not layerStatus == "available":
      if counter > MAX_ITERATIONS:
        print('waited too long for '+layer)
        success = 0
        break
      counter += 1
      c.wait(TIMEOUT)
      layerStatus = get(stateid=layer+"::data")
  #print('waiting for '+layer)

  return success
