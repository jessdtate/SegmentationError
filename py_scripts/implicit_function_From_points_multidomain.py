import os
import threading
#import numpy as np


def implicit_model(vol_fname, points, outputpath, out_name):

    # for running multiple files without restarting.
    # this doesn't work with command line scripting calls.
    
#    active_layer = get(stateid = "layermanager::active_layer")
#    print(active_layer)
#    while not active_layer=='[]':
#        print("deleting: ", active_layer)
#        deletelayers(layers = [active_layer])
#        active_layer = get(stateid = "layermanager::active_layer")
        
    
    layers = importlayer(filename=vol_fname, importer='[Teem Importer]', mode='data')
    wait_for_layer(layers[0])

    # assume the regions are split: LV 1/4 , RV 1/4, epi 1/2
    num_pnts = len(points)
    print(num_pnts)
    points_lv = points[:int(num_pnts/4)]
    points_rv = points[int(num_pnts/4):int(num_pnts/2)]
    points_epi = points[int(num_pnts/2):]
    
    pnt_str_lv = points2str(points_lv)
    pnt_str_rv = points2str(points_rv)
    pnt_str_epi = points2str(points_epi)
    

    
  
#    implicitmodel(layerid='layer_0',vertices='[[-17.4015,15,41.4717],[-22.6839,15,-13.1235],[35.4226,15,-32.496],[61.8347,15,0.378543],[47.1613,15,24.4474]]',view_modes='[[coronal],[coronal],[coronal],[coronal],[coronal]]',normal_offset='2',convex_hull_2D='false',invert_seed_order='false',kernel='thin_plate')
    
#    "implicitmodel(layerid='layer_0',vertices=[[-31.0487,-4.20975,-32.0863],[16.1315,-9.55522,-40.3893],[-46.5695,-23.3927,-16.5319],[-30.1696,11.2075,-30.2825],[3.62445,-22.3768,-38.5409],[30.0292,-29.6416,-40.7586],[-34.2162,-26.5468,20.8768],[-15.4794,31.8742,33.2607],[-38.0409,-28.3396,5.51522],[53.7664,2.50005,-35.6904]], view_modes='[[coronal]]',normal_offset='2',convex_hull_2D='false',invert_seed_order='false',kernel='thin_plate')"
  
    command_lv = "im_model = implicitmodel(layerid=layers[0], vertices=" + pnt_str_lv + ",  view_modes='[[coronal]]', normal_offset= 2)"
    command_rv = "im_model = implicitmodel(layerid=layers[0], vertices=" + pnt_str_rv + ",  view_modes='[[coronal]]', normal_offset= 2)"
    command_epi = "im_model = implicitmodel(layerid=layers[0], vertices=" + pnt_str_epi + ",  view_modes='[[coronal]]', normal_offset= 2)"
#    exec(command)

    print(outputpath)
    print(out_name)
    
    [im_model_lv, seg_im_lv, seg_vol_lv] = run_implicit_model(layers[0], points_lv, outputpath, out_name+'_lv')
    [im_model_rv, seg_im_rv, seg_vol_rv] = run_implicit_model(layers[0], points_rv, outputpath, out_name+'_rv')
    [im_model_epi, seg_im_epi, seg_vol_epi] = run_implicit_model(layers[0], points_epi, outputpath, out_name+'_epi')
    
    #im_model_air = invert(layerid = im_model_epi, replace='false')
    # invert layer not working as expected, use arithmatic filter
    im_model_air = arithmeticfilter(layerids=[im_model_epi], expressions="[RESULT = -A;]", output_type='data', replace='false',preserve_data_format='false')
    wait_for_layer(im_model_air)
    set(stateid=im_model_air+'::name',value=out_name+'_air_vol')
    result=exportlayer(layer=im_model_air, file_path = os.path.join(outputpath,out_name+'_air_vol'), exporter='[NRRD Exporter]', extension='.nrrd')
    
    if not check_for_crossing(im_model_air):
        raise ValueError("Warning:  air layer not valid")
    
    [im_model, seg_im, seg_vol] = combine_domains(im_model_lv, im_model_rv, im_model_epi, outputpath, out_name)
    
    
    
    
    
    
#    deletelayers(layers = [seg_im_lv,seg_im_rv,seg_im_epi,seg_im])
    
    return im_model_lv, im_model_rv, im_model_epi, im_model
def check_for_crossing(im_model):
    im_max = get(stateid = im_model+'::max')
    im_min = get(stateid = im_model+'::min')
    
    if im_min<0 and im_max>0:
        return True
    else:
        return False
    
def combine_domains(im_model_lv, im_model_rv, im_model_epi, outputpath, out_name):

    union_text = "[A = -A;\nB = -B;\nRESULT = - select(A<0,select(B<0,select(A<B,A,B),A),select(B<0,B,select(B==0,B,select(A==0,A,select(A<B,A,B)))));]"
    remove_text = "[A = -A;\nB = -B;\nRESULT = - select(A<0,select(B<=0,-B,select(A>-B,A,-B)),select(A==0,select(B<0,-B,A),select(B<0,select(A<-B,A,-B),A)));]"
    
    im_model_endo = arithmeticfilter(layerids=[im_model_lv, im_model_rv], expressions=union_text, output_type='data', replace='false',preserve_data_format='false')
    
    wait_for_layer(im_model_endo)
    
    set(stateid=im_model_endo+'::name',value=out_name+'_endo_vol')
    
    im_model_rough = arithmeticfilter(layerids=[im_model_epi, im_model_endo], expressions=remove_text,output_type='data',replace='false',preserve_data_format='false')
    
    wait_for_layer(im_model_rough)
    
    set(stateid=im_model_rough+'::name',value=out_name+'_vent_vol')
    result=exportlayer(layer=im_model_rough, file_path = os.path.join(outputpath,out_name+'_vent_vol'), exporter='[NRRD Exporter]', extension='.nrrd')
    
    print(outputpath)
    print(out_name)
    
    upper_threshold = get(stateid = im_model_rough+'::max')
    lower_threshold = 0
    
    seg_im_rough = threshold(layerid=im_model_rough, lower_threshold='0', upper_threshold=upper_threshold)
    wait_for_layer(seg_im_rough)
    seg_vol = calculatemaskvolume(mask = seg_im_rough)

    if get(stateid = seg_im_rough+'::calculated_volume') == 0:
      print("WARNING: No segmentation")
      
    computeisosurface(layerid=seg_im_rough)
    exportisosurface(layer=seg_im_rough, file_path=outputpath+out_name+'_vent.fac', binary='false')
    
#    im_model = medianfilter(layerid = im_model_rough, radius=1)
#
#    wait_for_layer(im_model)
#
#    set(stateid=im_model+'::name',value=out_name+'_vent_smth_vol')
#    result=exportlayer(layer=im_model, file_path = os.path.join(outputpath,out_name+'_vent_smth_vol'), exporter='[NRRD Exporter]', extension='.nrrd')
#
#    upper_threshold = get(stateid = im_model+'::max')
#    lower_threshold = 0
#
#    seg_im = threshold(layerid=im_model, lower_threshold='0', upper_threshold=upper_threshold)
#    wait_for_layer(seg_im)
#    seg_vol = calculatemaskvolume(mask = seg_im)
#
#    if get(stateid = seg_im+'::calculated_volume') == 0:
#      print("WARNING: No segmentation")
#
#    computeisosurface(layerid=seg_im)
#    exportisosurface(layer=seg_im, file_path=outputpath+out_name+'_vent_smth.fac', binary='false')
    
    return im_model_rough, seg_im_rough, seg_vol
    
    
  
def load_points(filename):

    fid = open(filename, "r")
    str_dump = fid.read()
    fid.close()

    lines = str_dump.split('\n')

#    pnts = [ [float(p) for p in l.split('   ')] for l in lines[:-1]]
    pnts = [ [float(p) for p in l.split(' ')] for l in lines[:-1]]

#    pnts=np.loadtxt(filename)

#    return pnts.tolist()
    return pnts

def points2str(points):
    pnt_str_l = [ '['+str(p[0]) + ','+str(p[1]) + ',' + str(p[2]) + ']' for p in points ]
    pnt_str = '[' + ','.join(pnt_str_l)+']'
    return pnt_str
    
def run_implicit_model(layer_id, points, outputpath, out_name):

    im_model = implicitmodel(layerid=layer_id,vertices=points, view_modes='[[coronal]]', normal_offset='2',convex_hull_2D='false',invert_seed_order='false',kernel='thin_plate')
  
    print(im_model)
    wait_for_layer(im_model)
    set(stateid=im_model+'::name',value=out_name+'_vol')
    result=exportlayer(layer=im_model, file_path = os.path.join(outputpath,out_name+'_vol'), exporter='[NRRD Exporter]', extension='.nrrd')
    
    upper_threshold = get(stateid = im_model+'::max')
    lower_threshold = 0

    seg_im = threshold(layerid=im_model, lower_threshold='0', upper_threshold=upper_threshold)
    wait_for_layer(seg_im)
    seg_vol = calculatemaskvolume(mask = seg_im)

    if get(stateid = seg_im+'::calculated_volume') == 0:
      print("WARNING: No segmentation")
      
    computeisosurface(layerid=seg_im)
    exportisosurface(layer=seg_im, file_path=outputpath+out_name+'.fac', binary='false')
    
    return im_model, seg_im, seg_vol

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
