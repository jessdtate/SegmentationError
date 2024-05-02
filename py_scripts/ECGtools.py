import numpy as np
import math
import scipy.signal as sig
from scipy.fft import fft, fftfreq, ifft, fftshift

import argparse

import plotly.express as px
import plotly.offline as pyo
from plotly.subplots import make_subplots
import plotly.graph_objects as go
  
  
  
def dominantFrequency(pots, sampling_rate=1):

  # ref: https://www.ahajournals.org/doi/full/10.1161/01.RES.86.4.408
  
  DF = []
  FS = []
  
#    print(pots.shape)
  N = pots.shape[1]
  xf = fftfreq(N, 1/sampling_rate)
  
  for p in pots:
#        print(p.shape)
    fspec = fft(p)
#        print(np.argmax(np.abs(fspec[1:]))+1)
#        print(fspec[np.argmax(np.abs(fspec[1:]))+1])
    df = xf[np.argmax(np.abs(fspec[1:]))+1]
    
#        print(df)
    
    DF.append(df)
    FS.append(fspec)
    
  return np.array(DF), np.array(FS), xf
  
def DVDT(signals, **kwargs):
  defaultKwargs = { "sampling_rate" : 1,
                    "window" : 20,
                    "deg" : 3
  }
  kwargs = { **defaultKwargs, **kwargs }
  
  deg = kwargs["deg"]
  window = kwargs["window"]
  

  T = signals.shape
  if len(T)>1:
    raise ValueError("signal should be 1D")

  cent = math.ceil(window/2)
  X = np.zeros((window,deg+1))
  L = np.arange(-cent,cent)
  for p in range(1,deg+2):
    X[:,p-1] = L**((deg+1)-p)
      
  E = np.dot(np.linalg.inv(np.dot(X.T,X)),X.T)
  
  temp = np.hstack((signals,signals[-1]*np.ones((cent))))
  
  a = sig.lfilter(E[deg-1,np.arange(window-1,-1,-1)],1,temp)
  
  dy = a[cent-1:-1]/kwargs["sampling_rate"]
  
  return dy
  
def findLAT(signals, **kwargs):

  defaultKwargs = { "sampling_rate" : 1,
                    "window" : 20,
                    "deg" : 3,
                    "ignore_first" : 100,
                    "n_LAT" : 1,
                    "mode" : "min",
                    "fids" : {}
  }
  kwargs = { **defaultKwargs, **kwargs }
  
  
  print(kwargs)

  [M,T] = signals.shape
  
  print(M,T)
  
  tau = []
  dy = np.zeros(signals.shape)
  
  ig_factor = 0.5
  
  sampling_rate = kwargs["sampling_rate"]
  ignore_first = kwargs["ignore_first"]
  window = kwargs["window"]
  deg = kwargs["deg"]
  mx_n_peaks = 0
  
  if kwargs["fids"]:
    fids = kwargs["fids"]
    qoff = int(fids["qoff"]*fids["sf"][0,0]/1000 - ignore_first*ig_factor)
  else:
    qoff = T-int(ignore_first*ig_factor)
      
  for k in range(M):
    dy[k,:] = DVDT(signals[k,:], **kwargs)
    peaks = findNPeaksDVDT(dy[k,ignore_first:qoff], **kwargs)
    n_peaks = len(peaks)
    
    if n_peaks > mx_n_peaks:
      mx_n_peaks = n_peaks
      
    tau.append(peaks+ignore_first)
    
#    print(k)

  print(mx_n_peaks)
      
  if kwargs["n_LAT"]>0:
    tau = np.array(tau)
      
  return tau, dy, mx_n_peaks
  
def findNPeaksDVDT(dvdt, **kwargs):

  defaultKwargs = { "window" : 20,
                    "n_LAT" : 1,
                    "mode" : "min"
  }
  kwargs = { **defaultKwargs, **kwargs }

  
  if kwargs["mode"] == "min":
    check = 0
  elif kwargs["mode"] == "max":
    check = -1
  else:
    raise ValueError("mode not recognized")
  
  n_LAT = kwargs["n_LAT"]
  if n_LAT<1:
    n_LAT = len(dvdt)
    
  std_thresh = 0.5

  dv = dvdt
  k = 0
  peaks = []
  
  sorted = np.argsort(dv)
  while k<n_LAT:
    
    peak = sorted[check]
    
    if kwargs["mode"] == "max" and (dvdt[peak] < (np.mean(dvdt) + std_thresh*np.std(dvdt))):
      break
      
    if kwargs["mode"] == "min" and (dvdt[peak] > (np.mean(dvdt) - std_thresh*np.std(dvdt))):
      break
    
    
    peaks.append(peak)
    
    wb = np.max([ peak-math.floor(kwargs["window"]/2), 0])
    we = np.min([ peak+math.ceil(kwargs["window"]/2), len(dvdt)])
    
    sorted = sorted[(sorted>we) |  (sorted<wb)]
    
#    print(sorted)
#    print(dvdt[sorted])
    
    if len(sorted) == 0:
      break
    
#    dv[wb:we] = 0
    k+=1
      
  return np.array(peaks)
  
  
def phaseAnalysis(pots, sampling_rate = 1):

  #ref:
  #https://www.ahajournals.org/doi/full/10.1161/circep.110.853804
  
  HT =[]
  phases = []
  cycles = []
  
  mx_cycles = 0
  
  for p in pots:
    ht = sig.hilbert(p)
    HT.append(ht)
    
    phase = np.arctan2(ht.imag, ht.real)
    phases.append(phase)
    
    cycle = zeroCrossings(phase)
    n_cycles = len(cycle[0][0])
    cycles.append(cycle[0])

    mean_phase = np.mean(phase)
    stdev_phase = np.std(phase)
    
    if n_cycles > mx_cycles:
      mx_cycles = n_cycles
    
#        print(mean_phase, stdev_phase)
      
  
  return np.array(HT), np.array(phases), cycles, mx_cycles
  
  
  
def zeroCrossings(signal):

  d_kern = np.array([1, 0, -1])
  d_s = np.convolve(signal, d_kern, mode = "same")
  
  sign_ch = np.where((signal[:-1]*signal[1:])<= 0)
  
  
  neg_cross = np.where(((signal[:-1]*signal[1:])<= 0) * d_s[:-1]<0)
  pos_cross = np.where(((signal[:-1]*signal[1:])<= 0) * d_s[:-1]>0)
  
  return neg_cross, pos_cross
  
def cycle_lengths(cycles):
# one spatiotemporal data

  mean_cylen = []
  for cyc in cycles:
    print(len(cyc))
    print(cyc.shape)
    print(cyc)
    
    cy = np.sort(cyc[cyc>0])
    
    print(cy)
    
    cy_itv = cy[1:]-cy[:-1]
    
    mean_cylen.append(np.mean(cy_itv))
    
  return np.array(mean_cylen)

def plot_PM_check(pots, **kwargs):
  defaultKwargs = { "HT" : np.array([]),
                    "phase": np.array([]),
                    "cycles" : np.array([]),
                    "mx_cycles" : 0,
                    "channels" : [0],
                    "sampling_rate" : 1
  }
  kwargs = { **defaultKwargs, **kwargs }
  
  HT = kwargs["HT"]
  phase = kwargs["phase"]
  cycles = kwargs["cycles"]
  mx_cycles = kwargs["mx_cycles"]
  sampling_rate = kwargs["sampling_rate"]
  channels = kwargs["channels"]
  
  num_col = 1
  num_row = 4

  if len(HT)==0 or len(phase)==0 or len(cycles)==0:
    HT, phase, cycles, mx_cycles = phaseAnalysis(pots, sampling_rate)

  if mx_cycles <= 0:
    mx_cycles = 0
    for cy in cycles:
      n_cy = len(cy[0])
      mx_cycles = max([mx_cycles, n_cy])
  
  
#  check = checkPhase(pots, HT, phase, cycles, mx_cycles)
  
  data_shape=np.shape(pots)
  t = list(range(data_shape[1]))
  
  figs = []
  for ch in channels:
    t_cyc = [t[int(cy)] for cy in cycles[ch][0]]
        
    fig = make_subplots(rows=num_row, cols=num_col)
    fig.update_layout(height=600, width=600,
                  title_text="Channel "+str(ch))
                  
    fig.add_trace(go.Scatter(x=t, y=pots[ch,:], name = "Pots"),
                  row=1, col=1 )
    
    fig.add_trace(go.Scatter(x=t, y=HT[ch,:].real, name = "v(t+T)"  ),
                  row=2, col=1 )
    fig.add_trace(go.Scatter(x=t, y=HT[ch,:].imag, name = "v(t)" ),
                  row=2, col=1 )
    fig.add_trace(go.Scatter(x=t, y=np.abs(HT[ch,:]), name = "envelope" ),
                  row=2, col=1 )
    
    fig.add_trace(go.Scatter(x=t, y=phase[ch,:], name = "phase"),
                  row=3, col=1 )
    fig.add_trace(go.Scatter(x=t_cyc, y=[0]*len(t_cyc), mode= "markers", name = "cycles"),
                  row=3, col=1 )

    fig.add_trace(go.Scatter(x=HT[ch,:].real, y=HT[ch,:].imag, name = "phase space"),
                  row=4, col=1 )
    
    fig.update_xaxes(title_text="t (ms)", row=1, col=1)
    fig.update_yaxes(title_text="pots (mV)", row=1, col=1)
    fig.update_xaxes(title_text="t (ms)", row=2, col=1)
    fig.update_yaxes(title_text="v(t) (mV)", row=2, col=1)
    fig.update_xaxes(title_text="t (ms)", row=3, col=1)
    fig.update_yaxes(title_text="phase (rads)", row=3, col=1)
    fig.update_xaxes(title_text="v(t+T) (mV)", row=4, col=1)
    fig.update_yaxes(title_text="v(t) (mV)", row=4, col=1)
    
    fig.show()
    
    figs.append(fig)
        
  return figs
        
    


#def plot_DF_check(pots, **kwargs):
#  defaultKwargs = { "DF" : np.array([]),
#                    "spec": np.array([]),
#                    "cycles" : np.array([]),
#                    "mx_cycles" : 0,
#                    "channels" : [0],
#                    "sampling_rate" : 1
#  }
#  kwargs = { **defaultKwargs, **kwargs }
#
#  HT = kwargs["HT"]
#  phase = kwargs["phase"]
#  cycles = kwargs["cycles"]
#  mx_cycles = kwargs["mx_cycles"]
#  sampling_rate = kwargs["sampling_rate"]
#
#  num_col = 2
#  num_row = 3
#
#  if not (HT and phase and cycles):
#    HT, phase, cycles, mx_cycles = phaseAnalysis(pots, sampling_rate)
#
#  if mx_cycles <= 0:
#    mx_cycles = 0
#    for cy in cycles:
#      n_cy = len(cy[0])
#      mx_cycles = max([mx_cycles, n_cy])
#
#  data_shape=np.shape(data)
#  t = list(range(data_shape[1]))/sampling_rate
#  for ch in channels:
#
#    t_cyc = [t[int(c)] for c in cycle[ch][0]]
#
#        fig = make_subplots(rows=num_row, cols=num_col)
#        fig.update_layout(height=600, width=600,
#                      title_text="Channel "+str(ch))
#
#        fig.add_trace(go.Scatter(x=t, y=data[ch,:], name = "Pots"),
#                      row=1, col=1 )
#
#        fig.add_trace(go.Scatter(x=LATs[ch],
#                                 y=data[ch,LATs[ch].astype(int)],
#                                 mode= "markers",
#                                 name = "LATs"),
#                      row=1, col=1 )
#
#
#        fig.add_trace(go.Scatter(x=xf,
#                                 y=np.abs(fftshift(spectrum[ch,:])),
#                                 name = "fft"),
#                      row=1, col=2 )
#
#        fig.add_trace(go.Scatter(x=[DF[ch], DF[ch]],
#                                 y=[smn[ch], smx[ch]],
#                                 text="DF = {:.2f}".format(DF[ch]),
#                                 textposition="top center",
#                                 name = "DF = {:.2f}".format(DF[ch]) ),
#                      row=1, col=2 )
#
#        fig.add_trace(go.Scatter(x=t,
#                                 y=HT[ch,:].real,
#                                 name = "v(t+T)"  ),
#                      row=2, col=1 )
#        fig.add_trace(go.Scatter(x=t,
#                                 y=HT[ch,:].imag,
#                                 name = "v(t)" ),
#                      row=2, col=1 )
#        fig.add_trace(go.Scatter(x=t,
#                                 y=np.abs(HT[ch,:]),
#                                 name = "envelope" ),
#                      row=2, col=1 )
#
#        fig.add_trace(go.Scatter(x=t,
#                                 y=phase[ch,:],
#                                 name = "phase"),
#                      row=2, col=2 )
#        fig.add_trace(go.Scatter(x=t_cyc,
#                                 y=[0]*len(t_cyc),
#                                 mode= "markers",
#                                 name = "cycles"),
#                      row=2, col=2 )
#
#
#        fig.add_trace(go.Scatter(x=t,
#                                 y=dys[ch,:]*1000,
#                                 name = "dv/dt"),
#                      row=3, col=1 )
#
#        fig.add_trace(go.Scatter(x=LATs[ch],
#                                 y=dys[ch,LATs[ch].astype(int)]*1000,
#                                 mode= "markers",
#                                 name = "max"),
#                      row=3, col=1 )
#
#
#        fig.add_trace(go.Scatter(x=HT[ch,t_vt:].real,
#                                 y=HT[ch,t_vt:].imag,
#                                 name = "phase space"),
#                      row=3, col=2 )
#
#
#
#        fig.update_xaxes(title_text="t (ms)", row=1, col=1)
#        fig.update_yaxes(title_text="TMP (mV)", row=1, col=1)
#        fig.update_xaxes(title_text="freq (Hz)", range=[-20, 20], row=1, col=2)
#        fig.update_yaxes(title_text="power", row=1, col=2)
#        fig.update_xaxes(title_text="t (ms)", row=2, col=1)
#        fig.update_yaxes(title_text="phase (rads)", row=2, col=1)
#        fig.update_xaxes(title_text="v(t+T) (mV)", row=2, col=2)
#        fig.update_yaxes(title_text="v(t) (mV)", row=2, col=2)
