import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats
import numpy as np
import math
import operator
from matplotlib import rc
import matplotlib.lines as mlines
import os.path

rc('font', **{'family':'serif','serif':['sffamily']})
rc('text', usetex=True)
params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

def set_xaxis(axis,limits,label,units):
  '''
  Set limits, scientific notation and offsets.
  '''
  axis.set_xlim(limits)
  offset = math.floor(np.log10(max(abs(limits[0]),abs(limits[1]))))
  axis.ticklabel_format(axis='x', style='sci',scilimits=(offset,offset))
  if offset==0:
    axis.set_xlabel(label+' ('+units+')')
  else:
    axis.set_xlabel(label+' ('+r"$\times10^{{{}}}$".format(offset)+' '+units+')')
  axis.xaxis.offsetText.set_visible(False)

def set_yaxis(axis,limits,label,units):
  '''
  Set limits, scientific notation and offsets.
  '''
  axis.set_ylim(limits)
  offset = math.floor(np.log10(max(abs(limits[0]),abs(limits[1]))))
  axis.ticklabel_format(axis='y', style='sci',scilimits=(offset,offset))
  if offset==0:
    axis.set_ylabel(label+' ('+units+')')
  else:
    axis.set_ylabel(label+' ('+r"$\times10^{{{}}}$".format(offset)+' '+units+')')
  axis.yaxis.offsetText.set_visible(False)

def main():
  # Define some static data.
  colours = {
    'turquoise':[102/255,194/255,165/255],
    'orange'   :[252/255,141/255, 98/255],
    'blue'     :[141/255,160/255,203/255],
    'purple'   :[231/255,138/255,195/255],
    'green'    :[166/255,216/255, 84/255],
    'yellow'   :[255/255,217/255, 47/255],
    'beige'    :[229/255,196/255,148/255],
    'grey'     :[179/255,179/255,179/255]}
  
  # Read in data.
  names  = ['Interpolated effective harmonic', 'VSCF', 'Harmonic', 'Uninterpolated effective harmonic', 'vscha', 'interpolated vscf']
  fnames = ['interpolated_vscha_', 'vscf_', '', 'vscha_', 'vscha_vscf_', 'interpolated_vscf_']
  dashes = [[8,2], [1,6,1,6], [10,0], [1,1,1,1,1,1,1,1,1,1,1,1,1,1], [2,5,2,5], [10,0]]
  widths = [2,3,2,2,2,1]
  
  calculating_stress = False
  
  data = []
  for name,fname,dash,width in zip(names,fnames,dashes,widths):
    file_name = fname+'thermodynamic_variables.dat'
    if os.path.isfile(file_name):
      data.append({'name':name})
      variables_file = [line.rstrip('\n').split() for line in open(file_name)]
      if 'stress' in variables_file[0]:
        calculating_stress = True
      data[-1]['dashes'] = dash
      data[-1]['width']  = width
      data[-1]['thermal energies'] = []
      data[-1]['energies'] = []
      data[-1]['free energies'] = []
      data[-1]['entropies'] = []
      if calculating_stress:
        data[-1]['stresses'] = []
        data[-1]['pressures'] = []
        data[-1]['volumes'] = []
        data[-1]['enthalpies'] = []
        data[-1]['gibbs free energies'] = []
      for line in variables_file[1:]:
        data[-1]['thermal energies'].append(float(line[0]))
        data[-1]['energies'].append(float(line[1]))
        data[-1]['free energies'].append(float(line[2]))
        data[-1]['entropies'].append(float(line[3]))
        if calculating_stress:
          stress = np.reshape([float(x) for x in line[4:13]],[3,3])
          data[-1]['stresses'].append(stress)
          data[-1]['pressures'].append(np.trace(stress)/3)
          data[-1]['volumes'].append(float(line[13]))
          data[-1]['enthalpies'].append(float(line[14]))
          data[-1]['gibbs free energies'].append(float(line[15]))
  # Plot everything.
  fig, ax_grid = plt.subplots(4,
                              1,
                              sharex=True,
                              gridspec_kw={'height_ratios':[3,3,3,1]})
  
  axes = {'energy' : {'l':ax_grid[0],
                      'b':ax_grid[0],
                      'r':ax_grid[0].twinx(),
                      't':ax_grid[0].twiny()},
          'entropy': {'l':ax_grid[1],
                      'b':ax_grid[1],
                      'r':ax_grid[1].twinx(),
                      't':ax_grid[1].twiny()},
          'stress' : {'l':ax_grid[2],
                      'b':ax_grid[2],
                      'r':ax_grid[2].twinx(),
                      't':ax_grid[2].twiny()},
          'legend' :  ax_grid[3]}
  
  xmin_hartree = None
  xmax_hartree = None
  ymin_hartree = None
  ymax_hartree = None
  ymin_shannon = 0
  ymax_shannon = None
  
  for datum in data:
    axes['energy']['l'].plot(datum['thermal energies'],
                             datum['energies'],
                             linewidth=datum['width'],
                             dashes=datum['dashes'],
                             color=colours['turquoise'])
    axes['energy']['r'].plot(datum['thermal energies'],
                             datum['free energies'],
                             linewidth=datum['width'],
                             dashes=datum['dashes'],
                             color=colours['orange'])
    axes['entropy']['l'].plot(datum['thermal energies'],
                              datum['entropies'],
                              linewidth=datum['width'],
                              dashes=datum['dashes'],
                              color=colours['blue'])
    if calculating_stress:
      axes['energy']['l'].plot(datum['thermal energies'],
                               datum['enthalpies'],
                               linewidth=datum['width'],
                               dashes=datum['dashes'],
                               color=colours['grey'])
      axes['energy']['r'].plot(datum['thermal energies'],
                               datum['gibbs free energies'],
                               linewidth=datum['width'],
                               dashes=datum['dashes'],
                               color=colours['purple'])
      axes['stress']['l'].plot(datum['thermal energies'],
                               datum['pressures'],
                               linewidth=datum['width'],
                               dashes=datum['dashes'],
                               color=colours['green'])
    if xmin_hartree == None:
      xmin_hartree = datum['thermal energies'][0]
      xmax_hartree = datum['thermal energies'][-1]
      ymin_hartree = min(datum['free energies'])
      ymax_hartree = max(datum['energies'])
      ymax_shannon = max(datum['entropies'])
      ymin_pressure = min(datum['pressures'])
      ymax_pressure = max(datum['pressures'])
    else:
      xmin_hartree = min(xmin_hartree, datum['thermal energies'][0])
      xmax_hartree = max(xmax_hartree, datum['thermal energies'][-1])
      ymin_hartree = min(ymin_hartree, min(datum['free energies']))
      ymax_hartree = max(ymax_hartree, max(datum['energies']))
      ymax_shannon = max(ymax_shannon, max(datum['entropies']))
      ymin_pressure = min(ymin_pressure, min(datum['pressures']))
      ymax_pressure = max(ymax_pressure, max(datum['pressures']))
    if calculating_stress:
      ymin_hartree = min(ymin_hartree, min(datum['gibbs free energies']))
      ymax_hartree = max(ymax_hartree, max(datum['enthalpies']))
  
  ydiff_hartree = ymax_hartree-ymin_hartree
  ymin_hartree = ymin_hartree-0.05*ydiff_hartree
  ymax_hartree = ymax_hartree+0.05*ydiff_hartree
  
  ymax_shannon = 1.1*ymax_shannon
  
  ydiff_pressure = ymax_pressure-ymin_pressure
  ymin_pressure = ymin_pressure-0.05*ydiff_pressure
  ymax_pressure = ymax_pressure+0.05*ydiff_pressure
  
  kb_in_ev_per_k = 8.6173303e-5
  ev_per_hartree = 27.21138602
  kb_in_au = kb_in_ev_per_k / ev_per_hartree
  
  hartree_per_cubic_bohr_to_pascals = 29421.02648438959e9
  
  xmin_kelvin  = xmin_hartree/kb_in_au
  xmax_kelvin  = xmax_hartree/kb_in_au
  axes['energy']['b'].set_xlim(xmin_hartree,xmax_hartree)
  axes['energy']['t'].set_xlim(xmin_kelvin,xmax_kelvin)
  axes['energy']['t'].set_xlabel(r"$T$ (K)")
  axes['energy']['t'].tick_params(labelbottom=False,labeltop=True)
  axes['energy']['b'].tick_params(labelbottom=False,labeltop=False)
  axes['entropy']['b'].set_xlim(xmin_hartree,xmax_hartree)
  axes['entropy']['t'].set_xlim(xmin_kelvin,xmax_kelvin)
  axes['entropy']['t'].tick_params(labelbottom=False,labeltop=False)
  axes['entropy']['b'].tick_params(labelbottom=False,labeltop=False)
  axes['stress']['b'].set_xlim(xmin_hartree,xmax_hartree)
  axes['stress']['t'].set_xlim(xmin_kelvin,xmax_kelvin)
  axes['stress']['t'].tick_params(labelbottom=False,labeltop=False)
  axes['stress']['b'].tick_params(labelbottom=True,labeltop=False)
  axes['stress']['b'].set_xlabel(r"$k_BT$ (Ha)")
  set_xaxis(axes['stress']['b'],
            (xmin_hartree,xmax_hartree),
            r"$k_BT$",
            r"Ha")
  
  set_yaxis(axes['energy']['l'],
            (ymin_hartree,ymax_hartree),
            r"",
            r"Ha/cell")
  set_yaxis(axes['energy']['r'],
            (ymin_hartree,ymax_hartree),
            r"",
            r"Ha/cell")
  
  ymin_gibbs = ymin_shannon * kb_in_au
  ymax_gibbs = ymax_shannon * kb_in_au
  set_yaxis(axes['entropy']['l'],
            (ymin_shannon,ymax_shannon),
            r"$S/k_B$",
            r"arb. units/cell")
  set_yaxis(axes['entropy']['r'],
            (ymin_gibbs,ymax_gibbs),
            r"$S$",
            r"Ha/K/cell")
  
  ymin_pascals = ymin_pressure * hartree_per_cubic_bohr_to_pascals
  ymax_pascals = ymax_pressure * hartree_per_cubic_bohr_to_pascals
  set_yaxis(axes['stress']['l'],
            (ymin_pressure,ymax_pressure),
            r"$p$",
            r"Ha/bohr$^3$")
  set_yaxis(axes['stress']['r'],
            (ymin_pascals,ymax_pascals),
            r"$p$",
            r"Pa")
  
  
  # Format legend.
  handles = []
  labels = []
  
  things = [('energy','turquoise'),
            ('free energy','orange'),
            ('enthalpy','grey'),
            ('gibbs free energy','purple'),
            ('entropy','blue'),
            ('pressure','green')]
  if not calculating_stress:
    things = things[1,2,5]
  for label,colour in things:
    line = mlines.Line2D([],[],color=colours[colour])
    handles.append(line)
    labels.append(label)
  for datum in data:
    line = mlines.Line2D([],[],color=[0,0,0],dashes=datum['dashes'])
    handles.append(line)
    labels.append(datum['name'])
  
  axes['legend'].axis('off')
  axes['legend'].legend(handles=handles, labels=labels, loc='center', ncol=3)
  
  # Show plot.
  plt.show()

main()
