import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats
import numpy as np
import operator
from matplotlib import rc
import matplotlib.lines as mlines
import os.path

rc('font', **{'family':'serif','serif':['sffamily']})
rc('text', usetex=True)
params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

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
  names  = ['Interpolated effective harmonic', 'VSCF', 'Harmonic', 'Uninterpolated effective harmonic', 'vscha']
  fnames = ['interpolated_vscha_', 'vscf_', '', 'vscha_', 'vscha_vscf_']
  dashes = [[8,2], [10,0], [10,0], [1,1,1,1,1,1,1,1,1,1,1,1,1,1], [2,5,2,5]]
  widths = [2,1,1,2,2]
  
  data = []
  for name,fname,dash,width in zip(names,fnames,dashes,widths):
    file_name = fname+'thermodynamic_variables.dat'
    if os.path.isfile(file_name):
      data.append({'name':name})
      variables_file = [line.rstrip('\n').split() for line in open(file_name)]
      data[-1]['dashes'] = dash
      data[-1]['width']  = width
      data[-1]['thermal energies'] = []
      data[-1]['energies'] = []
      data[-1]['free energies'] = []
      data[-1]['entropies'] = []
      for line in variables_file[1:]:
        data[-1]['thermal energies'].append(float(line[0]))
        data[-1]['energies'].append(float(line[1]))
        data[-1]['free energies'].append(float(line[2]))
        data[-1]['entropies'].append(float(line[3]))
  # Plot everything.
  fig, ax_grid = plt.subplots(3,
                              1,
                              sharex=True,
                              gridspec_kw={'height_ratios':[3,3,1]})
  
  axes = {'energy' : {'l':ax_grid[0],
                      'b':ax_grid[0],
                      'r':ax_grid[0].twinx(),
                      't':ax_grid[0].twiny()},
          'entropy': {'l':ax_grid[1],
                      'b':ax_grid[1],
                      'r':ax_grid[1].twinx(),
                      't':ax_grid[1].twiny()},
          'legend':  ax_grid[2]}
  
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
                              color=colours['purple'])
    if xmin_hartree == None:
      xmin_hartree = datum['thermal energies'][0]
      xmax_hartree = datum['thermal energies'][-1]
      ymin_hartree = min(datum['free energies'])
      ymax_hartree = max(datum['energies'])
      ymax_shannon = max(datum['entropies'])
    else:
      xmin_hartree = min(xmin_hartree, datum['thermal energies'][0])
      xmax_hartree = max(xmax_hartree, datum['thermal energies'][-1])
      ymin_hartree = min(ymin_hartree, min(datum['free energies']))
      ymax_hartree = max(ymax_hartree, max(datum['energies']))
      ymax_shannon = max(ymax_shannon, max(datum['entropies']))
  
  ydiff_hartree = ymax_hartree-ymin_hartree
  ymin_hartree = ymin_hartree-0.05*ydiff_hartree
  ymax_hartree = ymax_hartree+0.05*ydiff_hartree
  
  ymax_shannon = 1.1*ymax_shannon
  
  kb_in_ev_per_k = 8.6173303e-5
  ev_per_hartree = 27.21138602
  kb_in_au = kb_in_ev_per_k / ev_per_hartree
  
  xmin_kelvin  = xmin_hartree/kb_in_au
  xmax_kelvin  = xmax_hartree/kb_in_au
  axes['energy']['b'].set_xlim(xmin_hartree,xmax_hartree)
  axes['energy']['t'].set_xlim(xmin_kelvin,xmax_kelvin)
  axes['energy']['l'].set_ylim(ymin_hartree,ymax_hartree)
  axes['energy']['r'].set_ylim(ymin_hartree,ymax_hartree)
  axes['energy']['t'].set_xlabel(r"Temperature, $T$, (K)")
  axes['energy']['b'].tick_params(labelbottom='off')
  axes['entropy']['b'].set_xlim(xmin_hartree,xmax_hartree)
  axes['entropy']['t'].set_xlim(xmin_kelvin,xmax_kelvin)
  axes['entropy']['t'].tick_params(labeltop='off')
  axes['entropy']['b'].set_xlabel(r"$k_BT$ (Ha)")
  
  axes['energy']['l'].set_ylabel(r"Vibrational Energy per cell (Ha)")
  axes['energy']['l'].tick_params(axis='y',
                                  width=2,
                                  color=colours['turquoise'])
  axes['energy']['l'].spines['left'].set_color(colours['turquoise'])
  axes['energy']['r'].set_ylabel(r"Vibrational Free Energy per cell (Ha)")
  axes['energy']['r'].tick_params(axis='y',
                                  width=2,
                                  color=colours['orange'])
  for _,ax in axes['energy'].items():
    ax.spines['left'].set_color(colours['turquoise'])
    ax.spines['right'].set_color(colours['orange'])
  
  ymin_gibbs   = ymin_shannon * kb_in_au
  ymax_gibbs   = ymax_shannon * kb_in_au
  axes['entropy']['l'].set_ylim(ymin_shannon,ymax_shannon)
  axes['entropy']['l'].set_ylabel(r"Vibrational Shannon Entropy per cell, $S/k_B$, (arb. units)")
  axes['entropy']['r'].set_ylim(ymin_gibbs,ymax_gibbs)
  axes['entropy']['r'].set_ylabel(r"Vibrational Gibbs Entropy per cell, $S$, (Ha per K)")
  
  # Format legend.
  handles = []
  labels = []
  
  things = [('energy','turquoise'),
            ('free energy','orange'),
            ('entropy','purple')]
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
