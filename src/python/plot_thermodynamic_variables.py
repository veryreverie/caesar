import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats
import numpy as np
import operator
from matplotlib import rc
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
  names = ['vscha_', 'vscf_', 'uninterpolated_', '']
  dashes = [[5,5], [10,2], [6,2,2,2], [1,0]]
  
  data = []
  for name,dash in zip(names,dashes):
    file_name = name+'thermodynamic_variables.dat'
    if os.path.isfile(file_name):
      data.append({'name':name})
      variables_file = [line.rstrip('\n').split() for line in open(file_name)]
      data[-1]['dashes'] = dash
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
  fig, ax_grid = plt.subplots(2,
                              1,
                              sharex=True)
  
  axes = {'energy' : {'l':ax_grid[0],
                      'b':ax_grid[0],
                      'r':ax_grid[0].twinx(),
                      't':ax_grid[0].twiny()},
          'entropy': {'l':ax_grid[1],
                      'b':ax_grid[1],
                      'r':ax_grid[1].twinx(),
                      't':ax_grid[1].twiny()} }
  
  xmin_hartree = None
  xmax_hartree = None
  ymin_shannon = None
  ymax_shannon = None
  
  for datum in data:
    axes['energy']['l'].plot(datum['thermal energies'],
                             datum['energies'],
                             linewidth=1,
                             dashes=datum['dashes'],
                             color=colours['turquoise'])
    axes['energy']['l'].plot(datum['thermal energies'],
                             datum['free energies'],
                             linewidth=1,
                             dashes=datum['dashes'],
                             color=colours['orange'])
    axes['entropy']['l'].plot(datum['thermal energies'],
                              datum['entropies'],
                              linewidth=1,
                              dashes=datum['dashes'],
                              color=colours['purple'])
    if xmin_hartree == None:
      xmin_hartree = datum['thermal energies'][0]
      xmax_hartree = datum['thermal energies'][-1]
      ymin_shannon = min(datum['entropies'])
      ymax_shannon = max(datum['entropies'])
    else:
      xmin_hartree = min(xmin_hartree, datum['thermal energies'][0])
      xmax_hartree = max(xmax_hartree, datum['thermal energies'][-1])
      ymin_shannon = min(ymin_shannon, min(datum['entropies']))
      ymax_shannon = max(ymax_shannon, max(datum['entropies']))
  
  kb_in_ev_per_k = 8.6173303e-5
  ev_per_hartree = 27.21138602
  kb_in_au = kb_in_ev_per_k / ev_per_hartree
  
  xmin_kelvin  = xmin_hartree/kb_in_au
  xmax_kelvin  = xmax_hartree/kb_in_au
  axes['energy']['b'].set_xlim(xmin_hartree,xmax_hartree)
  axes['energy']['t'].set_xlim(xmin_kelvin,xmax_kelvin)
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
  
  plt.show()

main()
