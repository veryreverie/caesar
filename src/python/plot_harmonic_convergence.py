import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.ticker as tick
import scipy.stats
import numpy as np
import operator
from matplotlib import rc
import os.path
import sys

from utils import read_file, split_into_sections

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
  
  # Read in and parse convergence file.
  file_name = 'convergence.dat'
  convergence_file = read_file(file_name)
  data = split_into_sections(convergence_file)
  
  if data[0][0][0]=='No.':
    no_atoms = int(data[0][0][3])
    data = data[1:]
  
  x_axes = []
  y_axes = []
  
  for i,entry in enumerate(data):
    data[i] = {}
    data[i]['header'] = entry[0]
    data[i]['x'] = [float(line[0]) for line in entry[1:]]
    data[i]['y'] = [[float(x) for x in line[1:]] for line in entry[1:]]
    data[i]['y'] = list(zip(*data[i]['y']))
    data[i]['converged'] = [line[-1] for line in data[i]['y']]
    data[i]['dy'] = [[y-line[-1] for y in line] for line in data[i]['y']]
    
    if data[i]['header'][0]=='k-point':
      data[i]['x type'] = 'k-point'
      data[i]['x'] = [1/x for x in data[i]['x']]
    elif data[i]['header'][0]=='q-point':
      data[i]['x type'] = 'q-point'
      data[i]['x'] = [1/x for x in data[i]['x']]
    elif data[i]['header'][0]=='Cutoff':
      data[i]['x type'] = 'cutoff'
    elif data[i]['header'][0]=='Electronic':
      data[i]['x type'] = 'smearing'
    else:
      raise ValueError('Unexpected header: '+' '.join(entry[0]))
    
    if data[i]['header'][4]=='Mode':
      data[i]['y type'] = 'frequencies'
      no_modes = len(entry[1])-1
      no_atoms = no_modes/3
    elif data[i]['header'][4]=='Free':
      data[i]['y type'] = 'free energies'
    else:
      print('Error: unexpected header:')
      print(data[i]['header'])
      sys.exit()
    
    if data[i]['x type'] not in x_axes:
      x_axes.append(data[i]['x type'])
    
    if data[i]['y type'] not in y_axes:
      y_axes.append(data[i]['y type'])
  
  # Read settings file.
  file_name = 'converge_harmonic_frequencies.used_settings'
  settings_file = [line.rstrip('\n').split() for line in open(file_name)]
  for line in settings_file:
    if line[0]=='min_temperature':
      min_temperature = float(line[1])
    if line[0]=='max_temperature':
      max_temperature = float(line[1])
    if line[0]=='no_temperature_steps':
      no_temperature_steps = int(line[1])
    if line[0]=='energy_tolerance':
      energy_tolerance = float(line[1])
  
  # Plot data.
  width_ratios = [x for _ in x_axes for x in [15,1]]
  fig, ax_grid = plt.subplots(len(y_axes),
                              2*len(x_axes),
                              gridspec_kw={'width_ratios':width_ratios},
                              squeeze=False)
  ax_grid = [[[x,y] for x,y in zip(ax[::2],ax[1::2])] for ax in ax_grid]
  
  for i,entry in enumerate(data):
    data[i]['axis'] = ax_grid[y_axes.index(entry['y type'])]\
                             [x_axes.index(entry['x type'])]\
                             [0]
    data[i]['cbar'] = ax_grid[y_axes.index(entry['y type'])]\
                             [x_axes.index(entry['x type'])]\
                             [1]
  
  cmap = plt.get_cmap('viridis')
  
  for entry in data:
    ax = entry['axis']
    cb = entry['cbar']
    
    min_x = min(entry['x'])
    max_x = max(entry['x'])
    
    min_y = min(entry['converged'])
    max_y = max(entry['converged'])
    
    # Plot data.
    for y,conv in zip(entry['dy'],entry['converged']):
      colour = cmap((conv-min_y)/(max_y-min_y))
      ax.plot(entry['x'],y,color=colour)
    
    # Configure x axis.
    if entry['x type']=='cutoff':
      ax.set_xlabel(r'cut-off energy (Ha)')
    elif entry['x type']=='k-point':
      ax.set_xlabel(r'1/$k$-point spacing (bohr)')
      
      ax2 = ax.twiny()
      ax2.set_xlabel(r'$k$-point spacing (bohr$^{-1}$)')
      ax2.set_xticks(entry['x'])
      round_to_2 = lambda x,n: round(x, n-1-int(np.floor(np.log10(abs(x)))))
      ax2.set_xticklabels([round_to_2(1/x,2) for x in entry['x']])
      ax2.set_xlim(min_x, max_x)
    elif entry['x type']=='q-point':
      ax.set_xlabel(r'1/$q$-point spacing (bohr)')
      
      ax2 = ax.twiny()
      ax2.set_xlabel(r'$q$-point spacing (bohr$^{-1}$)')
      ax2.set_xticks(entry['x'])
      round_to_2 = lambda x,n: round(x, n-1-int(np.floor(np.log10(abs(x)))))
      ax2.set_xticklabels([round_to_2(1/x,2) for x in entry['x']])
      ax2.set_xlim(min_x, max_x)
    elif entry['x type']=='smearing':
      ax.set_xlabel(r'electronic smearing (Ha)')
    ax.set_xlim(entry['x'][0],entry['x'][-1])
    
    # Configure y axis.
    ev_per_hartree = 27.2114
    if entry['y type']=='frequencies':
      ax.set_ylabel(r'error in mode frequencies (Ha)')
      cb.set_ylabel(r'Converged mode frequency (Ha)')
      ax.set_ylim(-1e-4,1e-4)
      ax.hlines([1e-3/(3*ev_per_hartree),-1e-3/(3*ev_per_hartree)],min_x,max_x)
    elif entry['y type']=='free energies':
      ax.set_ylabel(r'error in free energies (Ha per primitive cell)')
      cb.set_ylabel(r'Converged free energy (Ha per primitive cell)')
      ax.set_ylim(-1e-3,1e-3)
      ax.hlines([no_atoms*1e-3/ev_per_hartree,-no_atoms*1e-3/ev_per_hartree],min_x,max_x)
    
    # Add colourbar.
    norm = colors.Normalize(vmin=min_y,vmax=max_y)
    colorbar.ColorbarBase(cb, cmap=cmap, norm=norm)
    cb.yaxis.set_ticks_position('left')
  
  # Plot energy convergence.
  fig, ax_grid = plt.subplots(1,2,squeeze=False)
  
  temperatures = np.linspace(min_temperature,max_temperature,no_temperature_steps)
  cmap = plt.get_cmap('inferno')
  temp_colours = []
  for temperature in temperatures:
    fraction = 0.9*(temperature-temperatures[0])/ \
               (temperatures[-1]-temperatures[0])
    temp_colours.append(cmap(fraction))
  
  for i,entry in enumerate(data):
    if entry['x type']=='cutoff' and entry['y type']=='free energies':
      ax = ax_grid[0][0]
      for i,y in enumerate(entry['dy']):
        ax.plot(entry['x'],y,color=temp_colours[i])
      ax.set_yscale('symlog', linthreshy=energy_tolerance)
      ax.hlines(energy_tolerance,entry['x'][0],entry['x'][-1],linestyle='dotted')
      ax.hlines(-energy_tolerance,entry['x'][0],entry['x'][-1],linestyle='dotted')
    elif entry['x type']=='k-point' and entry['y type']=='free energies':
      ax = ax_grid[0][1]
      for i,y in enumerate(entry['dy']):
        ax.plot(entry['x'],y,color=temp_colours[i])
      ax.set_yscale('symlog', linthreshy=energy_tolerance)
      ax.hlines(energy_tolerance,entry['x'][0],entry['x'][-1],linestyle='dotted')
      ax.hlines(-energy_tolerance,entry['x'][0],entry['x'][-1],linestyle='dotted')
  
  plt.show()

if __name__=='__main__':
  main()
