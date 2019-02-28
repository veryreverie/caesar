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
    elif data[i]['header'][0]=='Cutoff':
      data[i]['x type'] = 'cutoff'
    else:
      print(data[i]['header'][0])
      sys.exit()
    
    if data[i]['header'][4]=='Mode':
      data[i]['y type'] = 'frequencies'
    elif data[i]['header'][4]=='Free':
      data[i]['y type'] = 'free energies'
    else:
      print('Error: unexpected header:')
      print(data[i]['header'])
      sys.exit()
  
  # Plot data.
  if len(data)==1:
    fig, ax_grid = plt.subplots(1,2,gridspec_kw={'width_ratios':[15,1]})
    data[0]['axis'] = ax_grid[0]
    data[0]['cbar'] = ax_grid[1]
  elif len(data)==2:
    if data[1]['y type']=='free energies':
      fig, ax_grid = plt.subplots(2,2,gridspec_kw={'width_ratios':[15,1]})
      data[0]['axis'] = ax_grid[0][0]
      data[0]['cbar'] = ax_grid[0][1]
      data[1]['axis'] = ax_grid[1][0]
      data[1]['cbar'] = ax_grid[1][1]
    else:
      fig, ax_grid = plt.subplots(1,4,gridspec_kw={'width_ratios':[15,1,15,1]})
      data[0]['axis'] = ax_grid[0]
      data[0]['cbar'] = ax_grid[1]
      data[1]['axis'] = ax_grid[2]
      data[1]['cbar'] = ax_grid[3]
  elif len(data)==4:
    fig, ax_grid = plt.subplots(2,4,gridspec_kw={'width_ratios':[15,1,15,1]})
    data[0]['axis'] = ax_grid[0][0]
    data[0]['cbar'] = ax_grid[0][1]
    data[1]['axis'] = ax_grid[1][0]
    data[1]['cbar'] = ax_grid[1][1]
    data[2]['axis'] = ax_grid[0][2]
    data[2]['cbar'] = ax_grid[0][3]
    data[3]['axis'] = ax_grid[1][2]
    data[3]['cbar'] = ax_grid[1][3]
  else:
    print('Error: '+len(data)+' entries.')
    sys.exit()
  
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
    ax.set_xlim(entry['x'][0],entry['x'][-1])
    
    # Configure y axis.
    if entry['y type']=='frequencies':
      ax.set_ylabel(r'error in mode frequencies (Ha)')
      cb.set_ylabel(r'Converged mode frequency (Ha)')
      ax.set_ylim(-1e-3,1e-3)
    elif entry['y type']=='free energies':
      ax.set_ylabel(r'error in free energies (Ha per primitive cell)')
      cb.set_ylabel(r'Converged free energy (Ha per primitive cell)')
      ax.set_ylim(-1e-2,1e-2)
    
    # Add colourbar.
    norm = colors.Normalize(vmin=min_y,vmax=max_y)
    colorbar.ColorbarBase(cb, cmap=cmap, norm=norm)
    cb.yaxis.set_ticks_position('left')
  
  plt.show()

if __name__=='__main__':
  main()
