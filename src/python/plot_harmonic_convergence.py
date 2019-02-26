import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
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
      data[i]['x label'] = r'$k$-point spacing (bohr$^{-1}$)'
    elif data[i]['header'][0]=='Cutoff':
      data[i]['x type'] = 'cutoff'
      data[i]['x label'] = r'cut-off energy (Ha)'
    else:
      print(data[i]['header'][0])
      sys.exit()
    
    if data[i]['header'][4]=='Mode':
      data[i]['y type'] = 'frequencies'
      data[i]['y label'] = r'error in mode frequencies (Ha)'
    elif data[i]['header'][4]=='Free':
      data[i]['y type'] = 'free energies'
      data[i]['y label'] = r'error in free energies (Ha per primitive cell)'
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
    
    min_y = min(entry['converged'])
    max_y = max(entry['converged'])
    
    for y,conv in zip(entry['dy'],entry['converged']):
      colour = cmap((conv-min_y)/(max_y-min_y))
      ax.plot(entry['x'],y,color=colour)
    ax.set_xlabel(entry['x label'])
    ax.set_ylabel(entry['y label'])
    ax.set_xlim(entry['x'][0],entry['x'][-1])
    ax.set_ylim(-1e-3,1e-3)
    
    norm = colors.Normalize(vmin=min_y,vmax=max_y)
    colorbar.ColorbarBase(cb, cmap=cmap, norm=norm)
    cb.yaxis.set_ticks_position('left')
  
  plt.show()

if __name__=='__main__':
  main()
