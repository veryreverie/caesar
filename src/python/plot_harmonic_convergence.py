import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
    data[i]['y'] = zip(*[[float(x) for x in line[1:]] for line in entry[1:]])
    data[i]['dy'] = [[y-line[-1] for y in line] for line in data[i]['y']]
    
    if data[i]['header'][0]=='k-point':
      data[i]['x label'] = r'$k$-point spacing (bohr$^{-1}$)'
    elif data[i]['header'][0]=='Cutoff':
      data[i]['x label'] = r'cut-off energy (Ha)'
    else:
      print(data[i]['header'][0])
      sys.exit()
    
    if data[i]['header'][4]=='Mode':
      data[i]['y label'] = r'error in mode frequencies (Ha)'
    elif data[i]['header'][0]=='Free':
      data[i]['y label'] = r'error in free energies (Ha per primitive cell)'
    else:
      print(data[i]['header'][4])
      sys.exit()
  
  # Plot data.
  fig, ax_grid = plt.subplots(1, len(data))
  
  if len(data)==1:
    ax_grid = [ax_grid]
  
  for entry, ax in zip(data, ax_grid):
    for y in entry['dy']:
      ax.plot(entry['x'],y)
    ax.set_xlabel(entry['x label'])
    ax.set_ylabel(entry['y label'])
    ax.set_xlim(entry['x'][0],entry['x'][-1])
    ax.set_ylim(-1e-4,1e-4)
  
  plt.show()
  

if __name__=='__main__':
  main()
