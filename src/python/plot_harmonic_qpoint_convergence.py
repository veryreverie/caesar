import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats
import numpy as np
import operator
from matplotlib import rc
import os.path
import sys

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
  
  # Read convergence.dat.
  file_name = 'convergence.dat'
  convergence_file = [line.rstrip('\n').split() for line in open(file_name)]
  
  # Parse the convergence file.
  threshold = float(convergence_file[2][4])
  min_temperature = float(convergence_file[3][4])
  max_temperature = float(convergence_file[3][5])
  qpoint_spacings = []
  free_energies = []
  for line in convergence_file[6:]:
    qpoint_spacings.append(float(line[0]))
    free_energies.append([float(f) for f in line[1:]])
  temperatures = np.linspace(min_temperature,max_temperature,len(free_energies[0]))
  
  # Generate colour map.
  cmap = plt.get_cmap('inferno')
  temp_colours = []
  for temperature in temperatures:
    fraction = 0.9*(temperature-temperatures[0])/ \
               (temperatures[-1]-temperatures[0])
    temp_colours.append(cmap(fraction))
  
  # Plot data.
  fig, ax_grid = plt.subplots(1, 1, squeeze=False)
  for i in range(len(temperatures)):
    ax_grid[0][0].plot(qpoint_spacings,
                       [f[i]-free_energies[-1][i] for f in free_energies],
                       color=temp_colours[i])
  
  # Configure axes.
  ax_grid[0][0].set_xlim((qpoint_spacings[0],qpoint_spacings[-1]))
  ax_grid[0][0].set_xlabel(r'q-point spacing (bohr$^{-1}$)')
  ax_grid[0][0].set_ylabel(r'free energy difference (Ha per cell)')
  ax_grid[0][0].set_yscale('symlog', linthreshy=threshold)
  ax_grid[0][0].hlines(threshold,qpoint_spacings[0],qpoint_spacings[-1],linestyle='dotted')
  ax_grid[0][0].hlines(-threshold,qpoint_spacings[0],qpoint_spacings[-1],linestyle='dotted')
  plt.show()

main()
