import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats
import numpy as np
import operator
from matplotlib import rc

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
  file_name = 'harmonic_observables/thermodynamic_variables.dat'
  variables_file = [line.rstrip('\n').split() for line in open(file_name)]
  thermal_energies = []
  energies = []
  free_energies = []
  entropies = []
  for line in variables_file[1:]:
    thermal_energies.append(float(line[0]))
    energies.append(float(line[1]))
    free_energies.append(float(line[2]))
    entropies.append(float(line[3]))
  
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
  axes['energy']['l'].plot(thermal_energies,
                           energies,
                           linewidth=2,
                           color=colours['turquoise'])
  axes['energy']['l'].plot(thermal_energies,
                           free_energies,
                           linewidth=2,
                           color=colours['orange'])
  axes['entropy']['l'].plot(thermal_energies,
                            entropies,
                            linewidth=2,
                            color=colours['purple'])
  
  kb_in_ev_per_k = 8.6173303e-5
  ev_per_rydberg = 13.605693009
  kb_in_au = kb_in_ev_per_k / ev_per_rydberg
  
  xmin_hartree = thermal_energies[0]
  xmax_hartree = thermal_energies[-1]
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
  
  ymin_shannon = min(entropies)
  ymax_shannon = max(entropies)
  ymin_gibbs   = ymin_shannon * kb_in_au
  ymax_gibbs   = ymax_shannon * kb_in_au
  axes['entropy']['l'].set_ylim(ymin_shannon,ymax_shannon)
  axes['entropy']['l'].set_ylabel(r"Vibrational Shannon Entropy per cell, $S/k_B$, (arb. units)")
  axes['entropy']['r'].set_ylim(ymin_gibbs,ymax_gibbs)
  axes['entropy']['r'].set_ylabel(r"Vibrational Gibbs Entropy per cell, $S$, (Ha per K)")
  
  plt.show()

main()
