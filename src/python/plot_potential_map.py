
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
import operator
from matplotlib import rc, gridspec

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
  
  # Read file.
  file_name = 'potential.dat'
  potential_file = [line.rstrip('\n').split() for line in open(file_name)]
  
  x_id = int(potential_file[0][3])
  y_id = int(potential_file[1][3])
  x_displacements = [float(d) for d in potential_file[3]]
  y_displacements = [float(d) for d in potential_file[5]]
  
  no_displacements = len(x_displacements)
  
  harmonic_energy = [[float(x) for x in line] for line in potential_file[7:7+no_displacements]]
  anharmonic_energy = [[float(x) for x in line] for line in potential_file[8+no_displacements:8+2*no_displacements]]
  if len(potential_file)==9+3*no_displacements:
    sampled_energy = [[float(x) for x in line] for line in potential_file[9+2*no_displacements:9+3*no_displacements]]
  
  harmonic_difference = [[x-y for x,y in zip(xs,ys)] for xs,ys in zip(harmonic_energy,sampled_energy)]
  anharmonic_difference = [[x-y for x,y in zip(xs,ys)] for xs,ys in zip(anharmonic_energy,sampled_energy)]
  
  # Calculate data ranges.
  min_energy = 0
  max_energy = 0
  min_difference = 0
  max_difference = 0
  
  min_energy = min(min_energy, np.min(harmonic_energy))
  min_energy = min(min_energy, np.min(anharmonic_energy))
  min_energy = min(min_energy, np.min(sampled_energy))
  max_energy = max(max_energy, np.max(harmonic_energy))
  max_energy = max(max_energy, np.max(anharmonic_energy))
  max_energy = max(max_energy, np.max(sampled_energy))
  
  min_difference = min(min_difference, np.min(harmonic_difference))
  min_difference = min(min_difference, np.min(anharmonic_difference))
  max_difference = max(max_difference, np.max(harmonic_difference))
  max_difference = max(max_difference, np.max(anharmonic_difference))
  
  energy_levels = np.linspace(min_energy,max_energy,50)
  if max_difference>abs(min_difference):
    difference_levels = np.linspace(-max_difference,max_difference,50)
  else:
    difference_levels = np.linspace(min_difference,-min_difference,50)
  
  # Plot data.
  fig,axes = plt.subplots(squeeze=False)
  gs = gridspec.GridSpec(2, 3)
  axes = [[plt.subplot(gs[0,0]),
           plt.subplot(gs[0,1]),
           plt.subplot(gs[0,2])],
          [plt.subplot(gs[1,0]),
           plt.subplot(gs[1,1])]]
  
  eh = axes[0][0].contourf(x_displacements,y_displacements,harmonic_energy,
                           levels=energy_levels,cmap=plt.cm.viridis)
  ea = axes[0][1].contourf(x_displacements,y_displacements,anharmonic_energy,
                           levels=energy_levels,cmap=plt.cm.viridis)
  ec = axes[0][2].contourf(x_displacements,y_displacements,sampled_energy,
                           levels=energy_levels,cmap=plt.cm.viridis)
  
  dh = axes[1][0].contourf(x_displacements,y_displacements,harmonic_difference,
                           levels=difference_levels,cmap=plt.cm.bwr)
  da = axes[1][1].contourf(x_displacements,y_displacements,
                           anharmonic_difference, levels=difference_levels,
                           cmap=plt.cm.bwr)
  
  # Configure axes.  
  for row in axes:
    for ax in row:
      ax.set_aspect('equal')
  
  # Add colourbars.
  ecb_ax = fig.add_axes([0.6,0.05,0.1,0.45])
  ecb = plt.colorbar(eh, ax=ecb_ax, orientation='vertical')
  ecb_ax.set_axis_off()
  
  dcb_ax = fig.add_axes([0.71,0.05,0.1,0.45])
  dcb = plt.colorbar(dh, ax=dcb_ax, orientation='vertical')
  dcb_ax.set_axis_off()
  
  # Label axes.
  axes[0][0].set_title(r'Harmonic energy $E_h$ (Ha per cell)')
  axes[0][1].set_title(r'Anharmonic energy $E_a$ (Ha per cell)')
  axes[0][2].set_title(r'Sampled energy $E_s$ (Ha per cell)')
  
  axes[1][0].set_title(r'$E_h-E_s$')
  axes[1][1].set_title(r'$E_a-E_s$')
  
  plt.show()

main()
