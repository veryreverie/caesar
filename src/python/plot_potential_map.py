
import matplotlib.pyplot as plt
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
  
  # Plot data.
  fig = plt.subplots(squeeze=False)
  axes = [[plt.subplot2grid((2,3),(0,0)),
           plt.subplot2grid((2,3),(0,1)),
           plt.subplot2grid((2,3),(0,2))],
          [plt.subplot2grid((2,3),(1,0)),
           plt.subplot2grid((2,3),(1,1))]]
  
  axes[0][0].contour(x_displacements,y_displacements,harmonic_energy)
  axes[0][1].contour(x_displacements,y_displacements,anharmonic_energy)
  axes[0][2].contour(x_displacements,y_displacements,sampled_energy)
  
  axes[1][0].contour(x_displacements,y_displacements,harmonic_difference)
  axes[1][1].contour(x_displacements,y_displacements,anharmonic_difference)
  
  # Configure axes.  
  for row in axes:
    for ax in row:
      ax.set_aspect('equal')
  
  plt.show()

main()
