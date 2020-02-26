import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
import operator
from matplotlib import rc

from utils import parse_mode_maps

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
  
  # Read vscf modes file.
  file_name = 'vscf_mode_maps.dat'
  _, anharmonic, _, _, vscf_modes = parse_mode_maps(file_name)
  if not anharmonic:
    sys.exit()
  
  # Read anharmonic modes file.
  file_name = 'anharmonic_mode_maps.dat'
  frequencies, anharmonic, sampled, pressure, anharmonic_modes = parse_mode_maps(file_name)
  if not anharmonic:
    sys.exit()
  
  # Plot data.
  
  fig, axes = plt.subplots(1,len(vscf_modes))
  if len(vscf_modes)==1:
    axes = [axes]
  for vscf_mode,anharmonic_mode,ax in zip(vscf_modes,anharmonic_modes,axes):
    ax.set_xlabel('Mode '+str(vscf_mode['ID']))
    ax.xaxis.set_label_position('top')
    ax.plot(vscf_mode['Displacements'],
            vscf_mode['Anharmonic energies'])
    ax.plot(anharmonic_mode['Displacements'],
            anharmonic_mode['Anharmonic energies'])
  
  plt.show()

main()
