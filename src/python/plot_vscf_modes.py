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
  file_name = 'vscf_mode_maps.dat'
  map_file = [line.rstrip('\n').split() for line in open(file_name)]
  
  displacements = [float(x) for x in map_file[0][1:]]
  modes = []
  for line in map_file[1:]:
    modes.append({'ID':int(line[1]), 'Potential':[float(x) for x in line[2:]]})
  
  # Read anharmonic modes file.
  file_name = 'mode_maps.dat'
  anharmonic_modes_file = [line.rstrip('\n').split() for line in open(file_name)]
  anharmonic_modes = []
  for line in anharmonic_modes_file[2:]:
    if len(line)>0:
      if line[0]=='Mode':
        anharmonic_modes.append({'ID':int(line[2])})
      elif line[0]=='Harmonic':
        anharmonic_modes[-1]['Harmonic frequency'] = float(line[2])
      elif line[0]=='Displacement':
        anharmonic_modes[-1]['Displacements'] = []
        anharmonic_modes[-1]['Harmonic energies'] = []
        anharmonic_modes[-1]['Harmonic forces'] = []
        anharmonic_modes[-1]['Anharmonic cos energies'] = []
        anharmonic_modes[-1]['Anharmonic cos forces'] = []
        anharmonic_modes[-1]['Anharmonic sin energies'] = []
        anharmonic_modes[-1]['Anharmonic sin forces'] = []
        if line[-3]=='Sampled':
          sampling = True
          anharmonic_modes[-1]['Sampled cos energies'] = []
          anharmonic_modes[-1]['Sampled cos forces'] = []
          anharmonic_modes[-1]['Sampled sin energies'] = []
          anharmonic_modes[-1]['Sampled sin forces'] = []
      else:
        anharmonic_modes[-1]['Displacements'].append(float(line[0]))
        anharmonic_modes[-1]['Harmonic energies'].append(float(line[1]))
        anharmonic_modes[-1]['Harmonic forces'].append(float(line[2]))
        anharmonic_modes[-1]['Anharmonic cos energies'].append(float(line[3]))
        anharmonic_modes[-1]['Anharmonic cos forces'].append(float(line[4]))
        anharmonic_modes[-1]['Anharmonic sin energies'].append(float(line[5]))
        anharmonic_modes[-1]['Anharmonic sin forces'].append(float(line[6]))
        if sampling:
          anharmonic_modes[-1]['Sampled cos energies'].append(float(line[7]))
          anharmonic_modes[-1]['Sampled cos forces'].append(float(line[8]))
          anharmonic_modes[-1]['Sampled sin energies'].append(float(line[9]))
          anharmonic_modes[-1]['Sampled sin forces'].append(float(line[10]))
  
  # Plot data.
  
  fig, axes = plt.subplots(1,len(modes))
  if len(modes)==1:
    axes = [axes]
  for mode,anharmonic_mode,ax in zip(modes,anharmonic_modes,axes):
    ax.set_xlabel('Mode '+str(mode['ID']))
    ax.xaxis.set_label_position('top')
    ax.plot(displacements, mode['Potential'])
    ax.plot(anharmonic_mode['Displacements'],
            anharmonic_mode['Anharmonic cos energies'])
  
  plt.show()

main()
