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
  file_name = 'effective_frequencies.dat'
  frequencies_file = [line.rstrip('\n').split() for line in open(file_name)]
  
  # Split file into modes.
  modes = []
  for line in frequencies_file:
    if len(line)>0:
      if line[0]=='Mode':
        modes.append({'ID':int(line[2])})
      elif line[0]=='Harmonic':
        modes[-1]['Harmonic frequency'] = float(line[2])
      elif line[0]=='Effective':
        modes[-1]['Effective frequency'] = float(line[2])
      elif line[0]=='Displacement':
        modes[-1]['Displacements'] = []
        modes[-1]['Anharmonic energies'] = []
        modes[-1]['Harmonic energies'] = []
        modes[-1]['Effective energies'] = []
        if len(line)==13:
          modes[-1]['Sampled energies'] = []
      else:
        modes[-1]['Displacements'].append(float(line[0]))
        modes[-1]['Anharmonic energies'].append(float(line[1]))
        modes[-1]['Harmonic energies'].append(float(line[2]))
        modes[-1]['Effective energies'].append(float(line[3]))
        if 'Sampled energies' in modes[-1]:
          modes[-1]['Sampled energies'].append(float(line[4]))
  
  for mode in modes:
    mode['Anharmonic difference'] = []
    for har,anh in zip(mode['Harmonic energies'],mode['Anharmonic energies']):
       mode['Anharmonic difference'].append(anh-har)
    
    mode['Effective difference'] = []
    for har,eff in zip(mode['Harmonic energies'],mode['Effective energies']):
       mode['Effective difference'].append(eff-har)
    
    if 'Sampled energies' in mode:
      mode['Sampled difference'] = []
      for har,eff in zip(mode['Harmonic energies'],mode['Sampled energies']):
         mode['Sampled difference'].append(eff-har)
  
  # Plot data.
  
  fig, axes = plt.subplots(3,len(modes))
  
  for ax in axes[0]:
    for mode in modes:
      ax.hlines([mode['Harmonic frequency']], 0, 1, color=colours['turquoise'])
  
  for mode,ax in zip(modes,axes[0]):
      ax.hlines([mode['Harmonic frequency']], 0, 1, color=colours['orange'])
  
  for mode,ax in zip(modes,axes[1]):
    ax.plot(mode['Displacements'], mode['Anharmonic energies'], 
            color=colours['turquoise'])
    ax.plot(mode['Displacements'], mode['Harmonic energies'],
            color=colours['orange'])
    ax.plot(mode['Displacements'], mode['Effective energies'],
            color=colours['purple'])
    if 'Sampled energies' in mode:
      ax.plot(mode['Displacements'], mode['Sampled energies'],
              color=colours['green'])
  
  for mode,ax in zip(modes,axes[2]):
    ax.plot(mode['Displacements'], mode['Anharmonic difference'], 
            color=colours['turquoise'])
    ax.plot(mode['Displacements'], mode['Effective difference'],
            color=colours['purple'])
    if 'Sampled difference' in mode:
      ax.plot(mode['Displacements'], mode['Sampled difference'],
              color=colours['green'])
  
  # Configure top axes.
  min_frequency = modes[0]['Harmonic frequency']
  max_frequency = modes[-1]['Harmonic frequency']
  
  ymin = min_frequency - 0.1*(max_frequency-min_frequency)
  ymax = max_frequency + 0.1*(max_frequency-min_frequency)
  for ax in axes[0]:
    ax.set_ylim(ymin,ymax)
  
  # Configure middle y-axes.
  min_energy = 0
  max_energy = 0
  for mode in modes:
    min_energy = min(min_energy, min(mode['Anharmonic energies']))
    min_energy = min(min_energy, min(mode['Harmonic energies']))
    min_energy = min(min_energy, min(mode['Effective energies']))
    max_energy = max(max_energy, max(mode['Anharmonic energies']))
    max_energy = max(max_energy, max(mode['Harmonic energies']))
    max_energy = max(max_energy, max(mode['Effective energies']))
    if 'Sampled energies' in mode:
      min_energy = min(min_energy, min(mode['Sampled energies']))
      max_energy = max(max_energy, max(mode['Sampled energies']))
  
  ymin = min_energy - 0.1*(max_energy-min_energy)
  ymax = max_energy + 0.1*(max_energy-min_energy)
  for ax in axes[1]:
    ax.set_ylim(ymin,ymax)
  
  # Configure bottom y-axes.
  min_difference = 0
  max_difference = 0
  for mode in modes:
    min_difference = min(min_difference, min(mode['Anharmonic difference']))
    min_difference = min(min_difference, min(mode['Effective difference']))
    max_difference = max(max_difference, max(mode['Anharmonic difference']))
    max_difference = max(max_difference, max(mode['Effective difference']))
    if 'Sampled difference' in mode:
      min_difference = min(min_difference, min(mode['Sampled difference']))
      max_difference = max(max_difference, max(mode['Sampled difference']))
  
  ymin = min_difference - 0.1*(max_difference-min_difference)
  ymax = max_difference + 0.1*(max_difference-min_difference)
  for ax in axes[2]:
    ax.set_ylim(ymin,ymax)
  
  # Configure middle and bottom x-axes.
  min_displacement = 0
  max_displacement = 0
  for mode in modes:
    min_displacement = min(min_displacement, mode['Displacements'][0])
    max_displacement = max(max_displacement, mode['Displacements'][-1])
  
  xmin = min_displacement - 0.1*(max_displacement-min_displacement)
  xmax = max_displacement + 0.1*(max_displacement-min_displacement)
  for ax in axes[1]:
    ax.set_xlim(xmin,xmax)
  for ax in axes[2]:
    ax.set_xlim(xmin,xmax)
  
  plt.show()

main()
