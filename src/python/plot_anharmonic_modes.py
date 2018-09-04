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
  file_name = 'mode_maps.dat'
  frequencies_file = [line.rstrip('\n').split() for line in open(file_name)]
  
  sampling = False
  
  # Split file into modes.
  modes = []
  for line in frequencies_file:
    if len(line)>0:
      if line[0]=='Mode':
        modes.append({'ID':int(line[2])})
      elif line[0]=='Harmonic':
        modes[-1]['Harmonic frequency'] = float(line[2])
      elif line[0]=='Displacement':
        modes[-1]['Displacements'] = []
        modes[-1]['Harmonic energies'] = []
        modes[-1]['Harmonic forces'] = []
        modes[-1]['Anharmonic cos energies'] = []
        modes[-1]['Anharmonic cos forces'] = []
        modes[-1]['Anharmonic sin energies'] = []
        modes[-1]['Anharmonic sin forces'] = []
        if line[-3]=='Sampled':
          sampling = True
          modes[-1]['Sampled cos energies'] = []
          modes[-1]['Sampled cos forces'] = []
          modes[-1]['Sampled sin energies'] = []
          modes[-1]['Sampled sin forces'] = []
      else:
        modes[-1]['Displacements'].append(float(line[0]))
        modes[-1]['Harmonic energies'].append(float(line[1]))
        modes[-1]['Harmonic forces'].append(float(line[2]))
        modes[-1]['Anharmonic cos energies'].append(float(line[3]))
        modes[-1]['Anharmonic cos forces'].append(float(line[4]))
        modes[-1]['Anharmonic sin energies'].append(float(line[5]))
        modes[-1]['Anharmonic sin forces'].append(float(line[6]))
        if sampling:
          modes[-1]['Sampled cos energies'].append(float(line[7]))
          modes[-1]['Sampled cos forces'].append(float(line[8]))
          modes[-1]['Sampled sin energies'].append(float(line[9]))
          modes[-1]['Sampled sin forces'].append(float(line[10]))
  
  if sampling:
    for mode in modes:
      mode['Harmonic cos energy difference'] = []
      for harmonic,sampled in zip(mode['Harmonic energies'],
                                  mode['Sampled cos energies']):
         mode['Harmonic cos energy difference'].append(harmonic-sampled)
      
      mode['Harmonic cos force difference'] = []
      for harmonic,sampled in zip(mode['Harmonic forces'],
                                  mode['Sampled cos forces']):
         mode['Harmonic cos force difference'].append(harmonic-sampled)
      
      mode['Harmonic sin energy difference'] = []
      for harmonic,sampled in zip(mode['Harmonic energies'],
                                  mode['Sampled sin energies']):
         mode['Harmonic sin energy difference'].append(harmonic-sampled)
      
      mode['Harmonic sin force difference'] = []
      for harmonic,sampled in zip(mode['Harmonic forces'],
                                  mode['Sampled sin forces']):
         mode['Harmonic sin force difference'].append(harmonic-sampled)
      
      mode['Anharmonic cos energy difference'] = []
      for anharmonic,sampled in zip(mode['Anharmonic cos energies'],
                                    mode['Sampled cos energies']):
         mode['Anharmonic cos energy difference'].append(anharmonic-sampled)
      
      mode['Anharmonic cos force difference'] = []
      for anharmonic,sampled in zip(mode['Anharmonic cos forces'],
                                    mode['Sampled cos forces']):
         mode['Anharmonic cos force difference'].append(anharmonic-sampled)
      
      mode['Anharmonic sin energy difference'] = []
      for anharmonic,sampled in zip(mode['Anharmonic sin energies'],
                                    mode['Sampled sin energies']):
         mode['Anharmonic sin energy difference'].append(anharmonic-sampled)
      
      mode['Anharmonic sin force difference'] = []
      for anharmonic,sampled in zip(mode['Anharmonic sin forces'],
                                    mode['Sampled sin forces']):
         mode['Anharmonic sin force difference'].append(anharmonic-sampled)
  
  # Plot data.
  
  if sampling:
    fig, axes = plt.subplots(5,len(modes))
  else:
    fig, axes = plt.subplots(3,len(modes))
  
  for ax in axes[0]:
    for mode in modes:
      ax.hlines([mode['Harmonic frequency']], 0, 1, color=colours['turquoise'])
  for mode,ax in zip(modes,axes[0]):
      ax.hlines([mode['Harmonic frequency']], 0, 1, color=colours['orange'])
      ax.set_xlabel('Mode '+str(mode['ID']))
      ax.xaxis.set_label_position('top')
  
  for mode,ax in zip(modes,axes[1]):
    ax.plot(mode['Displacements'], mode['Harmonic energies'],
            color=colours['orange'], lw=2)
    ax.plot(mode['Displacements'], mode['Anharmonic cos energies'], 
            color=colours['turquoise'], lw=2)
    ax.plot(mode['Displacements'], mode['Anharmonic sin energies'], 
            color=colours['turquoise'], lw=2, linestyle='dashed')
    if sampling:
      ax.plot(mode['Displacements'], mode['Sampled cos energies'],
              color=colours['green'], lw=2)
      ax.plot(mode['Displacements'], mode['Sampled sin energies'],
              color=colours['green'], lw=2, linestyle='dashed')
  
  for mode,ax in zip(modes,axes[2]):
    ax.plot(mode['Displacements'], mode['Harmonic forces'],
            color=colours['orange'], lw=2)
    ax.plot(mode['Displacements'], mode['Anharmonic cos forces'], 
            color=colours['turquoise'], lw=2)
    ax.plot(mode['Displacements'], mode['Anharmonic sin forces'], 
            color=colours['turquoise'], lw=2, linestyle='dashed')
    if sampling:
      ax.plot(mode['Displacements'], mode['Sampled cos forces'],
              color=colours['green'], lw=2)
      ax.plot(mode['Displacements'], mode['Sampled sin forces'],
              color=colours['green'], lw=2, linestyle='dashed')
  
  if sampling:
    for mode,ax in zip(modes,axes[3]):
      ax.plot(mode['Displacements'], mode['Harmonic cos energy difference'],
              color=colours['orange'], lw=2)
      ax.plot(mode['Displacements'], mode['Harmonic sin energy difference'],
              color=colours['orange'], lw=2, linestyle='dashed')
      ax.plot(mode['Displacements'], mode['Anharmonic cos energy difference'], 
              color=colours['turquoise'], lw=2)
      ax.plot(mode['Displacements'], mode['Anharmonic sin energy difference'], 
              color=colours['turquoise'], lw=2, linestyle='dashed')
  
  if sampling:
    for mode,ax in zip(modes,axes[4]):
      ax.plot(mode['Displacements'], mode['Harmonic cos force difference'],
              color=colours['orange'], lw=2)
      ax.plot(mode['Displacements'], mode['Harmonic sin force difference'],
              color=colours['orange'], lw=2, linestyle='dashed')
      ax.plot(mode['Displacements'], mode['Anharmonic cos force difference'], 
              color=colours['turquoise'], lw=2)
      ax.plot(mode['Displacements'], mode['Anharmonic sin force difference'], 
              color=colours['turquoise'], lw=2, linestyle='dashed')
  
  # Configure frequency axes.
  min_frequency = modes[0]['Harmonic frequency']
  max_frequency = modes[-1]['Harmonic frequency']
  
  ymin = min_frequency - 0.1*(max_frequency-min_frequency)
  ymax = max_frequency + 0.1*(max_frequency-min_frequency)
  
  ymin = min(ymin, 0)
  ymax = max(ymax, 0)
  
  for ax in axes[0]:
    ax.set_ylim(ymin,ymax)
    ax.set_xticks([])
  axes[0][0].tick_params(direction='out')
  axes[0][0].set_ylabel('Energy (Ha)')
  for ax in axes[0][1:]:
    ax.tick_params(direction='out', labelleft='off')
  
  # Configure energy y-axes.
  min_energy = 0
  max_energy = 0
  for mode in modes:
    min_energy = min(min_energy, min(mode['Anharmonic cos energies']))
    min_energy = min(min_energy, min(mode['Anharmonic sin energies']))
    max_energy = max(max_energy, max(mode['Anharmonic cos energies']))
    max_energy = max(max_energy, max(mode['Anharmonic sin energies']))
    if sampling:
      min_energy = min(min_energy, min(mode['Sampled cos energies']))
      min_energy = min(min_energy, min(mode['Sampled sin energies']))
      max_energy = max(max_energy, max(mode['Sampled cos energies']))
      max_energy = max(max_energy, max(mode['Sampled sin energies']))
  
  ymin = min_energy - 0.1*(max_energy-min_energy)
  ymax = max_energy + 0.1*(max_energy-min_energy)
  for ax in axes[1]:
    ax.set_ylim(ymin,ymax)
  for ax in axes[1][1:]:
    ax.tick_params(labelleft='off')
  axes[1][0].set_ylabel('Energy (Ha)')
  
  # Configure force y-axes.
  min_force = 0
  max_force = 0
  for mode in modes:
    min_force = min(min_force, min(mode['Anharmonic cos forces']))
    min_force = min(min_force, min(mode['Anharmonic sin forces']))
    max_force = max(max_force, max(mode['Anharmonic cos forces']))
    max_force = max(max_force, max(mode['Anharmonic sin forces']))
    if sampling:
      min_force = min(min_force, min(mode['Sampled cos forces']))
      min_force = min(min_force, min(mode['Sampled sin forces']))
      max_force = max(max_force, max(mode['Sampled cos forces']))
      max_force = max(max_force, max(mode['Sampled sin forces']))
  
  ymin = min_force - 0.1*(max_force-min_force)
  ymax = max_force + 0.1*(max_force-min_force)
  for ax in axes[2]:
    ax.set_ylim(ymin,ymax)
  for ax in axes[2][1:]:
    ax.tick_params(labelleft='off')
  axes[2][0].set_ylabel('Force (Ha/Bohr)')
  
  # Configure energy difference y-axes.
  if sampling:
    min_difference = 0
    max_difference = 0
    for mode in modes:
      min_difference = min(min_difference,
                           min(mode['Anharmonic cos energy difference']))
      min_difference = min(min_difference,
                           min(mode['Anharmonic sin energy difference']))
      max_difference = max(max_difference,
                           max(mode['Anharmonic cos energy difference']))
      max_difference = max(max_difference,
                           max(mode['Anharmonic sin energy difference']))
    
    ymin = min_difference - 0.1*(max_difference-min_difference)
    ymax = max_difference + 0.1*(max_difference-min_difference)
    ymin = min(ymin,-ymax)
    ymax = max(ymax,-ymin)
    for ax in axes[3]:
      ax.set_ylim(ymin,ymax)
    for ax in axes[3][1:]:
      ax.tick_params(labelleft='off')
    axes[3][0].set_ylabel('Energy difference (Ha)')
  
  # Configure force difference y-axes.
  if sampling:
    min_difference = 0
    max_difference = 0
    for mode in modes:
      min_difference = min(min_difference,
                           min(mode['Anharmonic cos force difference']))
      min_difference = min(min_difference,
                           min(mode['Anharmonic sin force difference']))
      max_difference = max(max_difference,
                           max(mode['Anharmonic cos force difference']))
      max_difference = max(max_difference,
                           max(mode['Anharmonic sin force difference']))
    
    ymin = min_difference - 0.1*(max_difference-min_difference)
    ymax = max_difference + 0.1*(max_difference-min_difference)
    ymin = min(ymin,-ymax)
    ymax = max(ymax,-ymin)
    for ax in axes[4]:
      ax.set_ylim(ymin,ymax)
    for ax in axes[4][1:]:
      ax.tick_params(labelleft='off')
    axes[4][0].set_ylabel('Force difference (Ha/Bohr)')
  
  # Configure displacement x-axes.
  min_displacement = 0
  max_displacement = 0
  for mode in modes:
    min_displacement = min(min_displacement, mode['Displacements'][0])
    max_displacement = max(max_displacement, mode['Displacements'][-1])
  
  xmin = min_displacement - 0.1*(max_displacement-min_displacement)
  xmax = max_displacement + 0.1*(max_displacement-min_displacement)
  for axis in axes[1:]:
    for ax in axis:
      ax.set_xlim(xmin,xmax)
      ax.tick_params(direction='out', labelbottom='off')
  
  if sampling:
    for ax in axes[4]:
      ax.tick_params(labelbottom='on')
  else:
    for ax in axes[2]:
      ax.tick_params(labelbottom='on')
  
  plt.show()

main()
