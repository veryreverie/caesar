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
  
  # Read file.
  file_name = 'anharmonic_mode_maps.dat'
  
  frequencies, sampling, using_stress, modes = parse_mode_maps(file_name)
  
  keys = ['frequency',
          'energies',
          'forces']
  if using_stress:
    keys.append('pressures')
  if sampling:
    keys.append('energy differences')
    keys.append('force differences')
    if using_stress:
      keys.append('pressure differences')
  
  # Plot data.
  fig, axs = plt.subplots(len(keys), len(modes), squeeze=False)
  axes = {}
  for key,ax in zip(keys,axs):
    axes[key] = ax
  
  for ax in axes['frequency']:
    for frequency in frequencies:
      ax.hlines([frequency], 0, 1, color=colours['turquoise'])
  for mode,ax in zip(modes,axes['frequency']):
    ax.hlines([mode['Harmonic frequency']], 0, 1, color=colours['orange'])
    ax.set_xlabel('Mode '+str(mode['ID']))
    ax.xaxis.set_label_position('top')
  
  for mode,ax in zip(modes,axes['energies']):
    ax.plot(mode['Displacements'], mode['Harmonic energies'],
            color=colours['orange'], lw=2)
    ax.plot(mode['Displacements'], mode['Anharmonic energies'], 
            color=colours['turquoise'], lw=2)
    if sampling:
      ax.plot(mode['Displacements'], mode['Sampled energies'],
              color=colours['green'], lw=2)
  
  for mode,ax in zip(modes,axes['forces']):
    ax.plot(mode['Displacements'], mode['Harmonic forces'],
            color=colours['orange'], lw=2)
    ax.plot(mode['Displacements'], mode['Anharmonic forces'], 
            color=colours['turquoise'], lw=2)
    if sampling:
      ax.plot(mode['Displacements'], mode['Sampled forces'],
              color=colours['green'], lw=2)
  
  if using_stress:
    for mode,ax in zip(modes,axes['pressures']):
      ax.plot(mode['Displacements'], mode['Anharmonic pressures'], 
              color=colours['turquoise'], lw=2)
      if sampling:
        ax.plot(mode['Displacements'], mode['Sampled pressures'],
                color=colours['green'], lw=2)
  
  if sampling:
    for mode,ax in zip(modes,axes['energy differences']):
      ax.plot(mode['Displacements'], mode['Harmonic energy difference'],
              color=colours['orange'], lw=2)
      ax.plot(mode['Displacements'], mode['Anharmonic energy difference'], 
              color=colours['turquoise'], lw=2)
    
    for mode,ax in zip(modes,axes['force differences']):
      ax.plot(mode['Displacements'], mode['Harmonic force difference'],
              color=colours['orange'], lw=2)
      ax.plot(mode['Displacements'], mode['Anharmonic force difference'], 
              color=colours['turquoise'], lw=2)
    
    if using_stress:
      for mode,ax in zip(modes,axes['pressure differences']):
        ax.plot(mode['Displacements'], mode['Anharmonic pressure difference'], 
                color=colours['turquoise'], lw=2)
  
  # Configure frequency axes.
  min_frequency = min(min(frequencies), 0)
  max_frequency = max(max(frequencies), 0)
  
  ymin = min_frequency - 0.1*(max_frequency-min_frequency)
  ymax = max_frequency + 0.1*(max_frequency-min_frequency)
  
  for ax in axes['frequency']:
    ax.set_ylim(ymin,ymax)
    ax.set_xticks([])
  axes['frequency'][0].tick_params(direction='out')
  axes['frequency'][0].set_ylabel('Energy (Ha)')
  for ax in axes['frequency'][1:]:
    ax.tick_params(direction='out', labelleft='off')
  
  # Configure energy y-axes.
  min_energy = 0
  max_energy = 0
  for mode in modes:
    min_energy = min(min_energy, min(mode['Anharmonic energies']))
    max_energy = max(max_energy, max(mode['Anharmonic energies']))
    if sampling:
      min_energy = min(min_energy, min(mode['Sampled energies']))
      max_energy = max(max_energy, max(mode['Sampled energies']))
  
  ymin = min_energy - 0.1*(max_energy-min_energy)
  ymax = max_energy + 0.1*(max_energy-min_energy)
  for ax in axes['energies']:
    ax.set_ylim(ymin,ymax)
  for ax in axes['energies'][1:]:
    ax.tick_params(labelleft='off')
  axes['energies'][0].set_ylabel('Energy (Ha)')
  
  # Configure force y-axes.
  min_force = 0
  max_force = 0
  for mode in modes:
    min_force = min(min_force, min(mode['Anharmonic forces']))
    max_force = max(max_force, max(mode['Anharmonic forces']))
    if sampling:
      min_force = min(min_force, min(mode['Sampled forces']))
      max_force = max(max_force, max(mode['Sampled forces']))
  
  ymin = min_force - 0.1*(max_force-min_force)
  ymax = max_force + 0.1*(max_force-min_force)
  for ax in axes['forces']:
    ax.set_ylim(ymin,ymax)
  for ax in axes['forces'][1:]:
    ax.tick_params(labelleft='off')
  axes['forces'][0].set_ylabel('Force (Ha/Bohr)')
  
  if sampling:
    # Configure energy difference y-axes.
    min_difference = 0
    max_difference = 0
    for mode in modes:
      min_difference = min(min_difference,
                           min(mode['Anharmonic energy difference']))
      max_difference = max(max_difference,
                           max(mode['Anharmonic energy difference']))
    
    ymin = min_difference - 0.1*(max_difference-min_difference)
    ymax = max_difference + 0.1*(max_difference-min_difference)
    ymin = min(ymin,-ymax)
    ymax = max(ymax,-ymin)
    for ax in axes['energy differences']:
      ax.set_ylim(ymin,ymax)
    for ax in axes['energy differences'][1:]:
      ax.tick_params(labelleft='off')
    axes['energy differences'][0].set_ylabel('Energy difference (Ha)')
    
    # Configure force difference y-axes.
    min_difference = 0
    max_difference = 0
    for mode in modes:
      min_difference = min(min_difference,
                           min(mode['Anharmonic force difference']))
      max_difference = max(max_difference,
                           max(mode['Anharmonic force difference']))
    
    ymin = min_difference - 0.1*(max_difference-min_difference)
    ymax = max_difference + 0.1*(max_difference-min_difference)
    ymin = min(ymin,-ymax)
    ymax = max(ymax,-ymin)
    for ax in axes['force differences']:
      ax.set_ylim(ymin,ymax)
    for ax in axes['force differences'][1:]:
      ax.tick_params(labelleft='off')
    axes['force differences'][0].set_ylabel('Force difference (Ha/Bohr)')
  
  # Configure displacement x-axes.
  min_displacement = 0
  max_displacement = 0
  for mode in modes:
    min_displacement = min(min_displacement, mode['Displacements'][0])
    max_displacement = max(max_displacement, mode['Displacements'][-1])
  
  xmin = min_displacement - 0.1*(max_displacement-min_displacement)
  xmax = max_displacement + 0.1*(max_displacement-min_displacement)
  for key,axis in axes.items():
    if key=='frequency':
      continue
    for ax in axis:
      ax.set_xlim(xmin,xmax)
      ax.tick_params(direction='out', labelbottom='off')
  
  for ax in axes[keys[-1]]:
    ax.tick_params(labelbottom='on')
  
  plt.show()

main()
