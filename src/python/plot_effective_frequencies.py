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
    'Orange'   :[252/255,141/255, 98/255],
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
      else:
        modes[-1]['Displacements'].append(float(line[0]))
        modes[-1]['Anharmonic energies'].append(float(line[1]))
        modes[-1]['Harmonic energies'].append(float(line[2]))
        modes[-1]['Effective energies'].append(float(line[3]))
  
  fig, axes = plt.subplots(1,len(modes))
  
  for mode,ax in zip(modes,axes):
    ax.plot(mode['Displacements'],modes[0]['Anharmonic energies'])
    ax.plot(mode['Displacements'],modes[0]['Harmonic energies'])
    ax.plot(mode['Displacements'],modes[0]['Effective energies'])
  
  plt.show()

main()
