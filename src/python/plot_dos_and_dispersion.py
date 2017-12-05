
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
  
  fig = plt.figure()
  axes = {}
  
  # Read in high symmetry points.
  file_name = 'high_symmetry_points.dat'
  points_file = [line.rstrip('\n').split() for line in open(file_name)]
  points = []
  for line in points_file:
    if (len(line)>0 and line[0]=='q-point:'):
      points.append({'q-point':[]})
      for coord in line[1:]:
        points[-1]['q-point'].append(float(coord))
    elif (len(line)>0 and line[0]=='Fraction'):
      points[-1]['path_length'] = float(line[3])
  
  # Read in phonon dispersion.
  file_name = 'phonon_dispersion_curve.dat'
  dispersion_file = [line.rstrip('\n').split() for line in open(file_name)]
  dispersion = {'path_length':[], 'frequencies':[], 'bands':[]}
  for line in dispersion_file:
    if (line[0]=='Fraction'):
      dispersion['path_length'].append(float(line[3]))
    elif (line[0]=='Frequencies:'):
      dispersion['frequencies'].append([])
      for frequency in line[1:]:
        dispersion['frequencies'][-1].append(float(frequency))
  
  for band in zip(*dispersion['frequencies']):
    dispersion['bands'].append(band)
  
  # Read in dos.
  file_name = 'freq_dos.dat'
  
  # Plot everything.
  xmin = 0
  xmax = 1
  ymin = min(dispersion['bands'][0])
  ymax = max(dispersion['bands'][-1])
  
  axes['dispersion'] = fig.add_subplot(111)
  axes['dispersion'].set_xlim(xmin,xmax)
  axes['dispersion'].set_ylim(ymin,ymax)
  for band in dispersion['bands']:
    axes['dispersion'].plot(dispersion['path_length'],band)
  
  for point in points:
    axes['dispersion'].vlines(point['path_length'],ymin,ymax)
  
  plt.show()

main()
