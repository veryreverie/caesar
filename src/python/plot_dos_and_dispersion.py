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
    'Orange'   :[252/255,141/255, 98/255],
    'blue'     :[141/255,160/255,203/255],
    'purple'   :[231/255,138/255,195/255],
    'green'    :[166/255,216/255, 84/255],
    'yellow'   :[255/255,217/255, 47/255],
    'beige'    :[229/255,196/255,148/255],
    'grey'     :[179/255,179/255,179/255]}
  
  # Read in high symmetry points.
  file_name = 'high_symmetry_points.dat'
  points_file = [line.rstrip('\n').split() for line in open(file_name)]
  points = {'labels':[], 'q-points':[], 'path_lengths':[]}
  for line in points_file:
    if len(line)>0 and line[0]=='q-point:':
      points['labels'].append(line[1])
      points['q-points'].append([])
      for coord in line[2:]:
        points['q-points'][-1].append(float(coord))
    elif len(line)>0 and line[0]=='Fraction':
      points['path_lengths'].append(float(line[3]))
  
  # Replace the label 'G' with Gamma
  for i,label in enumerate(points['labels']):
    if label=='G':
      points['labels'][i] = r'\Gamma'
  
  # Read in phonon dispersion.
  file_name = 'phonon_dispersion_curve.dat'
  dispersion_file = [line.rstrip('\n').split() for line in open(file_name)]
  dispersion = {'path_length':[], 'frequencies':[], 'bands':[]}
  for line in dispersion_file:
    if line[0]=='Fraction':
      dispersion['path_length'].append(float(line[3]))
    elif line[0]=='Frequencies:':
      dispersion['frequencies'].append([])
      for frequency in line[1:]:
        dispersion['frequencies'][-1].append(float(frequency))
  
  # Identify bands.
  band_order = []
  for i,frequency in enumerate(dispersion['frequencies'][0]):
    dispersion['bands'].append([])
    band_order.append(i)
  
  for frequencies in dispersion['frequencies']:
    if len(dispersion['bands'][0])<2:
      # If less than two points in each band,
      #    assign to bands in order of frequency.
      for i,frequency in enumerate(frequencies):
        dispersion['bands'][i].append(frequency)
    else:
      # Assign to bands in terms of the order of a linear extrapolation
      #    of the last two points in each band.
      next_point = []
      band_order = []
      for i,band in enumerate(dispersion['bands']):
        next_point.append(2*band[-1]-band[-2])
        band_order.append(i)
      band_order = [x for _,x in sorted(zip(next_point,band_order))]
      for i,frequency in enumerate(frequencies):
        dispersion['bands'][band_order[i]].append(frequency)
    
  #for band in zip(*dispersion['frequencies']):
  #  dispersion['bands'].append(band)
  
  # Read in dos.
  file_name = 'freq_dos.dat'
  dos_file = [line.rstrip('\n').split() for line in open(file_name)]
  dos = {'bottoms':[], 'tops':[], 'middles':[], 'dos':[]}
  for line in dos_file:
    if len(line)>0 and line[0]=='Bin:':
      dos['bottoms'].append(float(line[1]))
      dos['tops'].append(float(line[3]))
      dos['middles'].append((float(line[1])+float(line[3]))/2)
    elif len(line)>0 and line[0]=='Bin':
      dos['dos'].append(float(line[2]))
  
  
  # Plot everything.
  xmin = 0
  xmax = 1
  ymin = dos['bottoms'][0]
  ymax = dos['tops'][-1]
  
  fig, ax_grid = plt.subplots( 1,
                               2,
                               sharey=True, 
                               gridspec_kw={ 'width_ratios':[3,1]})
  axes = {}
  
  axes['dispersion'] = ax_grid[0]
  axes['dispersion'].set_xlim(xmin,xmax)
  axes['dispersion'].set_ylim(ymin,ymax)
  for band in dispersion['bands']:
    axes['dispersion'].plot(dispersion['path_length'],band)
  
  axes['dispersion'].vlines(points['path_lengths'],ymin,ymax,linestyle=':')
  axes['dispersion'].set_xticks(points['path_lengths'])
  axes['dispersion'].set_xticklabels(points['labels'])
  axes['dispersion'].minorticks_off()
  axes['dispersion'].set_ylabel('Energy, Hartrees')
  
  axes['dos'] = ax_grid[1]
  axes['dos'].plot(dos['dos'],dos['middles'])
  axes['dos'].set_xticks([])
  axes['dos'].minorticks_off()
  
  plt.show()

main()