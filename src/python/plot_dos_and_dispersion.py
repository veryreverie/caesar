import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats
import numpy as np
import operator
from matplotlib import rc
import os.path

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
  
  # Check if a temperatures.dat file exists.
  temperatures = []
  directories = []
  if os.path.isfile('temperatures.dat'):
    file_name = 'temperatures.dat'
    temperature_file = [line.rstrip('\n').split() for line in open(file_name)]
    for i,line in enumerate(temperature_file[1:]):
      temperatures.append(float(line[0]))
      directories.append('temperature_' + \
                         str(i+1).zfill(len(str(len(temperature_file)))))
    cmap = plt.get_cmap('inferno')
    temp_colours = []
    for temperature in temperatures:
      fraction = 0.9*(temperature-temperatures[0])/ \
                 (temperatures[-1]-temperatures[0])
      temp_colours.append(cmap(fraction))
  else:
    temperatures = [0.0]
    directories = ['.']
  
  data = []
  for directory in directories:
    # Read in high symmetry points.
    file_name = directory + '/high_symmetry_points.dat'
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
    file_name = directory + '/phonon_dispersion_curve.dat'
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
    
    # Read in dos.
    file_name = directory + '/phonon_density_of_states.dat'
    dos_file = [line.rstrip('\n').split() for line in open(file_name)]
    dos = {'bottoms':[], 'tops':[], 'middles':[], 'dos':[]}
    for line in dos_file:
      dos['bottoms'].append(float(line[2]))
      dos['tops'].append(float(line[4]))
      dos['middles'].append((float(line[2])+float(line[4]))/2)
      dos['dos'].append(float(line[6]))
    
    data.append({'points':points, 'dispersion':dispersion, 'dos':dos})
  
  # Plot everything.
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 0
  for datum in data:
    ymin = min(ymin, datum['dos']['bottoms'][0])
    ymax = max(ymax, datum['dos']['tops'][-1])
  
  fig, ax_grid = plt.subplots( 1,
                               2,
                               sharey=True, 
                               gridspec_kw={ 'width_ratios':[3,1]})
  axes = {}
  
  axes['dispersion'] = ax_grid[0]
  axes['dispersion'].set_xlim(xmin,xmax)
  axes['dispersion'].set_ylim(ymin,ymax)
  for i,datum in reversed(list(enumerate(data))):
    xs = datum['dispersion']['path_length']
    for band in datum['dispersion']['bands']:
      ys = band
      
      if len(data)==1:
        # Plot each band segment by segment.
        positive_segments = []
        negative_segments = []
        if ys[0]>=0:
          positive_segments.append({'x':xs[:1], 'y':ys[:1]})
        else:
          negative_segments.append({'x':xs[:1], 'y':ys[:1]})
        for x1,y1,x2,y2 in zip(xs[:-1],ys[:-1],xs[2:],ys[2:]):
          if y1>=0:
            if y2>=0:
              positive_segments[-1]['x'].append(x2)
              positive_segments[-1]['y'].append(y2)
            else:
              x_mid = (x2*y1-x1*y2)/(y1-y2)
              positive_segments[-1]['x'].append(x_mid)
              positive_segments[-1]['y'].append(0)
              negative_segments.append({'x':[x_mid,x2], 'y':[0,y2]})
          else:
            if y2>=0:
              x_mid = (x2*y1-x1*y2)/(y1-y2)
              negative_segments[-1]['x'].append(x_mid)
              negative_segments[-1]['y'].append(0)
              positive_segments.append({'x':[x_mid,x2], 'y':[0,y2]})
            else:
              negative_segments[-1]['x'].append(x2)
              negative_segments[-1]['y'].append(y2)
        
        for segment in negative_segments:
          axes['dispersion'].plot(segment['x'],segment['y'],
                                  color=colours['orange'], lw=2)
        
        for segment in positive_segments:
          axes['dispersion'].plot(segment['x'],segment['y'],
                                  color=colours['turquoise'], lw=2)
      else:
        axes['dispersion'].plot(xs,ys,color=temp_colours[i], lw=2)
        
  
  axes['dispersion'].vlines(data[0]['points']['path_lengths'],
                            ymin,
                            ymax,
                            linestyle=':')
  axes['dispersion'].set_xticks(data[0]['points']['path_lengths'])
  axes['dispersion'].set_xticklabels(data[0]['points']['labels'])
  axes['dispersion'].minorticks_off()
  axes['dispersion'].set_ylabel('Energy, Hartrees')
  
  hartree_to_mev = 2.721138602e4
  axes['ev'] = axes['dispersion'].twinx()
  axes['ev'].set_ylim(ymin*hartree_to_mev, ymax*hartree_to_mev)
  axes['ev'].set_ylabel('Energy, meV')
  
  
  
  axes['dos'] = ax_grid[1]
  for i,datum in reversed(list(enumerate(data))):
    if len(data)==1:
      zero = next(i for i,e in enumerate(datum['dos']['middles']) if e>0)
      axes['dos'].plot(datum['dos']['dos'][:zero],
                       datum['dos']['middles'][:zero],
                       color=colours['orange'],lw=2)
      axes['dos'].plot(datum['dos']['dos'][zero:],
                       datum['dos']['middles'][zero:],
                       color=colours['turquoise'],lw=2)
    else:
      axes['dos'].plot(datum['dos']['dos'],datum['dos']['middles'],
                       color=temp_colours[i],lw=2)
  axes['dos'].set_xticks([])
  axes['dos'].minorticks_off()
  
  hartree_to_inverse_cm = 2.194746313702e5
  axes['cm'] = axes['dos'].twinx()
  axes['cm'].set_ylim(ymin*hartree_to_inverse_cm, ymax*hartree_to_inverse_cm)
  axes['cm'].set_ylabel(r'Frequency, cm$^{-1}$')
  
  plt.show()

main()
