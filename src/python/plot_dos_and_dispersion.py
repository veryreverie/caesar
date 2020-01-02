import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats
import numpy as np
import operator
from matplotlib import rc
import os.path
import sys

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
  file_name = 'temperatures.dat'
  if os.path.isfile(file_name):
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
  
  # Read files.
  data = []
  plotting_dispersion = False
  plotting_dos = False
  for directory in directories:
    # Read in high symmetry points.
    points = read_high_symmetry_points(directory)
    
    # Read in phonon dispersion.
    dispersion = read_phonon_dispersion(directory)
    
    if dispersion is not None:
      plotting_dispersion = True
    
    # Read in dos.
    dos = read_dos(directory)
    
    if dos is not None:
      plotting_dos = True
    
    data.append({'points':points, 'dispersion':dispersion, 'dos':dos})
  
  # --------------------------------------------------
  # Plot everything.
  # --------------------------------------------------
  # Initialise axes.
  if plotting_dispersion and plotting_dos:
    fig, ax_grid = plt.subplots( 1,
                                 2,
                                 sharey=True, 
                                 gridspec_kw={ 'width_ratios':[3,1]})
    axes = {'dispersion':ax_grid[0], 'dos':ax_grid[1]}
  elif plotting_dispersion:
    fig, ax_grid = plt.subplots(1,1)
    axes = {'dispersion':ax_grid}
  elif plotting_dos:
    fig, ax_grid = plt.subplots(1,1)
    axes = {'dos':ax_grid}
  else:
    print('Error: dos and dispersion files missing.')
    sys.exit()
  
  # Calculate y-axis range.
  ymin = 0
  ymax = 0
  for datum in data:
    if plotting_dos:
      ymin = min(ymin, datum['dos']['bottoms'][0])
      ymax = max(ymax, datum['dos']['tops'][-1])
    else:
      for frequencies in datum['dispersion']['frequencies']:
        ymin = min(ymin, frequencies[0])
        ymax = max(ymax, frequencies[-1])
  
  # Plot dispersion.
  if plotting_dispersion:
    axes['dispersion'].set_xlim(0,1)
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
                                    color=colours['orange'], lw=2, zorder=2)
          
          for segment in positive_segments:
            axes['dispersion'].plot(segment['x'],segment['y'],
                                    color=colours['turquoise'], lw=2, zorder=2)
        else:
          axes['dispersion'].plot(xs,ys,color=temp_colours[i], lw=2, zorder=2)
          
    
    dashed_lines = [y for x,y in zip(data[0]['points']['lines'],data[0]['points']['path_lengths']) if x=='dashed']
    solid_lines = [y for x,y in zip(data[0]['points']['lines'],data[0]['points']['path_lengths']) if x=='solid']
    axes['dispersion'].vlines(dashed_lines, ymin, ymax, linestyle=':')
    axes['dispersion'].vlines(solid_lines, ymin, ymax, zorder=3, lw=2)
    axes['dispersion'].set_xticks(data[0]['points']['path_lengths'])
    axes['dispersion'].set_xticklabels(data[0]['points']['labels'])
    axes['dispersion'].minorticks_off()
    axes['dispersion'].set_ylabel('Energy, Hartrees')
    axes['dispersion'].tick_params(bottom=False)
    
    hartree_to_mev = 2.721138602e4
    axes['ev'] = axes['dispersion'].twinx()
    axes['ev'].set_ylim(ymin*hartree_to_mev, ymax*hartree_to_mev)
    axes['ev'].set_ylabel('Energy, meV')
  
  # Plot DOS.
  if plotting_dos:
    max_dos = 0
    for i,datum in reversed(list(enumerate(data))):
      if len(data)==1:
        # zero is the first index where ['middles'] is > 0.
        zero = next(i for i,e in enumerate(datum['dos']['middles']) if e>0)
        # Calculated ['dos'] at the point w=0.
        midpoint = ( datum['dos']['dos'][zero-1] \
                   * datum['dos']['middles'][zero] \
                   - datum['dos']['dos'][zero] \
                   * datum['dos']['middles'][zero-1] ) \
                 / ( datum['dos']['middles'][zero] \
                   - datum['dos']['middles'][zero-1] )
        # Plot the DOS for w<=0 in orange.
        axes['dos'].plot(datum['dos']['dos'][:zero]+[midpoint],
                         datum['dos']['middles'][:zero]+[0],
                         color=colours['orange'],lw=2)
        # Plot the DOS for w>=0 in turquoise.
        axes['dos'].plot([midpoint]+datum['dos']['dos'][zero:],
                         [0]+datum['dos']['middles'][zero:],
                         color=colours['turquoise'],lw=2)
        max_dos = max(max_dos, max(datum['dos']['dos']))
      else:
        axes['dos'].plot(datum['dos']['dos'],datum['dos']['middles'],
                         color=temp_colours[i],lw=2)
        max_dos = max(max_dos, max(datum['dos']['dos']))
    axes['dos'].set_xticks([])
    axes['dos'].set_xlim(0,1.1*max_dos)
    axes['dos'].minorticks_off()
    
    hartree_to_inverse_cm = 2.194746313702e5
    axes['cm'] = axes['dos'].twinx()
    axes['cm'].set_ylim(ymin*hartree_to_inverse_cm, ymax*hartree_to_inverse_cm)
    axes['cm'].set_ylabel(r'Frequency, cm$^{-1}$')
  
  #plt.tight_layout()
  #plt.savefig('esdg_fig.pdf', bbox_inches='tight')
  
  plt.show()

def read_high_symmetry_points(directory):
  # Read in high symmetry points.
  file_name = directory + '/high_symmetry_points.dat'
  
  if not os.path.isfile(file_name):
    return None
  
  points_file = [line.rstrip('\n').split() for line in open(file_name)]
  points = {'labels':[], 'q-points':[], 'path_lengths':[], 'indices':[]}
  for line in points_file:
    
    # Old file layout.
    if len(line)==5 and line[0]=='q-point:':
      points['labels'].append(line[1])
      points['q-points'].append([float(x) for x in line[2:]])
    elif len(line)==4 and line[:3]==['Fraction','along','path:']:
      points['path_lengths'].append(float(line[3]))
    
    # New file layout.
    if len(line)==4 and line[:2]==['q-point','label']:
      points['labels'].append(line[3])
    elif len(line)==5 and line[0]=='q-point':
      points['q-points'].append([float(x) for x in line[2:]])
    elif len(line)==5 and line[:3]==['Fraction','along','path']:
      points['path_lengths'].append(float(line[4]))
    elif len(line)==5 and line[:3]==['Index','along','path']:
      points['indices'].append(int(line[4]))
  
  # Replace the label 'G' with Gamma
  for i,label in enumerate(points['labels']):
    if label=='G':
      points['labels'][i] = r'\Gamma'
  
  # Identify segment breaks.
  points['lines'] = ['dashed' for _ in points['labels']]
  for i in range(len(points['indices'])-1):
    if points['indices'][i+1] == points['indices'][i]+1:
      points['lines'][i:i+2] = ['solid','solid']
      points['labels'][i] = points['labels'][i]+' '+points['labels'][i+1]
      points['labels'][i+1] = ''
  
  return points

def read_phonon_dispersion(directory):
  # Read in phonon dispersion.
  file_name = directory + '/phonon_dispersion_curve.dat'
  
  if not os.path.isfile(file_name):
    return None
  
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
  
  return dispersion

def read_dos(directory):
  # Read in dos.
  file_name = directory + '/phonon_density_of_states.dat'
  
  if not os.path.isfile(file_name):
    return None
  
  dos_file = [line.rstrip('\n').split() for line in open(file_name)]
  dos = {'bottoms':[], 'tops':[], 'middles':[], 'dos':[]}
  for line in dos_file:
    dos['bottoms'].append(float(line[2]))
    dos['tops'].append(float(line[4]))
    dos['middles'].append((float(line[2])+float(line[4]))/2)
    dos['dos'].append(float(line[6]))
  
  return dos

main()
