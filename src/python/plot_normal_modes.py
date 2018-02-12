import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
import operator
from matplotlib import rc
import matplotlib.patches as pch
from matplotlib.collections import PatchCollection
import os

rc('font', **{'family':'serif','serif':['sffamily']})
rc('text', usetex=True)
params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

def dblecomplex(displacement):
  if displacement[0]=='-':
    re = float(displacement[:24])
    im = float(displacement[24:-1])
  else:
    re = float(displacement[:23])
    im = float(displacement[23:-1])
  return complex(re,im)

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
  
  filename = '../structure.dat'
  contents = [line.rstrip('\n').split() for line in open(filename)]
  reading_species = False
  species = []
  for line in contents:
    if len(line)>0 and line[0]=='Atoms':
      reading_species = True
    elif reading_species:
      if len(line)==1:
        reading_species = False
      else:
        species.append(line[0])
  
  modes = []
  for filename in sorted(os.listdir()):
    if filename[:5]=='mode_':
      mode = int(filename[5:-4])
      contents = [line.rstrip('\n').split() for line in open(filename)]
      reading_displacements = False
      for line in contents:
        if len(line)>0 and line[0]=='Frequency:':
          frequency = float(line[1])
          modes.append({ 'mode':mode,
                         'frequency':frequency,
                         'displacements':[]})
        elif len(line)>0 and line[0]=='Primitive':
          reading_displacements = not reading_displacements
        elif reading_displacements:
          modes[-1]['displacements'].append([])
          for displacement in line:
            modes[-1]['displacements'][-1].append(dblecomplex(displacement))
  
  no_atoms = len(modes[0]['displacements'])
  
  fig, axes_numbered = plt.subplots(4, len(modes))
  axes = {'f':axes_numbered[0],
          'xy':axes_numbered[1],
          'yz':axes_numbered[2],
          'zx':axes_numbered[3]}
  
  min_frequency = modes[0]['frequency']
  max_frequency = modes[-1]['frequency']
  
  max_displacement = 0
  for mode in modes:
    for displacement in mode['displacements']:
      for number in displacement:
        max_displacement = max(max_displacement, np.absolute(number))
  radius = max_displacement/10
  lim = max_displacement * 1.5
  
  eqm_pos = []
  for i in range(no_atoms):
    eqm_pos.append(lim*(2*(no_atoms-i)-1))
  
  for i,mode in enumerate(modes):
    # Plot frequencies.
    for mode2 in modes:
      axes['f'][i].hlines([mode2['frequency']],0,1,color=colours['turquoise'])
    axes['f'][i].hlines([mode['frequency']],0,1,color=colours['orange'])
    axes['f'][i].set_xlim(0,1)
    axes['f'][i].set_ylim(min_frequency-0.001,max_frequency+0.001)
    
    axes['f'][i].set_xticks([])
    axes['f'][i].minorticks_off()
    axes['f'][i].tick_params( direction='out', 
                            labelleft='off')
    
    for dirn in ['xy','yz','zx']:
      axes[dirn][i].set_xlim(-lim,lim)
      axes[dirn][i].set_ylim(0,no_atoms*lim*2)
      axes[dirn][i].set_xticks([])
      axes[dirn][i].set_aspect('equal')
      if i>0:
        axes[dirn][i].set_yticks([])
    
    # Plot displacements.
    for j,displacement in enumerate(mode['displacements']):
      # Calculate positions along path.
      no_plots = 24
      mid = no_plots//2
      points = {'xy':{'h':[],'v':[]},
                'yz':{'h':[],'v':[]},
                'zx':{'h':[],'v':[]}}
      for k in range(no_plots):
        phase = complex(0,np.pi*k*2/no_plots)
        x = (displacement[0]*np.exp(phase)).real
        y = (displacement[1]*np.exp(phase)).real
        z = (displacement[2]*np.exp(phase)).real
        scale = (max_displacement*2)/(max_displacement*2-z)
        points['xy']['h'].append(x)
        points['xy']['v'].append(y+eqm_pos[j])
        points['yz']['h'].append(y)
        points['yz']['v'].append(z+eqm_pos[j])
        points['zx']['h'].append(z)
        points['zx']['v'].append(x+eqm_pos[j])
      
      for k in ['xy','yz','zx']:
        # Plot equilibrium position.
        axes[k][i].add_patch(pch.Circle((0,eqm_pos[j]),
                             radius=radius*2,
                             fill=False,
                             edgecolor=colours['turquoise']))
        # Plot lines.
        axes[k][i].plot(points[k]['h'],points[k]['v'],color=colours['orange'])
        
        if (abs(points[k]['h'][0]-points[k]['h'][mid])>0.5*radius or
            abs(points[k]['v'][0]-points[k]['v'][mid])>0.5*radius):
          # Plot arrows.
          axes[k][i].arrow(points[k]['h'][-1],
                           points[k]['v'][-1],
                           points[k]['h'][0]-points[k]['h'][-1],
                           points[k]['v'][0]-points[k]['v'][-1],
                           width=0,
                           head_width=0.3,
                           head_length=0.15,
                           overhang=1.5,
                           length_includes_head=True,
                           color=colours['orange'])
        else:
          # Plot circles.
          axes[k][i].add_patch(pch.Circle((points[k]['h'][0],
                                           points[k]['v'][0]),
                                          radius=radius,
                                          color=colours['orange']))
        
        
  axes['f'][0].set_ylabel('Energy, Hartrees')
  axes['f'][0].tick_params(labelleft='on')
  axes['xy'][0].set_ylabel('x-y elevation')
  axes['yz'][0].set_ylabel('y-z elevation')
  axes['zx'][0].set_ylabel('z-x elevation')
  
  for i in ['xy','yz','zx']:
    axes[i][0].tick_params( direction='out', 
                            labelleft='on',
                            labelright='off')
    axes[i][0].set_yticks(eqm_pos)
    axes[i][0].set_yticklabels(species)
    axes[i][0].tick_params(length=0)
  
  hartree_to_inverse_cm = 2.194746313702e5
  rax = axes['f'][-1].twinx()
  axes['f'][-1].tick_params( direction='out', 
                             labelleft='off')
  rax.tick_params( direction='out', 
                   labelleft='off')
  rax.set_ylim( (min_frequency-0.001)*hartree_to_inverse_cm, 
                (max_frequency+0.001)*hartree_to_inverse_cm)
  rax.set_ylabel(r'Frequency, cm$^{-1}$')
  
  plt.show()

main()
