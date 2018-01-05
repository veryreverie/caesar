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
  
  fig, axes = plt.subplots(4, len(modes))
  
  min_frequency = modes[0]['frequency']
  max_frequency = modes[-1]['frequency']
  
  max_displacement = 0
  for mode in modes:
    for displacement in mode['displacements']:
      for number in displacement:
        max_displacement = max(max_displacement, np.absolute(number))
  radius = max_displacement/10
  lim = max_displacement * 1.2
  
  eqm_pos = []
  for i in range(no_atoms):
    eqm_pos.append(lim*(2*(no_atoms-i)-1))
  
  for i,mode in enumerate(modes):
    # Plot frequencies.
    for mode2 in modes:
      axes[0][i].hlines([mode2['frequency']],0,1,color=colours['turquoise'])
    axes[0][i].hlines([mode['frequency']],0,1,color=colours['orange'])
    axes[0][i].set_xlim(0,1)
    axes[0][i].set_ylim(min_frequency-0.001,max_frequency+0.001)
    
    axes[0][i].set_xticks([])
    axes[0][i].minorticks_off()
    axes[0][i].tick_params( direction='out', 
                            labelleft='off')
    
    for dirn in range(1,4):
      axes[dirn][i].set_xlim(-lim,lim)
      axes[dirn][i].set_ylim(0,no_atoms*lim*2)
      axes[dirn][i].set_xticks([])
      axes[dirn][i].set_aspect('equal')
      if i>0:
        axes[dirn][i].set_yticks([])
    
    # Plot displacements.
    for j,displacement in enumerate(mode['displacements']):
      for dirn in range(1,4):
        eqm = pch.Circle( (0,eqm_pos[j]),
                          radius=radius,
                          fill=False,
                          edgecolor=colours['turquoise'])
        axes[dirn][i].add_patch(eqm)
      
      no_plots = 12
      for k in range(no_plots):
        phase = complex(0,np.pi*k/(3*no_plots))
        x = (displacement[0]*np.exp(phase)).real
        y = (displacement[1]*np.exp(phase)).real
        z = (displacement[2]*np.exp(phase)).real
        scale = (max_displacement*2)/(max_displacement*2-z)
        circ = pch.Circle( (x,y+eqm_pos[j]),
                           radius=radius,
                           color=colours['orange'])
        axes[1][i].add_patch(circ)
        circ = pch.Circle( (y,z+eqm_pos[j]),
                           radius=radius,
                           color=colours['orange'])
        axes[2][i].add_patch(circ)
        circ = pch.Circle( (z,x+eqm_pos[j]),
                           radius=radius,
                           color=colours['orange'])
        axes[3][i].add_patch(circ)
        
  axes[0][0].set_ylabel('Energy, Hartrees')
  axes[0][0].tick_params(labelleft='on')
  axes[1][0].set_ylabel('x-y elevation')
  axes[2][0].set_ylabel('y-z elevation')
  axes[3][0].set_ylabel('z-x elevation')
  
  for i in [1,2,3]:
    axes[i][0].tick_params( direction='out', 
                            labelleft='on',
                            labelright='off')
    axes[i][0].set_yticks(eqm_pos)
    axes[i][0].set_yticklabels(species)
    axes[i][0].tick_params(length=0)
  
  hartree_to_inverse_cm = 2.194746313702e5
  rax = axes[0][-1].twinx()
  axes[0][-1].tick_params( direction='out', 
                           labelleft='off')
  rax.tick_params( direction='out', 
                   labelleft='off')
  rax.set_ylim( (min_frequency-0.001)*hartree_to_inverse_cm, 
                (max_frequency+0.001)*hartree_to_inverse_cm)
  rax.set_ylabel(r'Frequency, cm$^{-1}$')
  
  #axes['frequencies'].set_xlabel('Mode id')
  #axes['frequencies'].set_ylabel('Energy, Hartrees')
  
  plt.show()

main()
