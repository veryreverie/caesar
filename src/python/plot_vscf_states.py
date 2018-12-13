import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats
import numpy as np
import operator
from matplotlib import rc
import matplotlib.lines as mlines
import os.path

from utils import read_file,split_into_sections

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
  
  # Read in anharmonic data.
  filename = 'anharmonic_data.dat'
  anharmonic_data = read_file(filename)
  for i,line in enumerate(anharmonic_data):
    if line==['Maximum','weighted','displacement']:
      maximum_displacement = float(anharmonic_data[i+1][0])
    elif line==['Frequency','of','maximum','displacement']:
      min_frequency = float(anharmonic_data[i+1][0])
  
  # Read in degenerate subspaces.
  filename = 'degenerate_subspaces.dat'
  subspaces_file = split_into_sections(read_file(filename))
  subspaces = []
  for section in subspaces_file:
    subspaces.append({'ID'         : int(section[0][4]),
                      'Frequency'  : float(section[1][4]),
                      'Mode IDs'   : [int(x) for x in section[2][4:]],
                      'Paired IDs' : [int(x) for x in section[3][4:]]
                      })
  
  max_subspace_id = max([subspace['ID'] for subspace in subspaces])
  len_max_subspace_id = len(str(max_subspace_id))
  format_string = '{0:0'+str(len_max_subspace_id)+'}'
  
  # For each subspace, read in vscf_wavefunctions file.
  single_mode_wavefunctions = []
  for subspace in subspaces:
    padded_id = format_string.format(subspace['ID'])
    filename = 'anharmonic_observables/subspace_'+padded_id+'/vscf_wavefunctions.dat'
    wavefunctions_file = split_into_sections(read_file(filename))
    
    # Parse the header.
    header = {'ID'         : int(wavefunctions_file[0][0][2]),
              'Mode IDs'   : [int(x) for x in wavefunctions_file[0][1][3:]],
              'Paired IDs' : [int(x) for x in wavefunctions_file[0][2][4:]],
              '|0>'        : wavefunctions_file[0][3][5]
             }
    
    # Parse the wavefunctions.
    wavefunctions = []
    for section in wavefunctions_file[1:]:
      wavefunctions.append({'Energy'       : float(section[0][3]),
                            'Degeneracy'   : int(section[1][3]),
                            'Wavefunction' : ' '.join(section[2][2:])
                           })
    
    # Pick out the ground state.
    ground_state = wavefunctions[np.argmin([x['Energy'] for x in wavefunctions])]
    
    # Replace '^' and '|0>' with python-friendly '**' and 'ground_state'.
    header['|0>'] = header['|0>'].replace('^','**')
    ground_state['Wavefunction'] = ground_state['Wavefunction'].replace('^','**')
    ground_state['Wavefunction'] = ground_state['Wavefunction'].replace('|0>','*ground_state')
    
    # Calculate displacements.
    max_mode_displacement = maximum_displacement * min_frequency / max(min_frequency, subspace['Frequency'])
    
    # Add the value of e to the keys for eval(),
    #    and initialise each mode to zero.
    keys = {'e' : np.e}
    for mode_id in header['Mode IDs']:
      keys['u'+str(mode_id)] = float(0)
    
    # Loop over modes, calculating the wavefunction along each.
    for mode_id in header['Mode IDs']:
      displacements = np.linspace(-max_mode_displacement,
                                  max_mode_displacement,
                                  50
                                 ).tolist()
      single_mode_wavefunction = []
      for displacement in displacements:
        # Set the displacement along the mode.
        keys['u'+str(mode_id)] = displacement
        # Evaluate the harmonic ground state.
        keys['ground_state'] = eval(header['|0>'],keys)
        # Evaluate the ground state wavefunction.
        single_mode_wavefunction.append(eval(ground_state['Wavefunction'],keys))
      single_mode_wavefunctions.append({'Subspace'      : subspace['ID'],
                                        'Mode'          : mode_id,
                                        'Displacements' : displacements,
                                        'Wavefunction'  : single_mode_wavefunction
                                       })
  
  # Plot everything.
  fig, ax_grid = plt.subplots(1,len(single_mode_wavefunctions))
  for ax,mode in zip(ax_grid,single_mode_wavefunctions):
    print('')
    print(mode['Displacements'])
    print(mode['Wavefunction'])
    print(mode)
    ax.plot(mode['Displacements'],mode['Wavefunction'])
  
  plt.show()
  
main()
