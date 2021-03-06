#!/usr/bin/python3
import os
import sys
import numpy as np

import ase
import ase.io

from ase.calculators.gulp import GULP
from ase.units import Ha, Bohr

def main(directory,seedname):
  '''
  Read an extended .xyz file, run it through Gulp, and produce an
     electronic_structure.dat file.
  N.B. this requires ASE 3.19 or later, as earlier versions have problems
     reading forces from Gulp.
  N.B. the environment variable GULP_LIB must point to the Libraries directory
     in Gulp.
  '''
  
  # Set up the calculator.
  # This can be replaced with any other ASE calculator.
  # N.B. this reads from Caesar's root working directory,
  #   so only one settings file is needed.
  keywords, library, options, has_forces, has_stress = read_gulp_file('{0}.gin'.format(seedname))
  calculator = GULP(keywords=keywords,library=library,options=options)
  
  # Change to the directory where calculations should be run.
  os.chdir(directory)
  
  # Read the atomic positions and lattice from the extended .xyz file.
  atoms = ase.io.read('{0}.xyz'.format(seedname))
  
  # Run the calculation and get the results.
  atoms.set_calculator(calculator)
  
  energy = atoms.get_potential_energy()
  
  if has_forces:
    forces = atoms.get_forces()
  else:
    forces = None
  
  if has_stress:
    stress = atoms.get_stress()
  else:
    stress = None
  
  # Construct the electronic_structrue.dat file.
  write_electronic_structure('electronic_structure.dat',
                             energy,
                             forces=forces,
                             stress=stress)

def read_gulp_file(filename):
  '''
  Read seedname.gin to get settings for gulp.
  N.B. the Gulp file must not have atomic positions or lattice vectors;
     these will be generated by Caesar.
  N.B. the Gulp file must contain a library for the ASE interface to Gulp.
  '''
  
  # Check the file exists, split it into lines,
  #    and split each line into tokens.
  if os.path.isfile(filename):
    lines = [line.rstrip('\n').split() for line in open(filename)]
  else:
    raise ValueError('Settings file {0} does not exist.'.format(filename))
  
  # The Gulp keywords should be the first line.
  keywords = ' '.join(lines[0])
  
  # Check for if forces and/or stresses are being calculated.
  has_forces = False
  has_stress = False
  for keyword in lines[0]:
    key = keyword[:4].lower()
    if key=='grad':
      has_forces = True
    elif key=='stre':
      has_stress = True
  
  # Parse the remainder of the file for other options.
  library = None
  options = []
  for line in lines[1:]:
    if len(line)==0:
      continue
    key = line[0][:4].lower()
    if key=='libr':
      library = ' '.join(line[1:])
    elif key in ['cell','frac','cart']:
      raise ValueError('Settings file {0} contains key {1}. This should not be specified.'.format(filename,line[0]))
    else:
      options.append(' '.join(line))
  
  # Check that 'library' was present in the file.
  if library is None:
    print('Error: settings file {0} does not contain a library.'.format(filename))
    print('A Gulp library is required for the ASE interface to Gulp.')
    print('If a library is not required by the potential,')
    print('   please create a dummy library with blank entries, e.g.')
    print('')
    print('species')
    print('H core 0')
    print('')
    raise ValueError
  
  return keywords, library, options, has_forces, has_stress
  
def write_electronic_structure(filename,energy,forces=None,hessian=None,
   stress=None):
  '''
  Create an electronic_structure.dat file containing energy,
     forces, a hessian and a stress tensor.
  Energy is required, but forces, hessian and stress are optional.
  '''
  
  with open(filename, 'w') as f:
    f.write('Energy (Hartree)\n')
    f.write('{0: .18E}\n'.format(energy/Ha))
    
    if not forces is None:
      f.write('Forces (Hartree/Bohr)\n')
      np.savetxt(f,forces*Bohr/Ha,fmt='% .18E')
    
    if not hessian is None:
      f.write('Hessian (Hartree/Bohr^2)\n')
      np.savetxt(f,hessian*Bohr**2/Ha,fmt='% .18E')
    
    if not stress is None:
      f.write('Stress (Hartree/Bohr^3)\n')
      if stress.shape==(3,3):
        # Stress is in the format [[xx,xy,xz],[xy,yy,yz][xz,yz,zz]].
        np.savetxt(f,stress*Bohr**3/Ha,fmt='% .18E')
      elif stress.shape==(6,):
        # Stress is in the format [xx,yy,zz,yz,xz,xy].
        f.write('{0: .18E} {1: .18E} {2: .18E}\n'.format(stress[0]*Bohr**3/Ha,stress[5]*Bohr**3/Ha,stress[4]*Bohr**3/Ha))
        f.write('{0: .18E} {1: .18E} {2: .18E}\n'.format(stress[5]*Bohr**3/Ha,stress[1]*Bohr**3/Ha,stress[3]*Bohr**3/Ha))
        f.write('{0: .18E} {1: .18E} {2: .18E}\n'.format(stress[4]*Bohr**3/Ha,stress[3]*Bohr**3/Ha,stress[2]*Bohr**3/Ha))
      else:
        raise(ValueError('Stress in unexpected shape.'))

if __name__=='__main__':
  if len(sys.argv)>0 and sys.argv[1] in ['-h', '--help']:
    print('{0} file_type directory no_cores no_nodes seedname'.format(sys.argv[0]))
    sys.exit()
  elif len(sys.argv)!=6:
    raise ValueError('Wrong number of arguments. Call -h or --help for help.')
  
  file_type = sys.argv[1]
  directory = sys.argv[2]
  no_cores = int(sys.argv[3])
  no_nodes = int(sys.argv[4])
  seedname = sys.argv[5]
  
  if file_type!='xyz':
    raise ValueError('This script can only read extended .xyz files.')
  
  main(directory,seedname)
