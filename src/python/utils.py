# ======================================================================
# Shared functions between python plotting scripts.
# ======================================================================
import os.path

def read_file(filename):
  if os.path.isfile(filename):
    return [line.rstrip('\n').split() for line in open(filename)]
  else:
    raise ValueError('File '+filename+' does not exist.')

def split_into_sections(flines):
  sections = [[]]
  for line in flines:
    if len(line)==0:
      if sections[-1]:
        sections.append([])
    else:
      sections[-1].append(line)
  if not sections[-1]:
    sections = sections[:-1]
  return sections
  
def dblecomplex(displacement):
  if displacement[0]=='-':
    re = float(displacement[:25])
    im = float(displacement[25:-1])
  else:
    re = float(displacement[:24])
    im = float(displacement[24:-1])
  return complex(re,im)

def parse_mode_maps(filename):
  lines = read_file(filename)
  
  frequencies = [float(x) for x in lines[0][2:]]
  
  # Split file into modes.
  modes = []
  for line in lines[2:]:
    if len(line)>0:
      if line[0]=='Mode':
        modes.append({'ID':int(line[2])})
      elif line[0]=='Harmonic':
        modes[-1]['Harmonic frequency'] = float(line[2])
      elif line[0]=='Displacement':
        sampling = 'Sampled' in line
        using_stress = 'pressure' in line
        
        modes[-1]['Displacements'] = []
        modes[-1]['Harmonic energies'] = []
        modes[-1]['Harmonic forces'] = []
        modes[-1]['Anharmonic energies'] = []
        modes[-1]['Anharmonic forces'] = []
        if using_stress:
          modes[-1]['Anharmonic pressures'] = []
        if sampling:
          modes[-1]['Sampled energies'] = []
          modes[-1]['Sampled forces'] = []
          if using_stress:
            modes[-1]['Sampled pressures'] = []
      else:
        modes[-1]['Displacements'].append(float(line[0]))
        modes[-1]['Harmonic energies'].append(float(line[1]))
        modes[-1]['Harmonic forces'].append(float(line[2]))
        modes[-1]['Anharmonic energies'].append(float(line[3]))
        modes[-1]['Anharmonic forces'].append(float(line[4]))
        if using_stress:
          modes[-1]['Anharmonic pressures'].append(float(line[5]))
          if sampling:
            modes[-1]['Sampled energies'].append(float(line[6]))
            modes[-1]['Sampled forces'].append(float(line[7]))
            modes[-1]['Sampled pressures'].append(float(line[8]))
        else:
          if sampling:
            modes[-1]['Sampled energies'].append(float(line[5]))
            modes[-1]['Sampled forces'].append(float(line[6]))
  
  if sampling:
    for mode in modes:
      mode['Harmonic energy difference'] = []
      for harmonic,sampled in zip(mode['Harmonic energies'],
                                  mode['Sampled energies']):
        mode['Harmonic energy difference'].append(harmonic-sampled)
      
      mode['Harmonic force difference'] = []
      for harmonic,sampled in zip(mode['Harmonic forces'],
                                  mode['Sampled forces']):
        mode['Harmonic force difference'].append(harmonic-sampled)
      
      mode['Anharmonic energy difference'] = []
      for anharmonic,sampled in zip(mode['Anharmonic energies'],
                                    mode['Sampled energies']):
        mode['Anharmonic energy difference'].append(anharmonic-sampled)
      
      mode['Anharmonic force difference'] = []
      for anharmonic,sampled in zip(mode['Anharmonic forces'],
                                    mode['Sampled forces']):
        mode['Anharmonic force difference'].append(anharmonic-sampled)
      if using_stress:
        mode['Anharmonic pressure difference'] = []
        for anharmonic,sampled in zip(mode['Anharmonic pressures'],
                                      mode['Sampled pressures']):
          mode['Anharmonic pressure difference'].append(anharmonic-sampled)
  
  return frequencies, sampling, using_stress, modes

