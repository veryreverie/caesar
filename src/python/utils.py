# ======================================================================
# Shared functions between python plotting scripts.
# ======================================================================
def dblecomplex(displacement):
  if displacement[0]=='-':
    re = float(displacement[:24])
    im = float(displacement[24:-1])
  else:
    re = float(displacement[:23])
    im = float(displacement[23:-1])
  return complex(re,im)

def parse_mode_maps(filename):
  lines = [line.rstrip('\n').split() for line in open(filename)]
  
  sampling = False
  
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
        modes[-1]['Displacements'] = []
        modes[-1]['Harmonic energies'] = []
        modes[-1]['Harmonic forces'] = []
        modes[-1]['Anharmonic energies'] = []
        modes[-1]['Anharmonic forces'] = []
        if line[-2]=='Sampled':
          sampling = True
          modes[-1]['Sampled energies'] = []
          modes[-1]['Sampled forces'] = []
      else:
        modes[-1]['Displacements'].append(float(line[0]))
        modes[-1]['Harmonic energies'].append(float(line[1]))
        modes[-1]['Harmonic forces'].append(float(line[2]))
        modes[-1]['Anharmonic energies'].append(float(line[3]))
        modes[-1]['Anharmonic forces'].append(float(line[4]))
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
  
  return frequencies, sampling, modes

