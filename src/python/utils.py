# ======================================================================
# Shared functions between python plotting scripts.
# ======================================================================
import os.path
import math
import numpy as np
import matplotlib.ticker as mticker

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
  
  if len(lines[0])>0:
    frequencies = [float(x) for x in lines[0][2:]]
    lines = lines[2:]
  else:
    frequencies = []
    lines = lines[1:]
  
  # Split file into modes.
  modes = []
  for line in lines:
    if len(line)>0:
      if line[0]=='Mode' and line[1]=='ID:':
        modes.append({'ID':int(line[2])})
      elif line[0]=='Harmonic':
        modes[-1]['Harmonic frequency'] = float(line[2])
      elif line[0]=='Mode' and line[1]=='displacement':
        anharmonic = 'Anharmonic' in line
        sampled = 'Sampled' in line
        pressure = 'pressure' in line
        
        modes[-1]['Mode displacements'] = []
        modes[-1]['L2 Cartesian displacements'] = []
        modes[-1]['Harmonic energies'] = []
        modes[-1]['Harmonic forces'] = []
        if anharmonic:
          modes[-1]['Anharmonic energies'] = []
          modes[-1]['Anharmonic forces'] = []
          if pressure:
            modes[-1]['Anharmonic pressures'] = []
        if sampled:
          modes[-1]['Sampled energies'] = []
          modes[-1]['Sampled forces'] = []
          if pressure:
            modes[-1]['Sampled pressures'] = []
      else:
        modes[-1]['Mode displacements'].append(float(line[0]))
        modes[-1]['L2 Cartesian displacements'].append(float(line[1]))
        modes[-1]['Harmonic energies'].append(float(line[2]))
        modes[-1]['Harmonic forces'].append(float(line[3]))
        i=3
        if anharmonic:
          modes[-1]['Anharmonic energies'].append(float(line[i+1]))
          modes[-1]['Anharmonic forces'].append(float(line[i+2]))
          i += 2
          if pressure:
            modes[-1]['Anharmonic pressures'].append(float(line[i+1]))
            i += 1
        if sampled:
          modes[-1]['Sampled energies'].append(float(line[i+1]))
          modes[-1]['Sampled forces'].append(float(line[i+2]))
          i += 2
          if pressure:
            modes[-1]['Sampled pressures'].append(float(line[i+1]))
  
  for mode in modes:
    if sampled:
      mode['Harmonic energy difference'] = []
      for harmonic,sampled in zip(mode['Harmonic energies'],
                                  mode['Sampled energies']):
        mode['Harmonic energy difference'].append(harmonic-sampled)
      
      mode['Harmonic force difference'] = []
      for harmonic,sampled in zip(mode['Harmonic forces'],
                                  mode['Sampled forces']):
        mode['Harmonic force difference'].append(harmonic-sampled)
      
      if anharmonic:
        mode['Anharmonic energy difference'] = []
        for anharmonic,sampled in zip(mode['Anharmonic energies'],
                                      mode['Sampled energies']):
          mode['Anharmonic energy difference'].append(anharmonic-sampled)
        
        mode['Anharmonic force difference'] = []
        for anharmonic,sampled in zip(mode['Anharmonic forces'],
                                      mode['Sampled forces']):
          mode['Anharmonic force difference'].append(anharmonic-sampled)
        
        if pressure:
          mode['Anharmonic pressure difference'] = []
          for anharmonic,sampled in zip(mode['Anharmonic pressures'],
                                        mode['Sampled pressures']):
            mode['Anharmonic pressure difference'].append(anharmonic-sampled)
  
  return frequencies, anharmonic, sampled, pressure, modes

def set_xaxis(axis,limits,label,units):
  '''
  Set limits, scientific notation and offsets.
  '''
  axis.set_xlim(limits)
  offset = math.floor(np.log10(max(abs(limits[0]),abs(limits[1]))))
  axis.ticklabel_format(axis='x', style='sci',scilimits=(offset,offset))
  if offset==0:
    axis.set_xlabel(label+' ('+units+')')
  else:
    axis.set_xlabel(label+' ('+r"$\times10^{{{}}}$".format(offset)+' '+units+')')
  sf = mticker.ScalarFormatter()
  sf.set_useOffset(False)
  sf.set_scientific(True)
  sf.set_powerlimits((-1,1))
  axis.xaxis.set_major_formatter(sf)
  axis.xaxis.offsetText.set_visible(False)

def set_yaxis(axis,limits,label,units):
  '''
  Set limits, scientific notation and offsets.
  '''
  axis.set_ylim(limits)
  offset = math.floor(np.log10(max(abs(limits[0]),abs(limits[1]))))
  axis.ticklabel_format(axis='y', style='sci',scilimits=(offset,offset))
  if offset==0:
    axis.set_ylabel(label+' ('+units+')')
  else:
    axis.set_ylabel(label+' ('+r"$\times10^{{{}}}$".format(offset)+' '+units+')')
  sf = mticker.ScalarFormatter()
  sf.set_useOffset(False)
  sf.set_scientific(True)
  sf.set_powerlimits((-1,1))
  axis.yaxis.set_major_formatter(sf)
  axis.yaxis.offsetText.set_visible(False)

