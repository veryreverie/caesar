import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats
import numpy as np
import math
import operator
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
import matplotlib.lines as mlines
import os.path
import sys

rc('font', **{'family':'serif','serif':['sffamily']})
rc('text', usetex=True)
params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

def main():
  colours = {
    'turquoise':[102/255,194/255,165/255],
    'orange'   :[252/255,141/255, 98/255],
    'blue'     :[141/255,160/255,203/255],
    'purple'   :[231/255,138/255,195/255],
    'green'    :[166/255,216/255, 84/255],
    'yellow'   :[255/255,217/255, 47/255],
    'beige'    :[229/255,196/255,148/255],
    'grey'     :[179/255,179/255,179/255]}
  colours_list = list(colours.values())
  
  file_name = 'convergence.dat'
  lines = [line.rstrip('\n').split() for line in open(file_name)]
  xs = []
  ys = []
  threshold = float(lines[0][-1])
  for i,line in enumerate(lines):
    if len(line)==2 and line[0]=='Iteration':
      if len(lines)>=i+4:
        xs.append([float(x) for x in lines[i+2]])
        ys.append([float(y) for y in lines[i+4]])
  
  if len(xs)<2:
    print('Not enough iterations to plot.')
    sys.exit()
  
  dxs = [[] for _ in xs[0]]
  dys = [[] for _ in ys[0]]
  for x,y in zip(xs,ys):
    for j in range(len(x)):
      dxs[j].append(x[j]-xs[-1][j])
      dys[j].append(y[j]-xs[-1][j])
  
  fig, ax_grid = plt.subplots(2,1,gridspec_kw={'height_ratios':[5,1]})
  handles = []
  labels = []
  line = mlines.Line2D([],[],color=[0,0,0], dashes=[10,0])
  handles.append(line)
  labels.append('Input Coefficients')
  line = mlines.Line2D([],[],color=[0,0,0], dashes=[4,1,4,1])
  handles.append(line)
  labels.append('Output Coefficients')
  
  # Colour the least-converged coefficients.
  ds = [abs(dy[-1]-dx[-1]) for dx,dy in zip(dxs,dys)]
  sort_key = sorted(range(len(ds)), key=lambda x: ds[x], reverse=True)
  line_colours = [[0,0,0] for _ in ds]
  for i in range(min(8,len(sort_key))):
    if ds[sort_key[i]]*100>ds[sort_key[0]]:
      line_colours[sort_key[i]] = colours_list[i]
      line = mlines.Line2D([],[],color=colours_list[i])
      handles.append(line)
      labels.append('Coefficient '+str(sort_key[i]+1))
  
  ax = ax_grid[0]
  x_axis = [i+1 for i in range(len(dxs[0]))]
  for c,dx,dy in [(line_colours[i],dxs[i],dys[i]) for i in reversed(sort_key)]:
    ax.plot(x_axis, dx, dashes=[10,0], color=c)
    ax.plot(x_axis, dy, dashes=[4,1,4,1], color=c)
  
  ax.hlines(threshold,0,len(x_axis), linestyle='dotted')
  ax.hlines(-threshold,0,len(x_axis), linestyle='dotted')
  ax.set_xlim(1,len(x_axis))
  ax.set_yscale('symlog', linthreshy=threshold)
  ax.xaxis.set_major_locator(MaxNLocator(integer=True))
  
  ax.set_xlabel('Iteration')
  ax.set_ylabel('Error in Coefficients, (Ha)')
  
  legend = ax_grid[1]
  legend.axis('off')
  legend.legend(handles=handles, labels=labels, loc='center', ncol=5)
  
  plt.show()

main()
