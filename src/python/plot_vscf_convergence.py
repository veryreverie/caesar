import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats
import numpy as np
import math
import operator
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import os.path
import sys

rc('font', **{'family':'serif','serif':['sffamily']})
rc('text', usetex=True)
params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

def main():
  colour_names = ['turquoise',
                  'orange'   ,
                  'blue'     ,
                  'purple'   ,
                  'green'    ,
                  'yellow'   ,
                  'beige'    ,
                  'grey'     ]
  colour_values = [[102/255,194/255,165/255],
                   [252/255,141/255, 98/255],
                   [141/255,160/255,203/255],
                   [231/255,138/255,195/255],
                   [166/255,216/255, 84/255],
                   [255/255,217/255, 47/255],
                   [229/255,196/255,148/255],
                   [179/255,179/255,179/255]]
  colours = {name:value for name,value in zip(colour_names,colour_values)}
  
  file_name = 'convergence.dat'
  lines = [line.rstrip('\n').split() for line in open(file_name)]
  xs = []
  ys = []
  fs = []
  threshold = float(lines[0][-1])
  for i,line in enumerate(lines):
    if len(line)==2 and line[0]=='Iteration':
      if len(lines)>=i+6:
        xs.append([float(x) for x in lines[i+2]])
        ys.append([float(y) for y in lines[i+4]])
        fs.append(float(lines[i+6][0]))
  
  if len(xs)<2:
    print('Not enough iterations to plot.')
    sys.exit()
  
  min_i = np.argmin(fs)
  #min_i = np.argmin([np.linalg.norm(np.asarray(y)-np.asarray(x)) for x,y in zip(xs,ys)])
  dxs = [[] for _ in xs[0]]
  dys = [[] for _ in ys[0]]
  for x,y in zip(xs,ys):
    for j in range(len(x)):
      #dxs[j].append(x[j]-(xs[min_i][j]+ys[min_i][j])/2)
      #dys[j].append(y[j]-(xs[min_i][j]+ys[min_i][j])/2)
      dxs[j].append(x[j]-ys[min_i][j])
      dys[j].append(y[j]-ys[min_i][j])
  dfs = [f-fs[min_i] for f in fs]
  
  handles = []
  labels = []
  line = mlines.Line2D([],[],marker='s',color='none',markerfacecolor=[0,0,0])
  handles.append(line)
  labels.append('Input Coefficients')
  line = mlines.Line2D([],[],marker='D',color='none',markerfacecolor=[0,0,0])
  
  handles.append(line)
  labels.append('Output Coefficients')
  
  # Colour the least-converged coefficients.
  ds = [abs(dy[-1]-dx[-1]) for dx,dy in zip(dxs,dys)]
  sort_key = sorted(range(len(ds)), key=lambda x: ds[x], reverse=True)
  line_colours = [[0,0,0,0.2] for _ in ds]
  highlighted = []
  for i in range(min(8,len(sort_key))):
    if ds[sort_key[i]]*100>ds[sort_key[0]]:
      highlighted.append(sort_key[i])
      line_colours[sort_key[i]] = colour_values[i]
      line = mlines.Line2D([],[],color=colour_values[i])
      handles.append(line)
      labels.append('Coefficient '+str(sort_key[i]+1))
  
  if len(highlighted)>2:
    highlighted = highlighted[:2]
  fig, ax_grid = plt.subplots(2,
                              2+len(highlighted),
                              gridspec_kw={'height_ratios':[5,1]},
                              squeeze=False)
  
  x_axis = [i+1 for i in range(len(dxs[0]))]
  
  ax_grid[0][1].scatter(x_axis, dfs, s=3)
  
  ax = ax_grid[0][0]
  
  for c,dx,dy in [(line_colours[i],dxs[i],dys[i]) for i in reversed(sort_key)]:
    ax.scatter(x_axis, dx, s=3, marker='s', color=[c])
    ax.scatter(x_axis, dy, s=3, marker='D', color=[c])
  
  ax.hlines(threshold,0,len(x_axis)+1, linestyle='dotted')
  ax.hlines(-threshold,0,len(x_axis)+1, linestyle='dotted')
  ax.set_xlim(0,len(x_axis)+1)
  ax.set_yscale('symlog', linthreshy=threshold)
  ax.xaxis.set_major_locator(MaxNLocator(integer=True))
  
  ax.set_xlabel('Iteration')
  ax.set_ylabel('Error in Coefficients, (Ha)')
  
  for ax,i in zip(ax_grid[0][2:],highlighted):
    colour_map = []
    for j in range(len(dxs[i])):
      colour_map.append([line_colours[i][k]*(j+1)/len(dxs[i]) for k in range(3)])
    ax.scatter(dxs[i],dys[i],color=colour_map,s=3)
    axlim = 2*max([max(abs(x),abs(y)) for x,y in zip(dxs[i],dys[i])])
    ax.plot([-axlim,axlim],[-axlim,axlim],color=[0,0,0],linestyle='dotted')
    ax.set_xlim(-axlim,axlim)
    ax.set_ylim(-axlim,axlim)
    #ax.set_xscale('symlog', linthreshx=threshold)
    #ax.set_yscale('symlog', linthreshy=threshold)
  
  legend = ax_grid[1][0]
  legend.axis('off')
  legend.legend(handles=handles, labels=labels, loc='center', ncol=5)
  
  for ax in ax_grid[1][1:]:
    ax.axis('off')
  
  plt.show()

main()
