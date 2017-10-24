import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
import operator
from matplotlib import rc

rc('font', **{'family':'serif','serif':['sffamily']})
rc('text', usetex=True)
params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

def main():
  # Define some static data.
  colours = {
    'turquoise':[102/255,194/255,165/255],
    'Orange'   :[252/255,141/255, 98/255],
    'blue'     :[141/255,160/255,203/255],
    'purple'   :[231/255,138/255,195/255],
    'green'    :[166/255,216/255, 84/255],
    'yellow'   :[255/255,217/255, 47/255],
    'beige'    :[229/255,196/255,148/255],
    'grey'     :[179/255,179/255,179/255]}
  
  fig = plt.figure()
  axes = {}
  
  # Read cutoff convergence file into cutoff structure.
  file_name = 'cutoff_convergence.dat'
  cutoff_file = [line.rstrip('\n').split() for line in open(file_name)]
  cutoff = {'cutoff':[], 'energy difference':[], 'force difference':[]}
  cutoff['final cutoff'] = int(cutoff_file[0][4])
  cutoff['energy tolerance'] = float(cutoff_file[1][3])
  cutoff['force tolerance'] = float(cutoff_file[2][3])
  for line in cutoff_file[5:]:
    cutoff['cutoff'].append(int(line[0]))
    cutoff['energy difference'].append(float(line[1]))
    cutoff['force difference'].append(float(line[2]))
  
  # Read k-point convergence file into kpoint structure
  file_name = 'kpoints_convergence.dat'
  kpoints_file = [line.rstrip('\n').split() for line in open(file_name)]
  kpoints = {'kpoints':[], 'energy difference':[], 'force difference':[]}
  kpoints['final kpoints'] = int(kpoints_file[1][4])
  kpoints['energy tolerance'] = float(kpoints_file[3][3])
  kpoints['force tolerance'] = float(kpoints_file[4][3])
  for line in kpoints_file[7:]:
    kpoints['kpoints'].append(int(line[1]))
    kpoints['energy difference'].append(float(line[5]))
    kpoints['force difference'].append(float(line[6]))
  
  # Plot energy cutoff convergence.
  axes['cutoff'] = fig.add_subplot(211)
  axes['cutoff'].plot(
    cutoff['cutoff'],
    cutoff['energy difference'],
    color=colours['turquoise'])
  axes['cutoff'].set_xlabel(r'Cutoff energy, eV')
  axes['cutoff'].set_ylabel(r'Energy difference, eV/cell')
  axes['cutoff'].set_xlim(cutoff['cutoff'][0], cutoff['cutoff'][-1])
  axes['cutoff'].set_ylim(0, 4*cutoff['energy tolerance'])
  axes['cutoff'].hlines(
    cutoff['energy tolerance'],
    cutoff['cutoff'][0],
    cutoff['cutoff'][-1],
    color=[0,0,0])
  
  # Plot k-point convergence.
  axes['kpoints'] = fig.add_subplot(212)
  axes['kpoints'].plot(
    kpoints['kpoints'],
    kpoints['energy difference'],
    color=colours['turquoise'])
  axes['kpoints'].set_xlabel(r'No. k-points')
  axes['kpoints'].set_ylabel(r'Energy difference, eV/cell')
  axes['kpoints'].set_xlim(kpoints['kpoints'][0], kpoints['kpoints'][-1])
  axes['kpoints'].set_ylim(0, 4*kpoints['energy tolerance'])
  axes['kpoints'].hlines(
    kpoints['energy tolerance'],
    kpoints['kpoints'][0],
    kpoints['kpoints'][-1],
    color=[0,0,0])
  
  # Show result.
  plt.show()

def main_old():
  file_name = "/home/cdt1505/Documents/project.git/MCDFT/run_convergence_kpoints/result_file"
  contents = [line.rstrip('\n').split() for line in open(file_name)]
  
  speciess = ["C", "Ni"]
  variables = ["energy cutoff","unit cell","kpoint spacing"]
  
  energies = {}
  for variable in variables:
    energies[variable] = {}
    for species in speciess:
      energies[variable][species] = {"x":[],"energy":[]}
  
  species = ""
  variable = ""
  reading = False
  for i,line in enumerate(contents):
    if reading == False:
      if len(line) == 1:
        if line[0] == "C":
          species = "C"
        elif line[0] == "Ni":
          species = "Ni"
      elif len(line) >= 2:
        if line[-2:] == ["cluster","energy"]:
          reading = True
          if line[:2] == ["cutoff","energy"]:
            variable = "energy cutoff"
          elif line[:2] == ["unit","cell"]:
            variable = "unit cell"
          elif line[:2] == ["kpoint","spacing"]:
            variable = "kpoint spacing"
    else:
      if len(line) == 0:
        reading = False
      else:
        energies[variable][species]["x"].append(float(line[0]))
        energies[variable][species]["energy"].append(float(line[1]))
  
  fig = plt.figure(figsize=(8,8.5))
  
  axes = {}
  axes["energy cutoff"] = fig.add_subplot(311)
  axes["unit cell"] = fig.add_subplot(312)
  axes["kpoint spacing"] = fig.add_subplot(313)
  
  axes["energy cutoff"].set_xlabel(r"Energy Cutoff, eV")
  axes["energy cutoff"].set_ylabel(r"$\Delta H_f$ Discrepancy, eV")
  axes["unit cell"].set_xlabel(r"Unit cell dimension, $\,\text{\AA}$")
  axes["unit cell"].set_ylabel(r"$\Delta H_f$ Discrepancy, eV")
  axes["kpoint spacing"].set_xlabel(r"k-point spacing, $\text{\AA}^{-1}$")
  axes["kpoint spacing"].set_ylabel(r"$\Delta H_f$ Discrepancy, eV")
  colours = {
    "turquoise":[102/255,194/255,165/255],
    "Ni"       :[252/255,141/255, 98/255],
    "C"        :[141/255,160/255,203/255],
    "NiC"      :[231/255,138/255,195/255],
    "green"    :[166/255,216/255, 84/255],
    "S"        :[255/255,217/255, 47/255],
    "beige"    :[229/255,196/255,148/255],
    "grey"     :[179/255,179/255,179/255]}
  
  energy_final = {}
  for species in speciess:
    # could use "unit cell" or "kpoint spacing" here. the ["energy"][-1] values are all the same.
    energy_final[species] = energies["energy cutoff"][species]["energy"][-1]
  
  previous_energies = {"C":-148.09972410079998895 ,"Ni":-1351.21379329799992774}
  
  print()
  print("single atom energies:")
  for species in speciess:
    print(species," ",str.format('{0:.15e}',energy_final[species]+previous_energies[species]))
  print()
  
  error = 0.312571209055/100
  
  for variable in variables:
    for species in speciess:
      line = energies[variable][species]
      line["discrepancy"] = []
      line["log"] = []
      for energy in line["energy"]:
        discrepancy = abs(energy-energy_final[species])
        line["discrepancy"].append(discrepancy)
        line["log"].append(np.log(discrepancy))
      # y = ax**2+bx+c
      line["fit"] = np.poly1d(np.polyfit(line["x"][:-1],line["log"][:-1],1))
      line["c"] = line["fit"](0)
      line["m"] = line["fit"](1)-line["c"]
      line["x cross"] = (np.log(error)-line["c"])/line["m"]
      print(species," ",variable," x    :",str.format('{0:.15e}',line["x cross"]))
      print(species," ",variable," de/dx:",str.format('{0:.15e}',-line["m"]*error))
      
  
  for variable in variables:
    for species in speciess:
      line = energies[variable][species]
      axes[variable].scatter(
        line["x"],
        line["discrepancy"],
        color=colours[species],
        marker="x",
        label=species,
        s=36)
      if variable != "kpoint spacing":
        axes[variable].plot(
          line["x"],
          np.exp(line["fit"](line["x"])),
          color=colours[species],
          linewidth=2)
  
  ymin = {"energy cutoff":1e-4,"unit cell":1e-4,"kpoint spacing":1e-6}
  ymax = {"energy cutoff":1e3 ,"unit cell":1e3 ,"kpoint spacing":1e0 }
  xmin = {"energy cutoff":100 ,"unit cell":1   ,"kpoint spacing":0.3 }
  xmax = {"energy cutoff":1000,"unit cell":10  ,"kpoint spacing":0.03}
  
  for variable in variables:
    axis = axes[variable]
    axis.set_yscale('log')
    axis.set_xlim(xmin[variable],xmax[variable])
    axis.set_ylim(ymin[variable],ymax[variable])
    axis.legend(loc="upper center",ncol=2)
    axis.hlines(error,0,1000,color=[0,0,0])
    '''for species in speciess:
      if variable != "kpoint spacing":
        axis.vlines(
          energies[variable][species]["x cross"],
          ymin[variable],
          ymax[variable],
          color=colours[species],
          linewidth=2)'''
  
  plt.tight_layout()
  fig.savefig("/home/cdt1505/Documents/project.git/report/convergence.pdf")

main()
