project: caesar
summary: A utility for calculating the vibrational free energy of periodic crystals.
author: Mark Johnson
email: mj464@cam.ac.uk
project_github: https://github.com/veryreverie/caesar
src_dir: ../src
exclude_dir: ../src/old/old/harmonic/programs
             ../src/old/old/quadratic/programs
             ../src/old/quadratic
             ../src/old/coupled
             ../src/anharmonic/old
exclude: spglib_dummy.f90
         spglib_dummy_submodule.f90
include: ../src/utils/macros
output_dir: ./ford
graph_dir: ./ford/graphs
docmark: !
predocmark: >
display: public
         protected
         private
source: true
graph: true
coloured_edges: true
graph_maxdepth: 2
search: true
warn: false
parallel: 0
