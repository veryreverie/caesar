gfortran band_folding.f90 -o band_folding
gfortran calculate_anharmonic.f90 -o calculate_anharmonic
gfortran calculate_bs.f90 -o calculate_bs
gfortran calculate_gap.f90 -o calculate_gap
gfortran generate_amplitudes.f90 -o generate_amplitudes
gfortran generate_quadratic_configurations.f90 -o generate_quadratic_configurations
gfortran generate_sc_path.f90 -o generate_sc_path
gfortran quadratic_spline.f90 -o quadratic_spline
gfortran -lblas -llapack vscf_1d.f90 -o vscf_1d
