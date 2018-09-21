#!/bin/bash
#!
#! Example Slurm submission script.
#! Runs caesar as a single serial instance, but with enough resources to
#!    run electronic structure calculations in parallel.

CMD="caesar run_harmonic -f run_harmonic.input_file --no_cores $SLURM_NTASKS --no_nodes $SLURM_JOB_NUM_NODES"

cd $SLURM_SUBMIT_DIR
eval $CMD
