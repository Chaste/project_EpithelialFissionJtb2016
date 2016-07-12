#!/bin/bash
#
# Script to illustrate the running of batch jobs for
# the Almet et al. (2016) crypt fission model.
#
# This script assumes the following command has been run
# in the Chaste folder:
# scons b=GccOpt co=1 ts=projects/EpithelialFissionJtb2016/test/TestCryptFissionSweepsLiteratePaper.hpp

#Specify the first and last seeds over which the simulation is run.
initial_seed=$1;
final_seed=$2;


for ((i=initial_seed; i<=final_seed; i++))
do

echo "Beginning run ${i}."

# "nice" lets other more computationally intensive jobs move up the queue (provided there
# are at least 20 processes running on the machine).
# ">" directs std::cout to the file.
# "2>&1" directs std::cerr to the same place.
# "&" on the end lets the script carry on and not wait until this has finished
# (provided the machine on which this script is executed has multiple cores---the original
# simulations for this model were run on multi-core CPU servers).

nice ../build/optimised/TestCryptFissionSweepsLiteratePaperRunner -myintval ${i} > epithelial_layer_sweeps_${i}_output.txt 2>&1 &

done

echo "Jobs submitted"