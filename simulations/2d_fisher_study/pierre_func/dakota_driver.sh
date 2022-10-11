#! /bin/bash
# --------------
# PRE-PROCESSING
# --------------
# Incorporate the parameters from DAKOTA into the template
#/blue/michael.tonks/chaitanya.bhave/yellowjacket/first_paper/sensitivity_analysis

# getting workdir tag
name=${PWD##*/}

export num=`echo $name | cut -c 9-`

# copy any necessary files and data files into workdir
cd ../
cp 2_grain.template workdir.$num/
cp -p dprepro workdir.$num/
cp 2_grain.i workdir.$num/
cd workdir.$num/
#input params.in into the moose simulation
./dprepro --left-delimiter=% --right-delimiter=% params.in 2_grain.template 2_grain.i

# --------
# ANALYSIS
# --------
#submitting the moose job
mpirun -n 120 ~/projects/moose/modules/phase_field/phase_field-opt -i 2_grain.i > job_output.out


echo '1' > results.tmp

# write results.out and cleanup
cp results.tmp ./results.out
cd ../
