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
cp 1d_corrosion.template workdir.$num/
#cp mesh_file.e workdir.$num/
cp -p dprepro workdir.$num/
cp 1d_corrosion.i workdir.$num/
# cp /home/chaitanya.bhave/projects/moose/modules/phase_field/phase_field-opt workdir.$num/
# RUN the simulation from workdir.num
cd workdir.$num/
#input params.in into the moose simulation
./dprepro --left-delimiter=% --right-delimiter=% params.in 1d_corrosion.template 1d_corrosion.i

# --------
# ANALYSIS
# --------
#submitting the moose job
mpirun -n 10 ~/projects/moose/modules/phase_field/phase_field-opt -i 1d_corrosion.i > job_output.out

#sleep until execution of moose simulation finishes
# jobid=$(tail -1 sbatch.out | egrep -o '[0-9]+')
# while [ $(squeue -j $jobid | wc -l) -ne 0 ];
# do
#   sleep 300
# done

# ---------------
# POST-PROCESSING
# ---------------
# extract function value from the simulation output - Attention with the order!!
#1 line for each response function
#for testing
#tail -n+2 moose_out.csv|cut -f3 -d','

#initial_total_cr
echo '1' >> results.tmp
#head -n+3 2_component_gp/2_component_gp.csv | tail -n 1 | cut -f8 -d',' >> results.tmp
#final total Cr
#tail -n 1 2_component_gp/2_component_gp.csv | cut -f8 -d',' >> results.tmp
#initial metal Cr
#head -n+3 2_component_gp/2_component_gp.csv | tail -n 1 | cut -f2 -d',' >> results.tmp
#final metal Cr
#tail -n 1 2_component_gp/2_component_gp.csv | cut -f2 -d',' >> results.tmp
#inital metal thickness
#head -n+3 2_component_gp/2_component_gp.csv | tail -n 1 | cut -f6 -d',' >> results.tmp
#final metal thickness
#tail -n 1 2_component_gp/2_component_gp.csv | cut -f6 -d',' >> results.tmp

# write results.out and cleanup
cp results.tmp ./results.out
cd ../
