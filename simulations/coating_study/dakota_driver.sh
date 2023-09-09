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
# cp 1d_corrosion.template workdir.$num/
cp gen_microstructure_with_salt.py workdir.$num/gen_microstructure_with_salt.template
# cp -p dprepro workdir.$num/
cp ebsd_reader.i workdir.$num/
cp ni20cr_corr.i workdir.$num/ni20cr_corr.template

cp moose_job_script.sh workdir.$num/
# cp /home/chaitanya.bhave/projects/moose/modules/phase_field/phase_field-opt workdir.$num/
# RUN the simulation from workdir.num
cd workdir.$num/
#input params.in into the moose simulation
# ---------------
# PRE-PROCESSING
# ---------------
dprepro --left-delimiter=% --right-delimiter=% params.in gen_microstructure_with_salt.template gen_microstructure_with_salt.py
dprepro --left-delimiter=% --right-delimiter=% params.in ni20cr_corr.template ni20cr_corr.i

python gen_microstructure_with_salt.py

# #submitting the moose job
qsub -v job_site=`pwd` moose_job_script.sh > submit_moose.out
# #sleep until execution of moose simulation finishes
#
jobid=$(tail -1 submit_moose.out | egrep -o '[0-9]+')
status=`qstat | grep $(whoami)  | grep $jobid` # check to see if job is running
while [ -n "$status" ] # while $status is not empty
	do
		sleep 300
		status=`qstat | grep $(whoami)  | grep $jobid`
	done
# ---------------
# POST-PROCESSING
# ---------------
# extract function value from the simulation output - Attention with the order!!
#1 line for each response function
#for testing
#tail -n+2 moose_out.csv|cut -f3 -d','

#initial_total_cr
echo '1' > results.tmp
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
