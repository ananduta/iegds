#bin/sh
# Request one node with 4 free processor cores
#PBS -l nodes=1:ppn=4
# Set the name of the job
#PBS -N sim4_b5
# Mail me when the job ends for any reason
#PBS -m ae
# This is my email address
#PBS -M w.ananduta@tudelft.nl


# Go to the directory where I entered the qsub command
cd $PBS_O_WORKDIR


# Activate the Matlab version I want
module load matlab/2020a

# Run my M file and don't even try to display graphics
matlab -nodisplay -r sim4


