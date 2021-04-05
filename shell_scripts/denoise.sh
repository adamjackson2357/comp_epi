#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=24:mem=40gb
#PBS -N denoise

cd /rds/general/user/aj1520/home/Group8/Adam/comp_epi/code/

module load anaconda3/personal
source activate comp_epi_env

Rscript fast_denoise.R
