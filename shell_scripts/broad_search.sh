#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=1:mem=20gb
#PBS -N broad_search

module load anaconda3/personal
source activate comp_epi_env

step=10
sample=1000

cd /rds/general/user/aj1520/home/Group8/Adam/comp_epi/code/

Rscript variables_broad_search.R $step $sample
