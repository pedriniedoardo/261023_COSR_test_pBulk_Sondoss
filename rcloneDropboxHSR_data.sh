#!/bin/bash
#SBATCH --job-name=rclone
#SBATCH --account pedrini.edoardo
#SBATCH --mem=8GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=4  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="rcloneDropboxHSR_data.err"
#SBATCH --output="rcloneDropboxHSR_data.out"

echo "my job strart now" > rcloneDropboxHSR_data.log;

date >> rcloneDropboxHSR_data.log;

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_rclone;

rclone copy -PL /beegfs/scratch/ric.cosr/ric.bonanomi/edoardo_pedrini/misc/261023_test_pBulk_Sondoss/data \
dropbox_HSR:project_COSR/ric.bonanomi/edoardo_pedrini/misc/261023_test_pBulk_Sondoss/data

date >> rcloneDropboxHSR_data.log;