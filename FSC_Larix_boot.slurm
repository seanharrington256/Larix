#!/bin/bash

#SBATCH --job-name FSC_Larix_boot
#SBATCH -A inbreh
#SBATCH -t 5-00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/inbreh/larix_turck/scripts_larix/err_out/err_FSC_Larix_boot_%A_%a.err
#SBATCH -o /project/inbreh/larix_turck/scripts_larix/err_out/std_FSC_Larix_boot_%A_%a.out
#SBATCH --array=1-5000



# set the correct working directory - make sure that the array number is correct, too
cd /project/inbreh/larix_turck/boot_inputs/Reps   

for x in */*; do
  repdirs=(${repdirs[@]} "${x}")
done

cd ${repdirs[($SLURM_ARRAY_TASK_ID-1)]}   ## here, I subtract 1 from the $SLURM_ARRAY_TASK_ID because bash indexing starts at zero - I think it's less confusing to subtract 1 here than to remember to do it when specifying the number of jobs for the array

est=$(echo *.est)
tpl=$(echo *.tpl)

fsc26 -t $tpl -e $est -n 100000 -L 40 -m -M -u -x -q
