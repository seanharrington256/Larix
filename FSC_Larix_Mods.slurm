#!/bin/bash

#SBATCH --job-name FSC_Larix
#SBATCH -A inbreh
#SBATCH -t 7-00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=seanharrington256@gmail.com
#SBATCH -e /project/inbreh/larix_turck/scripts_larix/err_out/err_FSC_Larix_%A_%a.err
#SBATCH -o /project/inbreh/larix_turck/scripts_larix/err_out/std_FSC_Larix_%A_%a.out
#SBATCH --array=1-400

# Number of jobs for array is 50 * number of models



cd /project/inbreh/larix_turck/Reps   # set the correct working directory - make sure that the array number is correct, too

for x in */*; do
  repdirs=(${repdirs[@]} "${x}")
done

cd ${repdirs[($SLURM_ARRAY_TASK_ID-1)]}   ## here, I subtract 1 from the $SLURM_ARRAY_TASK_ID because bash indexing starts at zero - I think it's less confusing to subtract 1 here than to remember to do it when specifying the number of jobs for the array

est=$(echo *.est)
tpl=$(echo *.tpl)

echo "starting FSC on slurm job $SLURM_ARRAY_TASK_ID, $tpl"

fsc26 -t $tpl -e $est -n 100000 -L 40 -m -M -u -x  -q

echo "completed FSC on slurm job $SLURM_ARRAY_TASK_ID, $tpl"
