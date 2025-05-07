#!/bin/bash
#SBATCH --job-name=ESPIC-two-stream     
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=8               # total number of tasks across all nodes
#SBATCH --cpus-per-task=4        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=jk9075@princeton.edu
#SBATCH --error=job_two_stream.err
#SBATCH --output=job_two_stream.out

srun --ntasks=1 --exclusive --cpus-per-task=4 --mem-per-cpu=4G python3 pic_1d_two_stream.py --interpol "CIC" --method "midpoint"    &
srun --ntasks=1 --exclusive --cpus-per-task=4 --mem-per-cpu=4G python3 pic_1d_two_stream.py --interpol "CIC" --method "implicit"    &
srun --ntasks=1 --exclusive --cpus-per-task=4 --mem-per-cpu=4G python3 pic_1d_two_stream.py --interpol "CIC" --method "verlet"      &
srun --ntasks=1 --exclusive --cpus-per-task=4 --mem-per-cpu=4G python3 pic_1d_two_stream.py --interpol "CIC" --method "leapfrog"    &
srun --ntasks=1 --exclusive --cpus-per-task=4 --mem-per-cpu=4G python3 pic_1d_two_stream.py --interpol "TSC" --method "midpoint"    &
srun --ntasks=1 --exclusive --cpus-per-task=4 --mem-per-cpu=4G python3 pic_1d_two_stream.py --interpol "TSC" --method "implicit"    &
srun --ntasks=1 --exclusive --cpus-per-task=4 --mem-per-cpu=4G python3 pic_1d_two_stream.py --interpol "TSC" --method "verlet"      &
srun --ntasks=1 --exclusive --cpus-per-task=4 --mem-per-cpu=4G python3 pic_1d_two_stream.py --interpol "TSC" --method "leapfrog"    &
wait