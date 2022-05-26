#SBATCH --nodes=2
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH --partition=par
#SBATCH --output=job_%j.out.log
#SBATCH --error=job_%j.err.log
#SBATCH --job-name=nameandcase
