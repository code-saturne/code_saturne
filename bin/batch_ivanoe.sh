#
# BATCH OPTIONS FOR THE IVANOE CLUSTER (SLURM)
# ============================================
#
#SBATCH --ntasks=2
#SBATCH --time=0:30:00
#SBATCH --partition=para
#SBATCH --output nometcaso.%j 
#SBATCH --error nometcase.%j 
#SBATCH --job-name=nometcas
#
#  -t<time>, --time=<time>            : walltime in minutes, minutes:seconds,
#                                       hours:minutes:seconds, or
#                                       days-hours:minutes:seconds
#  -N, --nodes=<minnodes[-maxnodes]>  : number of allocated nodes
#  --ntasks=24, -n24                  : number of total tasks
#  --ntasks-per-node=<ntasks>         : number of tasks per node
#  --ntasks-per-socket=<ntasks>       : number of tasks per socket
#  --ntasks-per-core=<ntasks>         : number of tasks per core
#  --cpu-bind=cores, sockets          : bind CPUs to cores or sockets
#  --mem-bind=local
#  --contiguous                       : force use of contiguous nodes for
#                                       better latency and minimal latency
#  --partition=<name>, -p<name>       : partition (queue) (run "sinfo -s"
#                                       for list of partitions)
#  --reservation=<name>               : allocate resources from reservation
#  -o<pattern>, --output=<pattern>    : output file name pattern
#  -e<pattern>, --error=<pattern>     : error file name pattern
#  -J<jobname>, --job-name=<jobname>  : job name
#
