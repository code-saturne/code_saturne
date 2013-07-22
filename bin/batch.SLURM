#
# Batch options (SLURM)
# =====================
#
#SBATCH --ntasks=2
#SBATCH --time=0:30:00
#SBATCH --partition=par
#SBATCH --exclusive
#SBATCH --output nameandcaseo.%j 
#SBATCH --error nameandcasee.%j 
#SBATCH --job-name=nameandcase
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
#  --qos=<qos>                        : request given quality of service
#                                       (for example "debug")
#  --contiguous                       : use contiguous nodes for minimal latency
#  --exclusive                        : do not share nodes with other jobs
#  --partition=<name>, -p<name>       : partition (queue) (run "sinfo -s"
#                                       for list of partitions)
#  --reservation=<name>               : allocate resources from reservation
#  -o<pattern>, --output=<pattern>    : output file name pattern
#  -e<pattern>, --error=<pattern>     : error file name pattern
#  -J<jobname>, --job-name=<jobname>  : job name
#
# Warning: when coupling with SYRTHES, 1 processor will be reserved for each
#          instance of SYRTHES. The Kernel will be executed on the remaining
#          processors, so make sure to reserve a sufficiently high number
#          of processors.
#
