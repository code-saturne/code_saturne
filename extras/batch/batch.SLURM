#-------------------------------------------------------------------------------
#
# Batch options for SLURM (Simple Linux Utility for Resource Management) 
# =======================
#
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH --partition=par
#SBATCH --output=job_%j.out.log
#SBATCH --error=job_%j.err.log
#SBATCH --job-name=nameandcase
#
#  -t<time>, --time=<time>            : walltime in minutes, minutes:seconds,
#                                       hours:minutes:seconds, days-hours,
#                                       days-hours:minutes, or
#                                       days-hours:minutes:seconds
#  -N, --nodes=<minnodes[-maxnodes]>  : number of allocated nodes
#  --ntasks=24, -n24                  : number of total tasks
#  --ntasks-per-node=<ntasks>         : number of tasks per node
#  --ntasks-per-socket=<ntasks>       : number of tasks per socket
#  --ntasks-per-core=<ntasks>         : number of tasks per core
#  --cpus-per-task=<nthreads>         : number of threads per task
#  --cpu-bind=cores, sockets          : bind CPUs to cores or sockets
#  --mem-bind=local
#  --qos=<qos>                        : request given quality of service
#                                       (for example "debug")
#  --contiguous                       : use contiguous nodes for minimal latency
#  --exclusive                        : do not share nodes with other jobs
#  --switches=<count>[@max-delay]     : try to run on no more than count
#                                       switches (for better performance)
#  --partition=<name>, -p<name>       : partition (queue) (run "sinfo -s"
#                                       for list of partitions)
#  --reservation=<name>               : allocate resources from reservation
#  -o<pattern>, --output=<pattern>    : output file name pattern
#  -e<pattern>, --error=<pattern>     : error file name pattern
#  -J<jobname>, --job-name=<jobname>  : job name
#
#-------------------------------------------------------------------------------

# Change to submission directory
if test -n "$SLURM_SUBMIT_DIR" ; then cd $SLURM_SUBMIT_DIR ; fi

