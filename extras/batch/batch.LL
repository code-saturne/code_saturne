#-------------------------------------------------------------------------------
#
# Batch options for IBM LoadLeveler (example: IBM cluster under AIX)
# =================================
#
#@ shell = /bin/sh
#
#@ job_name = nameandcase
#
#@ job_type = parallel
#@ cpus = 128
#@ node_usage = not_shared
#
#@ network.MPI = csss,shared,US
#@ bulkxfer = yes
#
#@ wall_clock_limit = 00:20:00
#@ account_no = z001
#
#@ output = $(job_name).$(schedd_host).$(jobid).out
#@ error  = $(job_name).$(schedd_host).$(jobid).err
#@ notification = never
#
#@ queue
# suggested environment settings:
#  export MP_EAGER_LIMIT=65536
#  export MP_SHARED_MEMORY=yes
#  export MEMORY_AFFINITY=MCM
#  export MP_TASK_AFFINITY=MCM
#
#-------------------------------------------------------------------------------

# Change to submission drectory
cd $LOADL_STEP_INITDIR

# Launch runcase script
./scriptname
