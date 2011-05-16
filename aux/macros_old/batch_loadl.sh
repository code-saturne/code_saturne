#
# BATCH OPTIONS (Loadlever)
# =========================
#
#@ shell = /bin/bash
#
#@ job_name = nometcas
#
#@ job_type = mpich
#@ node = 2
#@ tasks_per_node = 6
#@ node_usage = not_shared
#@ class = standard
#
#@ wall_clock_limit = 00:20:00
#
#@ output = $(job_name).$(jobid).out
#@ error  = $(job_name).$(jobid).err
#@ notification = never
#
#@ queue
#
