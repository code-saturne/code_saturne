#
# BATCH OPTIONS (Sun Grid Engine)
# ===============================
#
# set the name of the job
##$ -N nometcas
#
# request between 2 and 4 slots
##$ -pe mpich 2-4
#
# Execute the job from the current working directory
# Job output will appear in this directory
##$ -cwd
#   can use -o dirname to redirect stdout
#   can use -e dirname to redirect stderr

#  Export these environment variables
##$ -v MPI_HOME
#
