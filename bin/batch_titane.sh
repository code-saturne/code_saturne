#
# BATCH OPTIONS FOR CEA CCRT's Titane cluster
# ===========================================
#
#MSUB -n 2
#MSUB -T 300
#MSUB -o nometcaso.%J 
#MSUB -e nometcase.%J 
#MSUB -r nometcas
#
#  -n : number of processors
#  -N : number of nodes
#  -T : walltime in seconds
#  -q : queue (run "class" for list of queues)
#  -o : output file name
#  -e : error file name
#  -r : job name
#
