#!/usr/bin/python
# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys,os
if(len(sys.argv)==1):
    print "Generate pbs script for an input file"
    print "Usage: pbsgen.py input.i"
    sys.exit(0)
outfilename=sys.argv[1].split('.')[0]+".pbs"

# <codecell>

template="""#!/bin/sh 
#  This is a pbs script for fusion at Argonne National Lab
#  See info at:  http://www.lcrc.anl.gov/info/BatchJobs 
#  
#  Fusion has two primary queues: shared and batch:
#  shared queue is limited to jobs less than 4 nodes and 1 hours 
#  please use #PBS -q shared or qsub -q shared for shared queue 
#  Submit this script using the command: qsub script.sample                                                                                                                                                           
#
#  Use the "qstat -a" command to check the status of a job.
#                                                                                                                                                                                                                     
#  To delete a script from the queue:
#       qdel <job_number>
#  or
#       qdel -W force <job_number>
#
# The following are embedded QSUB options. The syntax is #PBS (the # does
# _not_  denote that the lines are commented out so do not remove).
#
# resource limits  walltime: maximum wall clock time (hh:mm:ss)
#                  nodes: number of 2-processor nodes
#                  ppn: how many processors per node to use (1 or 2)
#                      (you are always charged for the entire node)
#                  pt: production Platinum cluster nodes
#PBS -l nodes=4:ppn=1
#PBS -l walltime=0:59:59
#                                                                                                                                                                                                                     
# export all my environment variables to the job
#PBS -V

# combine standard output and standard error (optional)
#PBS -j oe
#
# send mail when the job begins and ends (optional)
# PBS -m be
# End of embedded QSUB options
#

# set $NN to have the current number of nodes
export NN=`wc -l ${PBS_NODEFILE} | awk '{print $1}'`
echo '========================================'
# print out the list of nodes
echo 'NODES: '
cat $PBS_NODEFILE
echo '========================================'
echo ${NN}

# go to submission
cd $PBS_O_WORKDIR
"""

# <codecell>

F=os.getenv("F")
if F is None:
    F=raw_input("""executable name(/home/sgu/softwares-blues/herd-master/ferret-opt):""") or """/home/sgu/softwares-blues/herd-master/ferret-opt"""
with open(outfilename,'w') as f:
    f.write(template)
    f.write("export F="+F+"\n")
    f.write("mpiexec $F -i "+sys.argv[1]+"\n")

# <codecell>


