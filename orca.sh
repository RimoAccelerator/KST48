# A sample bash script to run ORCA tasks on an HPC cluster, which adapts with KST48. By using this script, the orca_comm can be set as: orca_comm = ./orca.sh 
export job=${1/\.inp/}

export PATH=/home/soft/openmpi411/bin:$PATH
export LD_LIBRARY_PATH=/home/soft/openmpi411/lib:$LD_LIBRARY_PATH
export orcadir=/home/soft/orca5
export RSH_COMMAND="/usr/bin/ssh -x"
export PATH=$orcadir:$PATH

if [ ! -d /tmp/$USER ]
then
  mkdir -p /tmp/$USER
fi

tdir=$(mktemp -d /tmp/$USER/orcajob_XXXX)
export WORKDIR=`pwd`
cd $tdir
mkdir JOBS
cp $WORKDIR/JOBS/* JOBS
$orcadir/orca $job.inp >> $WORKDIR/$job.out --allow-run-as-root

cp $tdir/$job.* $WORKDIR/JOBS
