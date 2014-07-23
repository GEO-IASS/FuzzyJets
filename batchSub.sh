#!/bin/bash

[ "$USER" == "chstan" ]   && WorkDir=/u/at/chstan/nfs/summer_2014/ForConrad/

SubFileLoc=`pwd`/_batchSingleSub.sh

DateSuffix=`date +%Y%m%d_%Hh%Mmin`

echo '#!/bin/bash
echo CD to $1
echo CMD is $2

cd $1
source setup.sh
cmd=$4

echo MAKING TEMP DIR $2
JOBFILEDIR=$2
mkdir $JOBFILEDIR
REALOUT=$3
echo MADE TEMP DIR $JOBFILEDIR
echo WILL COPY TO $REALOUT

shift
shift
echo Calling $cmd $*
$cmd $*
cp -r $JOBFILEDIR/*.root $REALOUT
echo COPYING to $REALOUT
rm -rf $JOBFILEDIR
' > $SubFileLoc
chmod u+x $SubFileLoc

#----------------
Process=2
pThatMin=200
pThatMax=400
BosonMass=800

for mu in 20; do
    Queue=long
    nevents=10000
    njobs=10
    [ "$mu" == "40" ]  && Queue=xlong  && njobs=200  && nevents=500
    [ "$mu" == "60" ]  && Queue=xlong  && njobs=200  && nevents=500
    [ "$mu" == "80" ]  && Queue=xlong  && njobs=200  && nevents=500
    LogPrefix=`pwd`/logs/${DateSuffix}/${DateSuffix}_bsub_${mu}_
    OutDirFinal=`pwd`/files/${DateSuffix}
    mkdir -p `dirname $LogPrefix`
    mkdir -p $OutDirFinal
    echo
    echo "Submitting $njobs jobs each with $nevents events to $Queue"
    echo $LogPrefix
    for (( ii=1; ii<=$njobs; ii++ )) ;  do
        echo $ii
        OutDir=/scratch/${DateSuffix}_${ii}/
        bsub -q ${Queue} -R 'select[(!preempt&&rhel60&&cvmfs&&inet)]' -o $LogPrefix${ii}.log $SubFileLoc           \
             ${WorkDir} ${OutDir} ${OutDirFinal} ./FCNC.exe  \
             --OutFile ${OutDir}/Sample_mu_${mu}_nevents_${nevents}_job_${ii}_Process_${Process}_${pThatMin}_${pThatMax}.root \
             --NEvents ${nevents} \

             done
    done
