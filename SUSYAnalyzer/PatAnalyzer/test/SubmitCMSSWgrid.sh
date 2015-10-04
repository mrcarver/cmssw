#!/bin/bash

mkdir -p HPCLogs/
mkdir -p HPCErrors/
cfg=SignalTriggers.py
count=0
input=$1

function usage() {
echo "sh SubmitCMSSW.sh <input_directory> <output_directory> <files_per_job> <events_per_job>"
exit;
}

if [ ! -n "$1" ] || [ ! -n "$2" ]  || [ ! -n "$3" ] ; then usage; fi;


INPUT=$1
if [ ! -d $INPUT ] ; then usage ; fi;
OUTPUT=$2
mkdir -p $OUTPUT
mkdir -p ${OUTPUT}/logs
mkdir -p ${OUTPUT}/errs

#FILES_PER_JOB=`expr $3 + 1`
FILES_PER_JOB=$3

TOTAL=`ls $INPUT | grep "lhe" | wc -l`;
JOBS=`expr $TOTAL / $FILES_PER_JOB + 1`;

EVENTS_PER_JOB=$4

echo "Submitting sample : $INPUT"
echo "Output will be stored to : $OUTPUT"
echo "LogFiles will be stored to : $OUTPUT/logs"
echo "total-- $JOBS -- jobs will be submitted to the HPC"

N=0
EVENTNUM=0
JOB=1
FILENAME=SMS-TChiWZ_mChi20-100to500_mChi10_0to200-MG_START53_V19_FSIM_AODSIM_UF
DIRNAME=SMS-TChiWZ_mChi20-100to500_mChi10_0to200-MG
for f in `ls -1d $INPUT/*lhe`; do
N=`expr $N + 1`

if [ -n "$FILES" ]; then FILES="$FILES, file:$f" ; else FILES="file:$f" ;fi;
if [ `expr $N % $FILES_PER_JOB` == 0 ] ; then

SKIPNOW=0

while [ $SKIPNOW -lt 150000 ]; do


mkdir dir_${JOB}
cd dir_${JOB}


echo "#!/bin/bash
#
#PBS -r n
##Job settings
#PBS -N PYTHIA
#PBS -o HPCLogs/Job_${JOB}.txt
#PBS -e HPCErrors/Job_${JOB}.txt

##Job Configuration
##Job Resources
#PBS -l walltime=00:06:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=3gb
# initialize environment for worker
BASE=/lfs/scratch/mrcarver/CMSSW_7_2_1_patch1/src/dir_${JOB}
cd \$BASE
mkdir HPCLogs
mkdir HPCErrors

export SCRAM_ARCH=slc6_amd64_gcc481/
export OSG_APP=/osg/app
export VO_CMS_SW_DIR=\${OSG_APP}/cmssoft/cms
export CMS_PATH=\${VO_CMS_SW_DIR}
. ${CMS_PATH}/cmsset_default.sh;
eval \`scramv1 runtime -sh\`


cd \$BASE
export LD_LIBRARY_PATH=\"/osg/app/cmssoft/cms/osg/app/GLITE/GLITE_3_2_7-0/d-cache/dcap/lib:$LD_LIBRARY_PATH\";

# enter working area
mkdir \${TMPDIR}/errs
mkdir \${TMPDIR}/logs

cmsRun ../${cfg} skipEventsIn=${SKIPNOW} eventNumber=${EVENTNUM} input=\"$FILES\" output=\${TMPDIR}/${FILENAME}_${JOB}.root > ${OUTPUT}/logs/Job_${JOB}.txt 2>  ${OUTPUT}/errs/Job_${JOB}.txt

(
export PATH=/sbin:/usr/sbin:/bin:/usr/bin:.
export LD_LIBRARY_PATH=/usr/lib:/usr/lib64
export X509_USER_PROXY=/scratch/osg/lesya/.cmssoft.proxy
lcg-cp -b -n 1 --vo cms -D srmv2 -T srmv2 -v file:\${TMPDIR}/${FILENAME}_${JOB}.root  srm://srm.ihepa.ufl.edu:8443/srm/v2/server?SFN=/cms/data/store/user/lshchuts/${DIRNAME}/${FILENAME}_${JOB}.root >> ${OUTPUT}/logs/Job_${JOB}.txt 2>> ${OUTPUT}/errs/Job_${JOB}.txt
)

line1=\`cksum \${TMPDIR}/${FILENAME}_${JOB}.root\`
line2=\`cksum /cms/data/store/user/lshchuts/${DIRNAME}/${FILENAME}_${JOB}.root\`

num1=\`echo \$line1 | cut -f 1 -d ' '\`
num2=\`echo \$line2 | cut -f 1 -d ' '\`

if [[ \$num1 == \$num2 ]]
then
echo \$line2 >> ../tchiwz_good.txt
rm \${TMPDIR}/${FILENAME}_${JOB}.root
else
cp \${TMPDIR}/${FILENAME}_${JOB}.root ${OUTPUT}/
echo \$line1 >> ../tchiwz_failed.txt
fi" > submit.sh ;

if [ $JOB -le 10000 ]; then
qsub submit.sh;
fi;

cd ..

#    echo $FILES
JOB=`expr $JOB + 1`
SKIPNOW=`expr $SKIPNOW + $EVENTS_PER_JOB`
EVENTNUM=`expr $EVENTNUM + $EVENTS_PER_JOB`
done;
FILES=""
N=0;
fi;
done;
