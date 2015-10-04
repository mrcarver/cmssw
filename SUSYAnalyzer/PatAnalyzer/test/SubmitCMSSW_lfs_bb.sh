	#!/bin/bash

	cfg=FakeMuonsYa.py
	count=0
	input=$1

	function usage() {
	    echo "sh SubmitCMSSW.sh <input_directory> <output_directory> <files_per_job>"
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

	TOTAL=`ls $INPUT | grep "root" | wc -l`;
	JOBS=`expr $TOTAL / $FILES_PER_JOB + 1`;


	echo "Submitting sample : $INPUT"
	echo "Output will be stored to : $OUTPUT"
	echo "LogFiles will be stored to : $OUTPUT/logs"
	echo "total-- $JOBS -- jobs will be submitted to the HPC"

	N=0
	JOB=1
	for f in `ls -1d $INPUT/*root`; do
	    N=`expr $N + 1`
	    
	    if [ -n "$FILES" ]; then FILES="$FILES, file:$f" ; else FILES="file:$f" ;fi;
	    if [ `expr $N % $FILES_PER_JOB` == 0 ] ; then  

	echo "#!/bin/bash              
	#
	#PBS -r n 
	##Job settings                                                                                                                                                                       
	#PBS -N ROOT                                                                                                                                                                       
	#PBS -o ${OUTPUT}/HPCLogs_${JOB}.txt
	#PBS -e ${OUTPUT}/HPCErrors_${JOB}.txt

	##Job Configuration                                                                                                                                                                  
	##Job Resources                                                                                                                                                                   
	#PBS -l walltime=00:02:00:00 
	#PBS -l nodes=1:ppn=1                                                                                                                                                                
	#PBS -l pmem=1200mb              
	# initialize environment for worker
	BASE=/scratch/osg/lesya/CMSSW_5_3_11_patch6/src/SUSYAnalyzer/PatAnalyzer/test
	cd \$BASE

	export SCRAM_ARCH=slc5_amd64_gcc462/
	export OSG_APP=/osg/app
	export VO_CMS_SW_DIR=\${OSG_APP}/cmssoft/cms
	export CMS_PATH=\${VO_CMS_SW_DIR}
	. ${CMS_PATH}/cmsset_default.sh;
	eval \`scramv1 runtime -sh\`


	cd \$BASE
	export LD_LIBRARY_PATH=\"/osg/app/cmssoft/cms/osg/app/GLITE/GLITE_3_2_7-0/d-cache/dcap/lib:$LD_LIBRARY_PATH\";

	# enter working area

	cmsRun FakeMuonsYa.py input=\"$FILES\" output=${OUTPUT}/Job_${JOB}.root > ${OUTPUT}/logs/Job_${JOB}.txt 2>  ${OUTPUT}/errs/Job_${JOB}.txt
	" > submit.sh ;


	    qsub submit.sh;

	#    echo $FILES
	    FILES=""
	    N=0;
	    JOB=`expr $JOB + 1`

	    fi;    

	done;

	if [ -n "$FILES" ]; then 

	echo "#!/bin/bash              
	#
	#PBS -r n 
	##Job settings                                                                                                                                                                       
	#PBS -N ROOT                                                                                                                                                                       
	#PBS -o ${OUTPUT}/HPCLogs_${JOB}.txt
	#PBS -e ${OUTPUT}/HPCErrors_${JOB}.txt

	##Job Configuration                                                                                                                                                                  
	##Job Resources                                                                                                                                                                   
	#PBS -l walltime=00:02:00:00 
	#PBS -l nodes=1:ppn=1                                                                                                                                                                
	#PBS -l pmem=1200mb              
	# initialize environment for worker
	BASE=/scratch/osg/lesya/CMSSW_5_3_11_patch6/src/SUSYAnalyzer/PatAnalyzer/test
	cd \$BASE

	export SCRAM_ARCH=slc5_amd64_gcc462/
	export OSG_APP=/osg/app
	export VO_CMS_SW_DIR=\${OSG_APP}/cmssoft/cms
	export CMS_PATH=\${VO_CMS_SW_DIR}
	. ${CMS_PATH}/cmsset_default.sh;
	eval \`scramv1 runtime -sh\`


	cd \$BASE
	export LD_LIBRARY_PATH=\"/osg/app/cmssoft/cms/osg/app/GLITE/GLITE_3_2_7-0/d-cache/dcap/lib:$LD_LIBRARY_PATH\";

	# enter working area

cmsRun FakeMuonsYa.py input=\"$FILES\" output=${OUTPUT}/Job_${JOB}.root > ${OUTPUT}/logs/Job_${JOB}.txt 2>  ${OUTPUT}/errs/Job_${JOB}.txt
" > submit.sh ;


    qsub submit.sh;

#    echo $FILES
    FILES=""
    N=0;
    JOB=`expr $JOB + 1`

fi;    


