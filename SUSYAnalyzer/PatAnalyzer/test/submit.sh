#!/bin/bash              
	#
	#PBS -r n 
	##Job settings                                                                                                                                                                       
	#PBS -N ROOT
	#PBS -m a
	#PBS -M mrcarver@phys.ufl.edu                                                                                                                                                                       
	#PBS -o /cms/data/store/user/t2/users/mrcarver/Run2Tuples/Data/ME_1lep/HPCLogs_15.txt
	#PBS -e /cms/data/store/user/t2/users/mrcarver/Run2Tuples/Data/ME_1lep/HPCErrors_15.txt

	##Job Configuration                                                                                                                                                                  
	##Job Resources                                                                                                                                                                   
	#PBS -l walltime=00:10:00:00 
	#PBS -l nodes=1:ppn=1                                                                                                                                                                
	#PBS -l pmem=5200mb              
	# initialize environment for worker
	BASE=/lfs/scratch/mrcarver/CMSSW_7_4_6/src/SUSYAnalyzer/PatAnalyzer/test
	cd $BASE

	export SCRAM_ARCH=slc6_amd64_gcc481/
	export OSG_APP=/osg/app
	export VO_CMS_SW_DIR=${OSG_APP}/cmssoft/cms
	export CMS_PATH=${VO_CMS_SW_DIR}
	. /osg/app/cmssoft/cms/cmsset_default.sh;
	eval `scramv1 runtime -sh`


	cd $BASE
	export LD_LIBRARY_PATH="/osg/app/cmssoft/cms/osg/app/GLITE/GLITE_3_2_7-0/d-cache/dcap/lib:/lfs/scratch/mrcarver/CMSSW_7_4_6/biglib/slc6_amd64_gcc491:/lfs/scratch/mrcarver/CMSSW_7_4_6/lib/slc6_amd64_gcc491:/lfs/scratch/mrcarver/CMSSW_7_4_6/external/slc6_amd64_gcc491/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_6/biglib/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_6/lib/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_6/external/slc6_amd64_gcc491/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/llvm/3.6/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/lib";

	# enter working area

	cmsRun Run2Data_Ntuplizer.py input="file:/cms/data/store/data/Run2015B/MuonEG/MINIAOD/PromptReco-v1//000/252/116/00000/5C48DCCC-5430-E511-ACB2-02163E01267F.root, file:/cms/data/store/data/Run2015B/MuonEG/MINIAOD/PromptReco-v1//000/252/126/00000/DCFBC64C-1031-E511-A70C-02163E01441A.root" output=/cms/data/store/user/t2/users/mrcarver/Run2Tuples/Data/ME_1lep/Job_15.root > /cms/data/store/user/t2/users/mrcarver/Run2Tuples/Data/ME_1lep/logs/Job_15.txt 2>  /cms/data/store/user/t2/users/mrcarver/Run2Tuples/Data/ME_1lep/errs/Job_15.txt
	
